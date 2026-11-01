import copy
from decimal import Decimal,localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym69_nearby_conditioned_locality as z


def decimal(value):
    value=F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class NearbyConditionedLocalityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=z.run()
        cls.eta=F(1,16000)

    def test_fresh_certificate_matches_saved_hash(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_python_312_is_required(self):
        with patch.object(z.native.previous.sys,'version_info',(3,11,0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_218_upstream_files_and_four_local_inputs_are_pinned(self):
        manifest=json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']),218)
        self.assertEqual(len(manifest['local_inputs']),4)
        self.assertEqual(manifest['base_tree'],'bf1636ba3aa0034bace2b6288f2d1c1c15908e6a')
        self.assertIn('papers/yang-mills-certified-benchmark/YM68_AVERAGED_EXTERIOR_LOCALITY.md',
                      manifest['upstream_sha256'])

    def test_upstream_changes_are_detected(self):
        with patch.object(z.native,'digest',return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()

    def test_rehashed_full_collar_claim_fails_fresh_comparison(self):
        with tempfile.TemporaryDirectory() as directory:
            result,pin=Path(directory)/'result',Path(directory)/'pin'
            wrong=copy.deepcopy(self.cert)
            wrong['model']['arbitrary_growing_near_collar_uniformity']=True
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError):
                z.check(self.cert,result,pin)

    def test_evidence_has_no_binary_floats(self):
        def walk(value):
            self.assertNotIsInstance(value,float)
            if isinstance(value,dict):
                for item in value.values(): walk(item)
            elif isinstance(value,(list,tuple)):
                for item in value: walk(item)
        walk(self.cert)

    def test_zero_coupling_is_independent_for_any_declared_near_size(self):
        for n in (0,4,1000):
            row=z.conditioned_budget(0,1,4,n,2)
            for field in ('uniform_near_expected_L1','uniform_near_expected_TV',
                          'conditional_mutual_information','conditional_relative_Fisher'):
                self.assertEqual(row[field],0)
            self.assertEqual(row['near_density_floor'],1)

    def test_domain_rejects_inexact_couplings_and_invalid_counts(self):
        for eta in (-1,F(1,8000),0.00001,True):
            with self.assertRaises(ValueError):
                z.conditioned_budget(eta,20,4,4,2)
        for args in ((self.eta,0,4,4,2),(self.eta,20,0,4,2),(self.eta,20,4,-1,2),
                     (self.eta,20,4,True,2),(self.eta,20,4,4,0)):
            with self.assertRaises(ValueError):
                z.conditioned_budget(*args)

    def test_unresolved_gap_near_endpoint_is_not_silently_used(self):
        with self.assertRaises(ArithmeticError):
            z.conditioned_budget(F(1,8000)-F(1,10**40),10,4,4,2)

    def test_logarithm_intervals_enclose_independent_decimal_values(self):
        with localcontext() as context:
            context.prec=90
            for value,iv in zip((2,40),z.logarithms()):
                actual=Decimal(value).ln()
                self.assertLessEqual(decimal(iv.lo),actual)
                self.assertGreaterEqual(decimal(iv.hi),actual)

    def test_combined_exponent_bound_dominates_unrounded_scalar_formula(self):
        row=z.conditioned_budget(self.eta,1024,4,4,2)
        joint=row['joint_pointwise']; s=8
        rho=z.pointwise.previous.parameters(self.eta)['rho'].lo
        cell=z.pointwise.parameters(self.eta)
        K=z.pointwise.power_upper(cell['kappa'].hi,s)
        r=joint['radius']
        with localcontext() as context:
            context.prec=90
            pref=decimal(1+F(2*2*s,3)+1/(3*rho*1024))
            raw=decimal(40**s)*pref*(-decimal(rho*1024)).exp()/decimal(r**(3*s))
            raw+=decimal(s*cell['density_second_ratio'].hi*K*r*r)
            self.assertGreaterEqual(decimal(joint['bound']),raw)

    def test_joint_pointwise_density_has_no_invalid_cap_two(self):
        joint=z.joint_pointwise_budget(self.eta,1,100,2)
        self.assertEqual(joint['method'],'density_range')
        self.assertGreater(joint['bound'],2)
        row=z.conditioned_budget(self.eta,1,4,96,2)
        self.assertLess(row['uniform_near_expected_L1'],F(15,100))

    def test_no_near_pins_recovers_or_improves_the_previous_averaged_bound(self):
        for R in (128,512,1024):
            old=z.previous.exterior_budget(self.eta,R,4,2)
            new=z.conditioned_budget(self.eta,R,4,0,2)
            self.assertLessEqual(new['uniform_near_expected_L1'],old['expected_L1'])
            self.assertEqual(new['posterior_normalization_upper'],1)
            self.assertEqual(new['analytic_exponent'],F(1,7))

    def test_near_values_are_uniform_but_near_cardinality_remains_in_the_budget(self):
        row=z.conditioned_budget(self.eta,1024,4,4,2)
        self.assertEqual(row['joint_links'],8)
        self.assertEqual(row['fine_joint_links'],16)
        self.assertEqual(row['analytic_exponent'],F(1,13))
        self.assertEqual(row['posterior_normalization_upper']*row['near_density_floor'],2)
        larger=z.conditioned_budget(self.eta,1024,4,8,2)
        self.assertGreater(larger['uniform_near_expected_L1'],row['uniform_near_expected_L1'])

    def test_Fisher_conversion_uses_target_links_only(self):
        cell=z.pointwise.parameters(self.eta)
        coefficient=24*cell['log_second'].hi+8*cell['score'].hi**2
        for n in (0,4,128):
            row=z.conditioned_budget(self.eta,1024,4,n,2)
            self.assertEqual(row['Fisher_coefficient'],z.pointwise.up(coefficient))

    def test_all_displayed_fixed_near_bounds_are_upward(self):
        row=z.conditioned_budget(self.eta,1024,4,4,2)
        self.assertLess(row['uniform_near_expected_TV'],F(210974,10**12))
        self.assertLess(row['conditional_mutual_information'],F(622141,10**13))
        self.assertLess(row['conditional_relative_Fisher'],F(150994,10**13))
        far=z.conditioned_budget(self.eta,2048,4,4,2)
        self.assertLess(far['uniform_near_expected_TV'],F(537202,10**19))
        values=[F(row['uniform_near_expected_L1']) for row in self.cert['fixed_near_budgets']]
        self.assertTrue(all(a>b for a,b in zip(values,values[1:])))

    def test_posterior_denominator_cancels_exactly_in_independent_finite_law(self):
        for row in self.cert['finite_posterior']['cases']:
            self.assertEqual(sum(map(F,row['posterior_exterior'])),1)
            self.assertEqual(F(row['conditional_L1']),F(row['canceled_identity']))
            self.assertLessEqual(F(row['conditional_L1']),2*F(row['joint_residual_over_near_density']))
        self.assertEqual(F(self.cert['finite_posterior']['cases'][0]['conditional_L1']),F(2,9))

    def test_omitting_normalization_residual_would_fail(self):
        row=self.cert['finite_posterior']['cases'][0]
        self.assertEqual(F(row['joint_residual_over_near_density']),F(1,5))
        self.assertGreater(F(row['conditional_L1']),F(row['joint_residual_over_near_density']))

    def test_unconditional_exterior_average_is_the_wrong_measure(self):
        row=self.cert['finite_posterior']['cases'][0]
        self.assertEqual(F(row['wrong_unconditional_exterior_average']),F(13,60))
        self.assertNotEqual(row['wrong_unconditional_exterior_average'],row['conditional_L1'])

    def test_native_posterior_has_correct_mean_and_positive_floor(self):
        data=self.cert['native_posterior']
        self.assertEqual(F(data['density_floor']),F(3,10))
        self.assertEqual(F(data['correct_conditional_target_mean']),0)
        self.assertEqual(F(data['wrong_unweighted_target_mean']),F(-1,216))
        self.assertEqual(F(data['exact_resolvent_integral']),F(10,9))
        self.assertLessEqual(F(data['resolvent_series_enclosure']['lo']),F(10,9))
        self.assertGreaterEqual(F(data['resolvent_series_enclosure']['hi']),F(10,9))

    def test_exact_native_pin_scores_retain_the_posterior_denominator(self):
        rows=self.cert['native_posterior']['cases']
        lookup={(F(r['near_coordinate']),F(r['remote_coordinate'])):r for r in rows}
        self.assertEqual(F(lookup[1,1]['score_at_target_i']),F(1,32))
        self.assertEqual(F(lookup[1,-1]['score_at_target_i']),F(-1,8))
        for remote in (-1,0,F(3,5),1):
            self.assertEqual(F(lookup[0,remote]['target_amplitude']),0)

    def test_square_root_growth_threshold_certifies_c_one(self):
        critical=z.square_root_growth_threshold(self.eta)
        self.assertGreater(critical.lo,F(348724,100000))
        self.assertLess(critical.hi,F(348725,100000))
        self.assertLess(1,critical.lo)
        with self.assertRaises(ValueError):
            z.square_root_growth_threshold(0)

    def test_growing_near_pins_have_decaying_certified_bounds(self):
        rows=self.cert['growing_near_pins']['cases']
        for row in rows:
            self.assertEqual(row['distance'],row['near_links']**2)
            self.assertEqual(row['target_links'],4)
        values=[F(row['uniform_near_expected_L1']) for row in rows]
        self.assertTrue(all(a>b for a,b in zip(values,values[1:])))
        self.assertLess(values[0],F(179606,10**9))
        self.assertLess(values[-1],F(208713,10**24))

    def test_combining_exponents_avoids_magnifying_the_old_rounding_floor(self):
        row=z.conditioned_budget(self.eta,16384,4,128,2)
        r=row['joint_pointwise']['radius']; s=132
        old=z.pointwise.covariance_upper(self.eta,16384,2*s)
        naive=40**s*old/r**(3*s)
        self.assertGreater(naive,2)
        self.assertLess(row['joint_pointwise']['bound'],F(1,10000))
        self.assertGreater(row['joint_pointwise']['exponent_lower'],0)

    def test_scope_preserves_full_collar_Schur_and_continuum_obligations(self):
        model=self.cert['model']
        for key in ('uniform_over_all_fixed_finite_near_pin_values',
                    'distant_exterior_cardinality_independent',
                    'posterior_null_set_common_to_all_near_values',
                    'sparse_growing_near_sets_controlled'):
            self.assertTrue(model[key])
        for key in ('arbitrary_growing_near_collar_uniformity',
                    'worst_case_distant_exterior_configuration_locality',
                    'full_Schur_memory_kernel_locality','time_zero_Markov_identified',
                    'continuum_R4_QFT'):
            self.assertFalse(model[key])


if __name__=='__main__':
    unittest.main()
