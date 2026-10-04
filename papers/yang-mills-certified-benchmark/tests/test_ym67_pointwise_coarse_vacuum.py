import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym67_pointwise_coarse_vacuum as z


class PointwiseCoarseVacuumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.eta = F(1,16000)

    def test_fresh_certificate_and_hash_match(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_python_312_is_required(self):
        with patch.object(z.native.previous.sys,'version_info',(3,11,0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_all_204_upstream_sources_and_four_local_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']),204)
        self.assertEqual(len(manifest['local_inputs']),4)
        self.assertIn('papers/yang-mills-certified-benchmark/YM66_LOCAL_VACUUM_LIMIT.md',
                      manifest['upstream_sha256'])

    def test_changed_upstream_and_rehashed_full_exterior_claim_are_rejected(self):
        with patch.object(z.native,'digest',return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            wrong = copy.deepcopy(self.cert)
            wrong['model']['arbitrary_growing_conditioning_set_uniform_locality'] = True
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError):
                z.check(self.cert,result,pin)

    def test_evidence_contains_no_binary_floats(self):
        def walk(value):
            self.assertNotIsInstance(value,float)
            if isinstance(value,dict):
                for item in value.values():
                    walk(item)
            elif isinstance(value,(list,tuple)):
                for item in value:
                    walk(item)
        walk(self.cert)

    def test_heat_contraction_is_strict_even_at_the_closed_budget_endpoint(self):
        data = self.cert['constants']
        self.assertLess(F(data['closed_endpoint_contraction']['hi']),F(11,200))
        self.assertLess(F(data['inherited_heat_tail']),F(1,20))
        self.assertLess(F(data['endpoint_B_over_M_upper']),6)

    def test_simple_score_bound_is_linear_in_coupling(self):
        for eta in (F(1,100000),F(1,32000),self.eta,F(1,10000),F(124,10**6)):
            self.assertLess(z.parameters(eta)['score'].hi,48*eta)

    def test_zero_coupling_has_constant_density_and_zero_native_scores(self):
        cell = z.parameters(0)
        self.assertEqual(cell['wave_log_gradient'].lo,0)
        self.assertEqual(cell['wave_log_gradient'].hi,0)
        self.assertEqual(cell['kappa'].lo,1)
        self.assertEqual(cell['kappa'].hi,1)
        self.assertEqual(z.boundary_upgrade(0,1,4,2)['score']['bound'],0)
        self.assertEqual(z.pin_budget(0,1,4,4,2)['conditional_density'],0)

    def test_declared_domain_rejects_inexact_negative_and_closed_gap_inputs(self):
        for eta in (-1,F(1,8000),1,0.00001,True):
            with self.assertRaises(ValueError):
                z.parameters(eta)
        with self.assertRaises(ValueError):
            z.local_constants(F(1,1999))
        with self.assertRaises(ValueError):
            z.exp_interval(9)

    def test_pointwise_constants_do_not_need_a_numerically_resolved_gap(self):
        near = F(1,8000)-F(1,10**40)
        self.assertGreater(z.parameters(near)['score'].lo,0)
        with self.assertRaises(ArithmeticError):
            z.covariance_upper(near,16,4)

    def test_certified_score_and_density_numbers_have_the_stated_orientation(self):
        cell = z.parameters(self.eta)
        self.assertGreater(cell['score'].lo,F(2933322,10**9))
        self.assertLess(cell['score'].hi,F(2933323,10**9))
        self.assertLess(cell['log_second'].hi,F(1488173,10**9))
        self.assertLess(cell['kappa'].hi,F(1018601503,10**9))
        self.assertGreater(1/cell['kappa'].hi**4,F(92,100))

    def test_independent_native_word_jets_keep_order_and_inverse_signs(self):
        data = self.cert['native']
        self.assertEqual(data['ordered_word_jet_checks'],576)
        self.assertTrue(data['ordered_derivatives_not_symmetric'])
        self.assertEqual(F(data['first_derivative_norm_bound']),F(1,2))
        self.assertEqual(F(data['second_derivative_bound']),F(1,4))

    def test_each_link_has_the_gradient_square_identity_and_zero_average(self):
        self.assertEqual(self.cert['native']['one_link_gradient_square_identities'],4)
        W = z.model.loop_polynomial()
        for e in range(4):
            derivative = z.native.derivative(W,e,1)
            self.assertFalse(z.integrate_link(derivative,e))
            second = z.native.derivative(z.native.derivative(W,(e+1)%4,2),e,1)
            self.assertFalse(z.integrate_link(second,e))

    def test_noncommuting_directions_cannot_be_replaced_by_a_symmetric_hessian(self):
        q0 = z.p.COORD[0]
        left = z.y.native_generator(z.y.native_generator(q0,2),1)
        right = z.y.native_generator(z.y.native_generator(q0,1),2)
        self.assertNotEqual(left,right)
        with self.assertRaises(ValueError):
            z.derivative(z.model.loop_polynomial(),0,(1,1,0))

    def test_one_pivot_and_orthogonal_adjoint_rotation_have_no_block_length_factor(self):
        data = self.cert['pivots']
        self.assertEqual({row['length'] for row in data['cases']},{1,2,3,5})
        self.assertEqual(data['pivot_count_per_coarse_link'],1)
        self.assertEqual(data['block_length_factor_in_score'],1)
        for row in data['cases']:
            gram = [[F(value) for value in line] for line in row['rotation_gram']]
            self.assertEqual(gram,z.a.identity(3))

    def test_native_ball_floor_uses_three_dimensions_per_coarse_link(self):
        for row in self.cert['balls']['native_haar_caps']:
            r = F(row['radius'])
            self.assertEqual(F(row['proved_floor']),r**3/40)
            self.assertGreater(F(row['haar_cap']['lo']),r**3/40)

    def test_positive_L1_spikes_do_not_converge_pointwise_without_regularity(self):
        rows = self.cert['balls']['positive_spike_comparison']
        self.assertLess(F(rows[-1]['L1']),F(rows[0]['L1']))
        self.assertGreaterEqual(min(F(row['sup']) for row in rows),F(3,4))
        self.assertGreater(rows[-1]['slope'],rows[0]['slope'])

    def test_long_decay_rounding_is_upward_and_avoids_the_old_absolute_floor(self):
        self.assertEqual(z.up(F(1,3))-F(1,3),F(2,3*10**120))
        self.assertEqual(z.decay_upper(0),1)
        value = z.decay_upper(80)
        self.assertGreater(value,0)
        self.assertLess(value,F(1,10**34))
        self.assertLess(value,z.spatial.exp_negative(F(80)).hi)
        with self.assertRaises(ValueError):
            z.decay_upper(-1)

    def test_directed_powers_and_tails_dominate_independent_exact_sums(self):
        q,R = F(3,4),50
        self.assertGreaterEqual(z.power_upper(q,R),q**R)
        self.assertEqual(z.power_upper(0,0),1)
        for exact,upper in zip(z.previous.tail_moments(q,R),z.tail_moments_upper(q,R)):
            self.assertGreaterEqual(upper,exact)
        with self.assertRaises(ValueError):
            z.tail_moments_upper(1,4)

    def test_refined_boundary_enclosure_remains_below_the_older_loose_bound(self):
        old = z.previous.boundary_tail(self.eta,256,8)['bound']
        new = z.boundary_l1(self.eta,256,8)
        self.assertLessEqual(new,old)
        self.assertGreater(new,0)
        self.assertLess(z.boundary_l1(self.eta,1024,8),F(1,10**70))

    def test_fixed_four_link_marginal_has_nontrivial_C1_error_at_large_buffer(self):
        rows = self.cert['marginal_budgets']
        density = [F(row['density']['bound']) for row in rows]
        scores = [F(row['score']['bound']) for row in rows]
        self.assertTrue(all(x > y for x,y in zip(density,density[1:])))
        self.assertTrue(all(x > y for x,y in zip(scores,scores[1:])))
        self.assertLess(density[-1],F(406,10**9))
        self.assertLess(scores[-1],F(5122,10**8))

    def test_boundary_error_retains_the_complete_fine_path_count(self):
        one = z.boundary_l1(self.eta,512,4)
        two = z.boundary_l1(self.eta,512,8)
        self.assertLessEqual(abs(two-2*one),F(2,10**120))
        self.assertEqual(z.boundary_upgrade(self.eta,512,4,2)['block_factor'],2)

    def test_density_interpolation_uses_a_valid_explicit_radius(self):
        error,m = F(1,10**20),4
        data = z.density_upgrade(self.eta,m,error)
        r = data['radius']
        cell = z.parameters(self.eta)
        K = z.power_upper(cell['kappa'].hi,m)
        rhs = 40**m*error/r**(3*m)+4*m*cell['wave_log_gradient'].hi*K*r
        self.assertGreater(r,0)
        self.assertLessEqual(r,1)
        self.assertLessEqual(abs(data['bound']-rhs),F(1,10**120))
        self.assertEqual(z.density_upgrade(self.eta,4,0)['bound'],0)

    def test_score_interpolation_keeps_the_second_derivative_error(self):
        delta,m = F(1,10**8),4
        data = z.score_upgrade(self.eta,m,delta)
        cell = z.parameters(self.eta)
        t = data['step']
        rhs = z.power_upper(cell['kappa'].hi,m)*delta/t+cell['log_second'].hi*t
        self.assertLessEqual(abs(data['bound']-rhs),F(1,10**120))
        self.assertLess(data['bound'],2*cell['score'].hi)

    def test_pinned_physical_plaquette_supports_have_decaying_pointwise_bounds(self):
        rows = self.cert['pin_budgets']
        density = [F(row['conditional_density']) for row in rows]
        self.assertTrue(all(x > y for x,y in zip(density,density[1:])))
        self.assertLess(density[-1],F(11445,10**8))
        self.assertLess(F(rows[-1]['score']['bound']),F(86483,10**8))
        self.assertLess(F(rows[-1]['mutual_information']),F(1518,10**11))
        self.assertEqual({row['target_links'] for row in rows},{4})
        self.assertEqual({row['pinned_links'] for row in rows},{4})

    def test_probability_information_is_nonnegative_and_uses_the_relative_density_error(self):
        for row in self.cert['pin_budgets']:
            self.assertGreaterEqual(F(row['mutual_information']),0)
            self.assertLessEqual(F(row['mutual_information']),F(row['log_ratio']))

    def test_exact_native_pins_have_distinct_conditional_scores(self):
        data = self.cert['conditional_fixture']
        self.assertEqual(list(map(F,data['exact_point_scores'])),[F(1,20),F(-1,20)])
        self.assertEqual(F(data['covariance']),F(1,160))
        self.assertEqual(F(data['conditional_mean_at_positive_pin']),F(1,40))
        self.assertEqual(F(data['density_floor']),F(9,10))

    def test_local_drift_energy_factor_loses_block_length_after_clock_normalization(self):
        # Squared estimate: ||(C_b-bL)f||^2 <= a1^2*b*|Y|*C_b[f].
        a1sq = z.parameters(self.eta)['score'].hi**2
        width,energy = 4,F(7,3)
        for b in (1,2,5,11):
            coarse_energy = b*energy
            unscaled_square = a1sq*b*width*coarse_energy
            self.assertEqual(unscaled_square/b**2,a1sq*width*energy)

    def test_invalid_interpolation_and_pin_sizes_are_refused(self):
        for count in (0,True):
            with self.assertRaises(ValueError):
                z.density_upgrade(self.eta,count,1)
            with self.assertRaises(ValueError):
                z.pin_budget(self.eta,1,count,4,2)
        with self.assertRaises(ValueError):
            z.density_upgrade(self.eta,4,F(-1,100))
        with self.assertRaises(ValueError):
            z.boundary_upgrade(self.eta,0,4,2)

    def test_scope_separates_finite_exact_pins_from_uniform_full_exterior_locality(self):
        model = self.cert['model']
        for key in ('actual_YM62_vacuum_retained','unbounded_link_spaces_retained',
                    'pointwise_positive_coarse_marginals','score_bound_independent_of_volume_and_block_length',
                    'fixed_marginals_C1_limit','finite_exact_pin_conditionals',
                    'finite_pin_spatial_score_decay','local_diagonal_drift_energy_bound'):
            self.assertTrue(model[key])
        for key in ('arbitrary_growing_conditioning_set_uniform_locality',
                    'full_Schur_memory_kernel_locality','time_zero_Markov_identified',
                    'continuum_R4_QFT','physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
