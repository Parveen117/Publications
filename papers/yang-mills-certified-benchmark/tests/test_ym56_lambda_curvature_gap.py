import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym56_lambda_curvature_gap as z


class LambdaGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_is_pinned(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_hessian_integrability_requires_scalar_strain(self):
        for lam in z.SAMPLES:
            d = z.family(lam)
            dx,dy = d['derivatives']
            self.assertEqual(dx[1][1],dy[0][1])
            self.assertEqual(dx[0][1],dy[0][0])
            self.assertEqual(d['shape'][0],z.native.K)
        self.assertNotEqual(z.native.K[1][1],z.a.scale(z.native.L,2)[0][1])

    def test_positive_chart_and_orientation_match_tangent_convention(self):
        self.assertEqual(self.cert['hessian']['positivity_corners'],48)
        for lam in z.SAMPLES:
            d = z.family(lam)
            self.assertEqual(d['F'][0][1],-lam/2)
            self.assertEqual(d['marker'],lam/2)

    def test_four_labels_keep_symmetry_at_nonzero_gap(self):
        self.assertGreater(z.rates(1)[0],0)
        self.assertEqual(self.cert['symmetry']['literal_four_label_moments'],9)
        self.assertEqual(self.cert['symmetry']['polynomial_symmetry_checks'],306)

    def test_squared_protocol_forgets_curvature_sign(self):
        for lam in (F(1,10),F(1),F(2)):
            left,right = z.family(lam),z.family(-lam)
            self.assertEqual(left['C'],right['C'])
            self.assertEqual(left['F'],z.a.scale(right['F'],-1))
            self.assertNotEqual(left['F'],right['F'])

    def test_closed_form_has_independent_sharp_polynomial_modes(self):
        self.assertEqual(self.cert['gap']['sharp_quadratic_modes'],18)
        self.assertEqual(self.cert['gap']['exact_spin_inequalities'],108)
        self.assertEqual(z.rates(F(1,10))[0],F(1,200))
        self.assertEqual(z.rates(1),(F(1,4),F(1,2),F(1,4)))
        self.assertEqual(z.rates(2),(F(1,2),F(1,2),F(5,8)))

    def test_zero_parameter_has_extra_kernel_not_zero_operator(self):
        C = z.family(0)['C']
        self.assertFalse(z.q.lap(z.witness(1),C))
        self.assertGreater(z.energy.variance(z.witness(1)),0)
        self.assertTrue(z.q.lap(z.p.COORD[0],C))
        self.assertEqual(z.rates(0),(F(0),F(0),F(1,8)))

    def test_asymptotic_quadratic_rate_is_not_valid_everywhere(self):
        for n in (2,3,10,100):
            self.assertEqual(z.rates(F(1,n))[0],F(1,2*n*n))
        self.assertNotEqual(z.rates(1)[0],F(1,2))

    def test_curvature_normalization_handles_rotated_and_rank_one_sources(self):
        self.assertEqual(self.cert['general_shape']['rotated_shape_families'],4)
        for lam in (F(3,5),F(1),F(5,3)):
            full = z.rates(lam)[0]
            self.assertEqual(16*full/(2*(1+lam*lam)),1)
        self.assertLess(16*z.rates(F(1,2))[0]/F(5,2),1)

    def test_same_dimensionless_shape_does_not_fix_clock(self):
        lam = F(2)
        chi = 2*lam/(1+lam*lam)
        inverse_chi = 2/lam/(1+1/(lam*lam))
        self.assertEqual(chi,inverse_chi)
        self.assertNotEqual(z.rates(lam)[0],z.rates(1/lam)[0])
        self.assertEqual(z.rates(lam)[0],lam*lam*z.rates(1/lam)[0])

    def test_direct_profile_floor_passes_where_generic_bound_abstains(self):
        c = self.cert['chain']
        self.assertEqual(F(c['beta']),F(1,8))
        self.assertEqual(F(c['curvature_only_beta']),F(1,40))
        self.assertGreater(F(c['theta']),F(c['curvature_only_beta'])/2400)
        self.assertLess(F(c['theta']),F(c['beta'])/2400)
        self.assertGreater(F(c['chain_gamma_lower']),0)
        self.assertEqual(c['profile_tensor_checks'],8)

    def test_profile_gate_refuses_missing_floor_and_strict_endpoint(self):
        for ell,bound in ((0,1),(-1,1),(2,1),(1,0)):
            with self.assertRaises(ValueError):z.profile_bounds(ell,bound)
        b = z.profile_bounds(F(1,2),2)
        with self.assertRaises(ValueError):
            z.joint.parameters(b['beta'],b['theta_ceiling'])
        with self.assertRaises(ValueError):
            z.joint.trace_profile([z.family(F(1,10))['C']],b['beta'],b['trace_max'])

    def test_each_site_gapped_does_not_imply_uniform_free_profile(self):
        rows = self.cert['chain']['nonuniform_free_profile']
        rates = [F(row['quotient']) for row in rows]
        self.assertTrue(all(r>0 for r in rates))
        self.assertEqual(rates,sorted(rates,reverse=True))
        self.assertEqual(rates[-1],F(1,20000))

    def test_exact_input_contract_rejects_float_bool_and_bad_axis(self):
        for value in (True,0.5,'1/2'):
            with self.assertRaises(ValueError):z.family(value)
            with self.assertRaises(ValueError):z.rates(value)
        for axis in (0,4,True):
            with self.assertRaises(ValueError):z.witness(axis)

    def test_python_312_is_required(self):
        z.runtime_check()
        with patch.object(z.sys,'version_info',(3,11,0)):
            with self.assertRaises(RuntimeError):z.runtime_check()

    def test_source_drift_and_certificate_tampering_are_rejected(self):
        with patch.object(z,'digest',return_value='0'*64):
            with self.assertRaises(ValueError):z.source_checks()
        with tempfile.TemporaryDirectory() as tmp:
            result,pin = Path(tmp)/'result',Path(tmp)/'pin'
            altered = copy.deepcopy(self.cert)
            altered['claim_status'] = 'CLAY_SOLVED'
            result.write_text(json.dumps(altered))
            pin.write_text(z.native.chain.old.canonical_sha(altered))
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)
            result.write_text(json.dumps(self.cert));pin.write_text('0'*64)
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)

    def test_free_chain_and_physical_scopes_remain_distinct(self):
        self.assertFalse(self.cert['chain']['physical_mass_identified'])
        self.assertTrue(self.cert['chain']['fixed_profile_required'])
        self.assertIn('4D_CLAY_OPEN',self.cert['claim_status'])
        self.assertIn('written',self.cert['evidence_scope']['general_proof'])


if __name__ == '__main__':
    unittest.main()
