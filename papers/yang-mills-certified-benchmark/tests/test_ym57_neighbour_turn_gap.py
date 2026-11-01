import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym57_neighbour_turn_gap as z


class NeighbourGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.cert=z.run()

    def test_fresh_exact_certificate(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_neighbour_echo_cancels_other_site_exactly(self):
        t=(F(3,5),F(0),F(4,5),F(0))
        self.assertEqual(z.echo(t),(z.y.qmul(t,t),z.y.ONE))
        self.assertNotEqual(z.echo(t,z.y.ONE)[1],z.y.ONE)

    def test_existing_scalar_potential_preserves_missing_direction(self):
        self.assertEqual(self.cert['obstruction']['potential_commutation_checks'],15)
        w=z.lift(z.previous.witness(1),0,2)
        self.assertFalse(z.operator(w,2,0))
        self.assertTrue(z.bond(w,0))

    def test_pointwise_witness_derivative_has_correct_factor(self):
        dw=z.y.native_generator(z.previous.witness(1),2)
        self.assertEqual(len(dw),2)
        self.assertTrue(all(abs(c)==1 for c in dw.values()))
        self.assertTrue(self.cert['obstruction']['weighted_variance_must_not_be_replaced'])

    def test_correlated_square_is_not_two_independent_site_squares(self):
        f=z.p.mul(z.lift(z.p.COORD[0],0,2),z.lift(z.p.COORD[0],1,2))
        mixed=z.p.add(z.operator(f,2,1),z.p.scale(z.operator(f,2,0),-1))
        separate=z.p.scale(z.p.add(z.derivative(z.derivative(f,0,2),0,2),
                                   z.derivative(z.derivative(f,1,2),1,2)),F(-1,2))
        self.assertNotEqual(mixed,separate)

    def test_independent_native_square_factorizations_and_gap(self):
        self.assertEqual(self.cert['matrices']['rational_gap_inequalities'],24)
        self.assertEqual(self.cert['matrices']['independent_square_factorizations'],24)
        self.assertEqual(self.cert['matrices']['max_block_dimension'],40)

    def test_actual_signed_records_reproduce_cross_site_covariance(self):
        self.assertEqual(self.cert['records']['exact_signed_record_second_moments'],9)
        self.assertTrue(self.cert['records']['singular_covariance_is_not_a_gaplessness_test'])

    def test_all_width_incidence_budget_is_uniform(self):
        for row in self.cert['graph']['rows']:
            self.assertLessEqual(row['max_indegree'],2)
            self.assertLessEqual(row['max_edge_uses'],2)
            self.assertEqual(z.lower_bound(row['width'],1),F(1,630))

    def test_small_edge_rate_is_not_hidden_by_fixed_constant(self):
        self.assertEqual(z.lower_bound(3,F(1,100)),F(1,42000))
        self.assertEqual(z.lower_bound(3,10),F(1,630))

    def test_zero_bond_or_isolated_site_has_stationary_centered_mode(self):
        for n,k in ((1,1),(2,0),(3,0)):
            w=z.lift(z.previous.witness(1),0,n)
            self.assertFalse(z.operator(w,n,k))
            self.assertEqual(z.variance(w,n),F(1,12))

    def test_local_lambda_addition_preserves_lower_form(self):
        w=z.lift(z.previous.witness(1),0,2)
        for lam in (F(-2),F(0),F(1,10)):
            ratio=z.phi(z.p.mul(w,z.operator(w,2,1,lam)),2)/z.variance(w,2)
            self.assertEqual(ratio,(1+lam*lam)/2)

    def test_global_half_turn_symmetry_survives(self):
        self.assertEqual(self.cert['symmetry_clock']['local_symmetry_checks'],11)

    def test_whole_chain_fixed_rate_clock_loses_uniformity(self):
        rows=self.cert['symmetry_clock']['endpoint_witnesses']
        for row in rows:
            n,k=row['width'],F(row['k'])
            self.assertEqual(F(row['normalized_endpoint_rayleigh']),k/(2*(n+k*(n-1))))
        upper=list(map(F,self.cert['symmetry_clock']['fixed_rate_clock_upper_bounds']))
        self.assertLess(upper[-1],F(1,630))

    def test_invalid_width_rate_and_inexact_inputs_are_refused(self):
        for n,k in ((1,1),(0,1),(True,1),(2,0),(2,-1),(2,0.1),(2,True)):
            with self.assertRaises(ValueError):z.lower_bound(n,k)
        with self.assertRaises(ValueError):z.lift(z.p.ONE,2,2)

    def test_pin_and_source_drift_are_refused(self):
        with patch.object(z,'digest',return_value='0'*64):
            with self.assertRaises(ValueError):z.source_checks()
        with tempfile.TemporaryDirectory() as tmp:
            result,pin=Path(tmp)/'result',Path(tmp)/'pin'
            bad=copy.deepcopy(self.cert);bad['claim_status']='CLAY_SOLVED'
            result.write_text(json.dumps(bad));pin.write_text(z.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)

    def test_model_and_proof_scope_remain_explicit(self):
        self.assertFalse(self.cert['symmetry_clock']['scalar_interaction_imported'])
        self.assertFalse(self.cert['symmetry_clock']['infinite_volume_constructed'])
        self.assertIn('4D_CLAY_OPEN',self.cert['claim_status'])
        self.assertIn('all-content',self.cert['proof_scope'])


if __name__=='__main__':unittest.main()
