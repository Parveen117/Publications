import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym52_energy_observer as z


class NativeEnergyObserverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert=z.run()

    def test_fresh_certificate_and_pin(self):
        self.assertEqual(z.check(self.cert),z.PIN.read_text().strip())

    def test_product_rule_and_both_energy_derivatives(self):
        self.assertEqual(self.cert['energy']['leibniz_and_integration_by_parts'],441)
        self.assertEqual(self.cert['energy']['norm_and_energy_derivatives'],63)
        self.assertEqual(z.energy(z.p.COORD[0],z.q.I3),F(3,16))
        self.assertNotEqual(-z.energy(z.p.COORD[0],z.q.I3),F(-3,8))

    def test_record_noise_is_independently_replayed(self):
        d=self.cert['discarded_records']
        self.assertEqual(d['independent_word_variance_checks'],12)
        self.assertEqual(d['largest_literal_word_bank'],64)
        self.assertTrue(any(F(x)>0 for x in d['integrated_noise']))

    def test_cofactor_does_not_need_inverse(self):
        C=z.a.diag((0,1,2))
        self.assertEqual(z.cofactor(C),z.a.diag((2,0,0)))
        self.assertEqual(z.cofactor(z.a.diag((0,0,3))),z.a.diag((0,0,0)))
        self.assertEqual(self.cert['bracket_tensor']['cofactor_factorizations'],5)

    def test_noncommuting_rank_two_relaxes(self):
        self.assertEqual(z.rates((0,F(3,2),F(3,2))),(F(3,4),F(3,2),F(3,4)))
        self.assertEqual(z.rates((0,0,3)),(0,0,F(3,4)))
        self.assertEqual(z.q.lap(z.EVEN,z.a.diag((0,0,3))),{})
        self.assertGreater(z.variance(z.EVEN),0)

    def test_exact_spin_energy_has_all_directions_and_parity(self):
        d=self.cert['spin_blocks']
        self.assertEqual(d['finite_spin_blocks'],12)
        self.assertEqual(d['gap_matrix_inequalities'],108)
        self.assertEqual(d['half_integer_axis_floor_checks'],18)

    def test_independent_polynomial_matrix_energy_bounds(self):
        self.assertEqual(self.cert['polynomial_gaps']['all_source_polynomial_matrix_checks'],40)

    def test_YM51_slow_bound_is_sharp_only_in_its_range(self):
        for eta in (F(1,1000),F(1,4),F(3,8),F(1)):
            self.assertEqual(z.rates((eta,eta,3-2*eta))[0],min(F(3,4),2*eta))
        self.assertNotEqual(z.rates((1,1,1))[0],F(2))

    def test_projection_changes_the_relaxation_target(self):
        f=z.p.add(z.p.COORD[0],z.EVEN)
        self.assertEqual(z.parity(f),z.EVEN)
        self.assertEqual(z.parity(f,False),z.p.COORD[0])
        self.assertEqual(z.parity(z.p.COORD[0]),{})
        self.assertEqual(z.rates((1,1,1)),(F(3,4),2,F(3,4)))
        self.assertEqual(self.cert['observer_selection']['parity_energy_checks'],9)

    def test_full_objective_does_not_select_isotropy(self):
        self.assertEqual(z.rates((0,F(3,2),F(3,2)))[0],z.rates((1,1,1))[0])
        self.assertLess(z.rates((0,F(3,2),F(3,2)))[1],z.rates((1,1,1))[1])
        self.assertTrue(self.cert['observer_selection']['full_gap_optimization_does_not_select_isotropy'])

    def test_even_selection_stability_constant_cannot_be_halved(self):
        eps=F(1,4);vals=(1-2*eps,1+eps,1+eps)
        self.assertEqual(2-z.rates(vals)[1],eps)
        self.assertEqual(max(abs(v-1) for v in vals),2*eps)
        self.assertGreater(max(abs(v-1) for v in vals),eps)

    def test_entropy_identity_and_rational_remainders(self):
        d=self.cert['entropy']
        self.assertEqual(d['entropy_differentiation_coefficients'],42)
        self.assertEqual(d['rational_entropy_production_enclosures'],6)
        for row in d['cases']:
            h0,h1=map(F,row['entropy_interval'])
            i0,i1=map(F,row['production_interval'])
            self.assertGreater(h0,0);self.assertGreaterEqual(h1,h0)
            self.assertGreaterEqual(i1,i0)
            self.assertGreaterEqual(i0,F(row['bounded_density_decay_lower'])*h1)

    def test_exact_domains_and_positivity_bounds_are_required(self):
        for values in ((1,0,2),(-1,1,2),(1,2),(0.0,1,2),(True,1,2)):
            with self.assertRaises(ValueError):z.rates(values)
        for eps in (0,1,-1,F(3,2)):
            with self.assertRaises(ValueError):z.entropy_series(z.p.COORD[0],z.q.I3,eps)
        with self.assertRaises(ValueError):
            z.entropy_series(z.p.scale(z.p.COORD[0],2),z.q.I3,F(3,4))
        self.assertEqual(z.quadratic_uniform_bound(z.p.scale(z.EVEN,2)),1)
        with self.assertRaises(ValueError):
            z.entropy_series(z.y.poly_power(z.p.COORD[0],3),z.q.I3,F(1,4))
        with self.assertRaises(ValueError):z.entropy_series(z.p.ONE,z.q.I3,F(1,4))
        with self.assertRaises(ValueError):z.spin_squares(-1)
        with self.assertRaises(ValueError):z.cofactor(z.a.diag((-1,1,2)))

    def test_clock_scaling_keeps_bracket_units_distinct(self):
        self.assertEqual(self.cert['clock']['linear_clock_and_quadratic_bracket_scalings'],21)
        C=z.a.diag((1,2,3))
        self.assertNotEqual(z.cofactor(z.a.scale(C,2)),z.a.scale(z.cofactor(C),2))

    def test_verdict_and_hash_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            result=Path(tmp)/'result.json';pin=Path(tmp)/'pin'
            changed=copy.deepcopy(self.cert)
            changed['polynomial_gaps']['rank_one_centered_rate']='3/4'
            result.write_text(json.dumps(changed))
            pin.write_text(z.y.canonical_sha(changed))
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)
            result.write_text(json.dumps(self.cert));pin.write_text('0'*64)
            with self.assertRaises(ValueError):z.check(self.cert,result,pin)

    def test_upstream_hash_mismatch_is_rejected(self):
        original=z.digest
        try:
            z.digest=lambda path:'0'*64
            with self.assertRaises(ValueError):z.source_checks()
        finally:z.digest=original


if __name__=='__main__':unittest.main()
