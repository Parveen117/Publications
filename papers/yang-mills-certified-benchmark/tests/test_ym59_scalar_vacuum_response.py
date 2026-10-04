import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym59_scalar_vacuum_response as z


class ScalarVacuumResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_exact_certificate_is_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_energy_pays_for_each_excited_site(self):
        self.assertEqual(self.cert['coercivity']['excited_site_form_inequalities'], 15)
        f = z.p.mul(z.native.lift(z.p.COORD[0], 0, 2), z.native.lift(z.p.COORD[0], 1, 2))
        for site in range(2):
            self.assertFalse(z.site_average(f, site, 2))
        self.assertEqual(sum(z.pairing(f, f, 2) for _ in range(2)), F(1, 8))

    def test_site_averaging_commutes_with_correlated_heat(self):
        self.assertEqual(self.cert['coercivity']['projection_commutations'], 132)
        f = z.p.mul(z.overlap(0, 1, 3), z.overlap(1, 2, 3))
        self.assertEqual(z.site_average(f, 1, 3), z.p.scale(z.overlap(0, 2, 3), F(1, 4)))

    def test_open_boundary_bond_rates_are_not_bulk_rates(self):
        self.assertEqual(z.omega(1), 1)
        self.assertEqual(z.omega(1, 1, 2), F(7, 8))
        self.assertEqual(z.omega(1, 1, 1), F(3, 4))
        self.assertEqual(self.cert['bonds']['local_bond_eigenvalue_checks'], 36)

    def test_mixed_native_derivative_has_the_required_positive_sign(self):
        v = z.overlap(0, 1, 2)
        mixed = z.native.derivative(z.native.derivative(v, 0, 2), 1, 2)
        self.assertEqual(mixed, z.p.scale(v, F(1, 4)))
        self.assertEqual(z.local_operator(v, (0, 1), 1), v)

    def test_poisson_inverse_is_complete_for_its_source(self):
        for distance, expected in enumerate((F(1, 30), F(31, 705), F(2, 63), F(1, 32))):
            data = z.fixture(distance)
            residual = z.p.add(z.local_operator(data['inverse'], data['support'], 1),
                               z.p.scale(data['source'], -1))
            self.assertEqual(z.pairing(residual, residual, len(data['support'])), 0)
            self.assertEqual(data['quadratic'], expected)
        self.assertEqual([r['dimension'] for r in self.cert['poisson']['exact_local_inverses']], [3, 4, 2, 1])

    def test_noncentered_poisson_source_is_refused(self):
        with self.assertRaises(ValueError):
            z.inverse_source(z.constant(1, 2), (0, 1), 1)
        self.assertEqual(z.inverse_source({}, (0, 1), 1)['dimension'], 0)

    def test_logarithmic_and_rayleigh_routes_agree(self):
        self.assertEqual(F(self.cert['poisson']['logarithmic_c4_single_bond']), F(7, 240))
        self.assertEqual(z.energy_coefficients(2, boundary='bulk'), (F(-1, 4), F(7, 240)))
        self.assertEqual(self.cert['poisson']['logarithmic_residual_identities'], 3)

    def test_square_factorization_retains_the_original_scalar_operator(self):
        support = (0, 1)
        v = z.overlap(0, 1, 2)
        Fp = z.p.add(z.p.scale(v, F(2, 7)), z.native.lift(z.p.COORD[1], 0, 2))
        theta = F(-1, 9)
        residual = z.p.add(z.local_operator(Fp, support, 1),
                          z.p.scale(z.gamma(Fp, Fp, support, 1), -1), z.p.scale(v, -theta))
        for f in (z.constant(1, 2), v, z.p.mul(v, v)):
            squares = []
            for weight, word in z.prior.directions(support, 1):
                D = lambda g: z.prior.apply_direction(g, word)
                dF = D(Fp)
                first = z.p.add(D(f), z.p.scale(z.p.mul(dF, f), -1))
                adjoint = z.p.add(z.p.scale(D(first), -1), z.p.scale(z.p.mul(dF, first), -1))
                squares.append(z.p.scale(adjoint, weight/2))
            actual = z.p.add(*squares, z.p.mul(residual, f))
            expected = z.p.add(z.local_operator(f, support, 1), z.p.scale(z.p.mul(v, f), -theta))
            self.assertEqual(actual, expected)

    def test_kinetic_bridge_has_a_nonzero_connected_contribution(self):
        bridge = F(self.cert['fourth_order']['connected_terms'][2])
        self.assertEqual(bridge, F(-1, 504))
        self.assertNotEqual(z.fixture(2)['quadratic'], z.fixture(3)['quadratic'])
        self.assertNotEqual(self.cert['fourth_order']['deleting_kinetic_bridge_gives'],
                            self.cert['fourth_order']['bulk_quartic'])

    def test_disconnected_pairs_cancel_with_unequal_rates(self):
        self.assertEqual(self.cert['fourth_order']['unequal_rate_disconnected_cancellations'], 9)
        self.assertEqual(F(self.cert['fourth_order']['normalization_removed_second_difference']), F(-1, 8))
        u, v = F(7, 8), F(1)
        subtract = (1/u+1/v)**2/(16*(u+v))
        restore = (u+v)/(16*u*u*v*v)
        self.assertEqual(subtract, restore)
        self.assertGreater(restore, 0)

    def test_bulk_quartic_is_sum_of_three_local_contributions(self):
        pieces = tuple(map(F, self.cert['fourth_order']['connected_terms']))
        self.assertEqual(pieces, (F(7, 240), F(-287, 5640), F(-1, 504)))
        self.assertEqual(sum(pieces), F(-1123, 47376))

    def test_assembled_four_site_source_checks_the_full_formula(self):
        self.assertTrue(self.cert['fourth_order']['assembled_full_source_inverse'])
        self.assertEqual(F(self.cert['fourth_order']['assembled_c4']), F(-3851, 236880))

    def test_open_and_marginal_boundaries_have_the_same_bulk_slope(self):
        bulk = F(self.cert['fourth_order']['bulk_quartic'])
        for boundary in ('open', 'bulk'):
            _, five = z.energy_coefficients(5, boundary=boundary)
            _, six = z.energy_coefficients(6, boundary=boundary)
            self.assertEqual(six-five, bulk)
        self.assertNotEqual(z.energy_coefficients(5), z.energy_coefficients(5, boundary='bulk'))

    def test_susceptibility_and_quadratic_normalization(self):
        self.assertEqual(F(self.cert['bonds']['bulk_susceptibility_at_unit_rate']), F(1, 2))
        self.assertEqual(z.energy_coefficients(2)[0], F(-1, 3))
        self.assertEqual(z.energy_coefficients(7)[0], F(-11, 7))

    def test_alternating_antipodal_symmetry_reverses_only_the_coupling(self):
        self.assertEqual(self.cert['bonds']['antipodal_intertwinings'], 12)
        self.assertEqual(self.cert['poisson']['logarithmic_c2'], '-1/4')
        self.assertTrue(self.cert['poisson']['independent_square_product_rule'])

    def test_excitation_incidence_and_young_bounds_are_uniform(self):
        self.assertEqual(self.cert['energy']['excitation_incidence_cases'], 508)
        for row in self.cert['energy']['weak_energy_cells']:
            self.assertEqual(F(row['excitation_upper_per_bond']), F(1, 36))
            self.assertLess(F(row['energy_lower_per_bond']), 0)
        self.assertEqual(z.weak_energy_budget(1, 0)['energy_lower_per_bond'], 0)

    def test_dimer_trial_factorizes_cross_bonds_but_retains_boundary_noise(self):
        n, c, theta = 4, F(1, 10), F(1, 10)
        v0, v1, v2 = [z.overlap(i, i+1, n) for i in range(3)]
        psi = z.p.mul(z.p.add(z.constant(1, n), z.p.scale(v0, c)),
                      z.p.add(z.constant(1, n), z.p.scale(v2, c)))
        self.assertEqual(z.pairing(psi, z.p.mul(v1, psi), n), 0)
        V = z.p.add(v0, v1, v2)
        Hpsi = z.p.add(z.native.operator(psi, n, 1), z.p.scale(z.p.mul(V, psi), -theta))
        energy = z.pairing(psi, Hpsi, n)/z.pairing(psi, psi, n)
        expected = 2*(F(7, 8)*c*c/4-theta*c/2)/(1+c*c/4)
        bulk_trial = 2*(c*c/4-theta*c/2)/(1+c*c/4)
        self.assertEqual(energy, expected)
        self.assertLess(energy, bulk_trial)
        self.assertLess(bulk_trial, 0)

    def test_fixed_volume_gap_budget_cannot_be_called_uniform(self):
        theta = z.native.lower_bound(2, 1)/8
        self.assertGreater(z.finite_gap_bound(2, 1, theta), 0)
        self.assertEqual(z.finite_gap_bound(5, 1, theta), 0)
        self.assertLess(z.finite_gap_bound(100, 1, theta), 0)
        self.assertFalse(self.cert['energy']['fixed_nonzero_coupling_uniform_gap_proved'])

    def test_wrong_boundary_inverse_is_detected(self):
        data = z.fixture(0)
        wrong = z.local_operator(data['inverse'], (0, 1), 1, (1, 1))
        defect = z.p.add(wrong, z.p.scale(data['source'], -1))
        self.assertGreater(z.pairing(defect, defect, 2), 0)

    def test_invalid_rates_degrees_and_inexact_inputs_are_refused(self):
        for k in (0, -1, True, 0.5):
            with self.assertRaises(ValueError):
                z.omega(k)
        for degrees in ((1,), (0, 2), (True, 2), (1.0, 2)):
            with self.assertRaises(ValueError):
                z.local_operator(z.overlap(0, 1, 2), (0, 1), 1, degrees)
        with self.assertRaises(ValueError):
            z.local_operator(z.constant(1, 3), (0, 1, 2), 1, (1, 1, 1))
        for theta in (z.native.lower_bound(2, 1)/2, 0.1):
            with self.assertRaises(ValueError):
                z.weak_energy_budget(1, theta)
        with self.assertRaises(ValueError):
            z.energy_coefficients(1)
        with self.assertRaises(ValueError):
            z.degrees_for(3, 'periodic')

    def test_runtime_and_assertions_policy_is_inherited(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_source_and_claim_tampering_are_refused(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            bad = copy.deepcopy(self.cert)
            bad['model']['interacting_uniform_gap'] = True
            result.write_text(json.dumps(bad))
            pin.write_text(z.native.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_uniform_asymptotics_do_not_assert_convergent_series_or_clay(self):
        self.assertTrue(self.cert['model']['original_scalar_interaction_retained'])
        self.assertTrue(self.cert['model']['thermodynamic_energy_exists'])
        self.assertTrue(self.cert['model']['uniform_finite_order_energy_asymptotics'])
        self.assertFalse(self.cert['model']['infinite_series_convergence'])
        self.assertFalse(self.cert['model']['interacting_history_limit'])
        self.assertFalse(self.cert['poisson']['numerical_order_five_remainder_constant_evaluated'])
        self.assertIn('4D_CLAY_OPEN', self.cert['claim_status'])


if __name__ == '__main__':
    unittest.main()
