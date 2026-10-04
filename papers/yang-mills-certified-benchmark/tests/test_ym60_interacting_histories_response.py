import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym60_interacting_histories_response as z


class InteractingHistoriesResponseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_exact_certificate_is_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_native_brackets_recover_each_site_including_endpoints(self):
        self.assertEqual(self.cert['brackets']['site_coordinate_bracket_checks'], 36)
        n, site = 2, 1
        f = z.native.lift(z.p.COORD[0], site, n)
        D = lambda g: z.native.derivative(g, site, 1)
        B = lambda g: z.native.bond(g, 0)
        bracket = z.p.add(D(B(f)), z.p.scale(B(D(f)), -1))
        self.assertEqual(bracket, z.p.scale(z.native.derivative(f, site, 3), -1))
        self.assertNotEqual(bracket, z.native.derivative(f, site, 3))

    def test_scalar_interaction_and_kinetic_symmetries_are_both_checked(self):
        self.assertEqual(self.cert['symmetry']['kinetic_symmetry_intertwinings'], 48)
        self.assertEqual(self.cert['symmetry']['bilinear_average_identities'], 64)
        self.assertEqual(self.cert['symmetry']['zero_mean_identities'], 20)

    def test_common_right_action_alone_does_not_force_scalar_cross_covariance(self):
        self.assertTrue(self.cert['symmetry']['right_symmetry_alone_insufficient'])
        J, L = z.y.qmatrix(z.UNITS[1]), z.y.qmatrix(z.UNITS[2])
        self.assertNotEqual(z.a.multiply(J, L), z.a.multiply(L, J))

    def test_bilinear_average_uses_both_left_and_right_actions(self):
        f = z.p.mul(z.native.lift(z.p.COORD[2], 0, 2), z.native.lift(z.p.COORD[2], 1, 2))
        self.assertEqual(z.quaternion_average(f, 2), z.p.scale(z.previous.overlap(0, 1, 2), F(1, 4)))
        off = z.p.mul(z.native.lift(z.p.COORD[2], 0, 2), z.native.lift(z.p.COORD[1], 1, 2))
        self.assertFalse(z.quaternion_average(off, 2))

    def test_weighted_ground_cancellation_and_double_commutator(self):
        self.assertEqual(self.cert['cancellation']['weighted_cancellation_identities'], 18)
        self.assertEqual(self.cert['cancellation']['double_commutator_identities'], 18)

    def test_a_chosen_positive_polynomial_is_not_the_true_ground_vector(self):
        self.assertGreater(F(self.cert['cancellation']['chosen_polynomial_is_not_ground_vector_residual_norm2']), 0)
        self.assertTrue(self.cert['forms']['scalar_vacuum_not_numerically_approximated'])

    def test_neighbour_form_entry_is_negative_for_positive_overlap(self):
        n = 2
        v = z.previous.overlap(0, 1, n)
        density = z.p.add(z.previous.constant(1, n), z.p.scale(v, F(1, 2)))
        mean = lambda f: z.native.phi(z.p.mul(density, f), n)
        f, g = z.native.lift(z.p.COORD[0], 0, n), z.native.lift(z.p.COORD[0], 1, n)
        expected = -mean(v)/32
        self.assertEqual(mean(z.gamma(f, g, (0, 1), 1)), expected)
        self.assertEqual(expected, F(-1, 256))

    def test_weighted_forms_include_cross_coordinates_and_disjoint_sites(self):
        self.assertEqual(self.cert['forms']['weighted_form_entries'], 108)
        f, g = z.native.lift(z.p.COORD[0], 0, 3), z.native.lift(z.p.COORD[0], 2, 3)
        self.assertFalse(z.gamma(f, g, (0, 1, 2), 1))

    def test_physical_boundary_correction_does_not_disappear_at_finite_width(self):
        self.assertEqual(z.finite_trial_rate(2, 1, 0), F(1, 4))
        self.assertEqual(z.finite_trial_rate(3, 1, 0), F(5, 16))
        self.assertEqual(z.trial_rate(1, 0), F(3, 8))
        self.assertNotEqual(z.finite_trial_rate(100, 1, 0), z.trial_rate(1, 0))

    def test_correlated_response_is_computed_from_the_three_local_inverse_sources(self):
        self.assertEqual(F(self.cert['response']['bulk_energy_quartic']), F(-1123, 47376))
        self.assertEqual(F(self.cert['response']['bulk_correlation_cubic']), F(1123, 11844))

    def test_response_to_trial_rate_series_has_the_exact_signs(self):
        series = tuple(map(F, self.cert['response']['trial_rate_series']))
        self.assertEqual(series, (F(3, 8), F(-1, 4), F(1, 8), F(-5207, 47376)))
        c = F(1, 7)
        self.assertEqual(z.trial_rate(1, c), F(5, 16))

    def test_formal_differentiation_of_one_remainder_is_not_enough(self):
        self.assertEqual(z.secant_error_order(6, 3), 3)
        self.assertEqual(z.secant_error_order(10, 5), 5)
        self.assertEqual(z.secant_error_order(10, 4), 4)
        with self.assertRaises(ValueError):
            z.secant_error_order(6, 6)

    def test_trial_upper_rate_strictly_decreases_with_aligned_response(self):
        for k in (F(1, 10), F(1), F(3)):
            x, y = F(1, 7), F(3, 5)
            difference = z.trial_rate(k, x)-z.trial_rate(k, y)
            self.assertEqual(difference, (1+3*k)*(y-x)/(8*(1+x)*(1+y)))
            self.assertGreater(difference, 0)

    def test_first_moment_does_not_bound_the_threshold_from_below(self):
        controls = self.cert['histories']['equal_first_moment_counterexamples']
        self.assertEqual([F(row['mean']) for row in controls], [F(1)]*4)
        self.assertEqual([F(row['threshold']) for row in controls], [F(1, 2), F(1, 10), F(1, 1000), F(0)])

    def test_positive_spectral_mixture_can_retain_a_zero_atom(self):
        result = z.spectral_summary(((0, F(3, 4)), (4, F(1, 4))))
        self.assertEqual(result['mean'], 1)
        self.assertEqual(result['threshold'], 0)
        self.assertTrue(result['positive_moment_matrix'])

    def test_reflection_positive_mixture_need_not_be_time_zero_markov(self):
        self.assertTrue(self.cert['histories']['reflected_mixture_matrix_psd'])
        self.assertEqual(F(self.cert['histories']['markov_composition_defect']), F(1, 64))
        self.assertNotEqual((F(1, 4)+F(1, 16))/2, ((F(1, 2)+F(1, 4))/2)**2)

    def test_local_time_estimate_is_independent_of_spatial_origin(self):
        f = z.p.add(z.native.lift(z.p.COORD[0], 0, 2), z.native.lift(z.p.COORD[0], 1, 2))
        self.assertEqual(z.increment_budget(f, (0, 1), 1, F(1, 7)), F(2, 7))
        self.assertEqual(z.increment_budget(f, (-300, -299), 1, F(1, 7)), F(2, 7))
        self.assertEqual(z.increment_budget(f, (0, 1), 1, 0), 0)

    def test_false_signed_permutation_is_refused(self):
        bad = z.a.identity(4)
        bad[0][0] = F(2)
        with self.assertRaises(ValueError):
            z.signed_transform(z.p.COORD[0], bad, 1)
        bad = [list(row) for row in z.a.identity(4)]
        bad[1] = list(bad[0])
        with self.assertRaises(ValueError):
            z.signed_transform(z.p.COORD[0], bad, 1)

    def test_unphysical_inputs_and_inexact_floats_are_refused(self):
        for k, x in ((0, 0), (1, -1), (1, F(3, 2)), (1.0, 0), (1, True)):
            with self.assertRaises(ValueError):
                z.trial_rate(k, x)
        with self.assertRaises(ValueError):
            z.finite_trial_rate(1, 1, 0)
        with self.assertRaises(ValueError):
            z.increment_budget(z.p.COORD[0], (0,), 1, -1)
        for atoms in ((), ((-1, 1),), ((1, -1),), ((1, 0),), ((1.0, 1),)):
            with self.assertRaises(ValueError):
                z.spectral_summary(atoms)

    def test_runtime_policy_rejects_non_python_312(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_upstream_and_claim_tampering_are_refused(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            bad = copy.deepcopy(self.cert)
            bad['model']['uniform_interacting_lower_gap'] = True
            result.write_text(json.dumps(bad))
            pin.write_text(z.native.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_historical_source_count_and_read_only_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 153)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn(str(z.previous.RESULT.relative_to(z.ROOT)), manifest['upstream_sha256'])

    def test_subsequential_history_result_does_not_claim_markov_uniqueness_or_lower_gap(self):
        model = self.cert['model']
        self.assertTrue(model['original_scalar_interaction_retained'])
        self.assertTrue(model['reflection_positive'])
        self.assertTrue(model['continuous_heat_time'])
        self.assertIn('Subsequential', model['interacting_history_limit'])
        for key in ('history_limit_unique', 'time_zero_markov', 'uniform_interacting_lower_gap', 'physical_4d_yang_mills'):
            self.assertFalse(model[key])
        self.assertFalse(self.cert['histories']['sample_path_continuity_claimed'])
        self.assertIn('LOWER_GAP_OPEN', self.cert['claim_status'])


if __name__ == '__main__':
    unittest.main()
