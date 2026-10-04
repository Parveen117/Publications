import copy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym58_infinite_volume_histories as z


class InfiniteVolumeHistoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_matches_pin(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_two_boundary_bonds_restore_singleton_noise(self):
        w = z.prior.previous.witness(1)
        self.assertFalse(z.prior.operator(w, 1, 1))
        self.assertEqual(z.marginal(w, (0,), 1), w)
        f = z.embed(w, (0,), (0, 1))
        self.assertEqual(z.prior.operator(f, 2, 1), z.p.scale(f, F(1, 2)))

    def test_missing_one_neighbour_layer_changes_heat_derivative(self):
        f = z.p.COORD[0]
        source, incomplete = (0,), (0, 1, 2)
        observed = z.embed(f, source, incomplete)
        actual = z.prior.operator(observed, 3, 1)
        target = z.embed(z.marginal(f, source, 1), source, incomplete)
        self.assertNotEqual(actual, target)
        full = (-1, 0, 1, 2)
        self.assertEqual(z.prior.operator(z.embed(f, source, full), 4, 1),
                         z.embed(z.marginal(f, source, 1), source, full))

    def test_arbitrary_disconnected_marginals_are_consistent(self):
        f = z.p.mul(z.prior.lift(z.p.COORD[0], 0, 2),
                    z.prior.lift(z.p.COORD[2], 1, 2))
        source, target = (-4, 1), (-5, -4, 0, 1, 2)
        self.assertEqual(z.marginal(z.embed(f, source, target), target, F(3, 7)),
                         z.embed(z.marginal(f, source, F(3, 7)), source, target))
        self.assertEqual(self.cert['locality']['arbitrary_subset_consistency_identities'], 204)

    def test_empty_support_is_constant_and_zero_generator(self):
        self.assertFalse(z.marginal({(): F(9)}, (), 1))
        f = z.embed({(): F(9)}, (), (-1, 0, 1))
        self.assertFalse(z.marginal(f, (-1, 0, 1), 1))
        self.assertEqual(z.prior.phi(f, 3), 9)

    def test_incident_edges_are_counted_once(self):
        self.assertEqual(z.incident_edges((0,)), ((-1, 0), (0, 1)))
        self.assertEqual(z.active_weight((0, 1), 1), 5)
        self.assertEqual(z.active_weight((0, 2), 1), 6)
        self.assertEqual(sum(w for w, word in z.directions((0, 1), 1)), 5)

    def test_support_and_native_squares_agree_in_all_checked_volumes(self):
        self.assertEqual(self.cert['locality']['open_volume_generator_identities'], 408)
        self.assertEqual(self.cert['locality']['power_identities'], 90)
        self.assertEqual(self.cert['energy']['independent_incident_square_factorizations'], 15)

    def test_gap_is_inherited_without_replacing_boundary_measure(self):
        self.assertEqual(self.cert['energy']['rational_gap_inequalities'], 15)
        self.assertEqual(self.cert['energy']['invariant_reference_identities'], 366)
        self.assertEqual(F(self.cert['energy']['unit_edge_gap_lower']), F(1, 630))
        self.assertEqual(z.prior.lower_bound(2, F(1, 100)), F(1, 42000))

    def test_zero_edge_rate_keeps_the_stationary_obstruction(self):
        w = z.prior.previous.witness(1)
        self.assertFalse(z.marginal(w, (0,), 0))
        self.assertEqual(z.prior.variance(w, 1), F(1, 12))
        with self.assertRaises(ValueError):
            z.prior.lower_bound(2, 0)

    def test_singleton_sharp_rate_uses_both_boundary_edges(self):
        rows = self.cert['boundary']['rows']
        self.assertEqual(F(next(r for r in rows if r['k'] == '1')['singleton_gap']), F(3, 8))
        self.assertEqual(F(next(r for r in rows if r['k'] == '1/10')['singleton_gap']), F(1, 10))

    def test_multitime_coefficients_have_no_remaining_volume_dependence(self):
        self.assertEqual(self.cert['histories']['multitime_coefficient_volume_identities'], 128)
        self.assertGreater(self.cert['histories']['nonzero_local_coefficients'], 0)

    def test_time_series_sign_and_factorials_match_exact_linear_heat(self):
        f = z.p.COORD[0]
        for r in range(7):
            coefficient = z.history_coefficient((f, f), (r,), (0,), 1)
            self.assertEqual(coefficient, F(-3, 8)**r/(4*factorial(r)))

    def test_adjacent_history_differs_from_separated_history(self):
        self.assertEqual(F(self.cert['histories']['adjacent_heat_correlation']), F(5, 512))
        self.assertEqual(F(self.cert['histories']['separated_heat_correlation']), F(1, 128))
        f = z.p.mul(z.prior.lift(z.p.COORD[0], 0, 2),
                    z.prior.lift(z.p.COORD[0], 1, 2))
        # Cross terms have zero initial expectation here, but the second
        # time derivative already detects the correlated kinetic source.
        adjacent = z.history_coefficient((f, f), (2,), (0, 1), 1)
        separated = z.history_coefficient((f, f), (2,), (0, 2), 1)
        self.assertEqual(adjacent-separated, F(1, 512))

    def test_reflected_heat_gram_is_positive_and_wrong_sign_is_refused(self):
        gram = [[F(v) for v in row] for row in self.cert['histories']['reflection_gram']]
        self.assertTrue(z.a.psd(gram))
        self.assertFalse(z.a.psd([[F(1, 4), F(1, 2)], [F(1, 2), F(1, 4)]]))
        self.assertFalse(self.cert['histories']['finite_Taylor_polynomials_claimed_positive'])

    def test_global_label_fourth_moment_retains_volume_factor(self):
        for row in self.cert['refinement']['global_label_moment_checks']:
            N, k = row['width'], F(row['k'])
            self.assertEqual(F(row['second_coefficient']), -(1+2*k)/8)
            self.assertEqual(F(row['fourth_coefficient']), (N+k*(N-1))*(1+2*k)/384)

    def test_refinement_gate_decreases_and_has_zero_degree_limit(self):
        rows = self.cert['refinement']['joint_gate_budgets']
        bounds = [F(r['upper']) for r in rows]
        self.assertEqual(bounds, sorted(bounds, reverse=True))
        self.assertLess(bounds[-1], F(1, 500))
        self.assertEqual(z.refinement_budget((0,), 1, 3, 10, 0, 1), 0)
        self.assertEqual(z.refinement_budget((0,), 1, 3, 10, 2, 0), 0)

    def test_fast_volume_sequence_stays_away_from_heat_target(self):
        lower = z.fast_volume_multiplier_lower(100)
        self.assertEqual(lower, F(19399, 19999))
        self.assertGreater(lower, F(96, 100))
        self.assertGreater(lower, F(self.cert['refinement']['exact_heat_multiplier_upper']))
        self.assertGreater(z.fast_volume_multiplier_lower(1000), lower)
        self.assertFalse(self.cert['refinement']['unrestricted_joint_limit'])
        self.assertFalse(self.cert['refinement']['volume_stabilization_at_finite_global_step'])

    def test_inexact_and_inconsistent_inputs_are_refused(self):
        for support in ((1, 0), (0, 0), (True,), (0.5,)):
            with self.assertRaises(ValueError):
                z.sites(support)
        for k in (-1, 0.5, True):
            with self.assertRaises(ValueError):
                z.marginal(z.p.COORD[0], (0,), k)
        with self.assertRaises(ValueError):
            z.embed(z.p.COORD[0], (0,), (1,))
        with self.assertRaises(ValueError):
            z.marginal(z.p.COORD[0], (0, 1), 1)
        for args in ((2, 1, 1, 1), (3, 0, 1, 1), (3, 1, -1, 1), (3, 1, 1, -1)):
            with self.assertRaises(ValueError):
                z.refinement_budget((0,), 1, *args)
        with self.assertRaises(ValueError):
            z.history_coefficient((z.p.COORD[0],), (1,), (0,), 1)

    def test_source_and_claim_tampering_are_refused(self):
        with patch.object(z.prior, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            bad = copy.deepcopy(self.cert)
            bad['refinement']['unrestricted_joint_limit'] = True
            result.write_text(json.dumps(bad))
            pin.write_text(z.prior.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_new_history_scope_does_not_import_the_scalar_model(self):
        self.assertTrue(self.cert['model']['infinite_volume_histories'])
        self.assertTrue(self.cert['model']['row_markov_for_this_model'])
        self.assertFalse(self.cert['model']['scalar_interaction_imported'])
        self.assertFalse(self.cert['model']['physical_clock_calibrated'])
        self.assertFalse(self.cert['model']['spatial_continuum_constructed'])
        self.assertIn('4D_CLAY_OPEN', self.cert['claim_status'])


if __name__ == '__main__':
    unittest.main()
