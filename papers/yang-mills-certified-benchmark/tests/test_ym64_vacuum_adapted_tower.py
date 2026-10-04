import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym64_vacuum_adapted_tower as z


class VacuumAdaptedTowerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.H, cls.C, cls.B, cls.D = z.fixture()

    def test_fresh_result_and_canonical_hash_are_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_is_enforced(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_183_upstream_sources_and_four_local_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 183)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn('papers/generalized-euler-evolution/GE4_MINIMAL_OBSERVER_MEMORY.md',
                      manifest['upstream_sha256'])

    def test_changed_upstream_sources_are_refused(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()

    def test_rehashed_continuum_claim_tampering_still_fails_fresh_comparison(self):
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            changed = copy.deepcopy(self.cert)
            changed['model']['continuum_R4_QFT'] = True
            result.write_text(json.dumps(changed))
            pin.write_text(z.p.canonical_sha(changed))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_conditioning_is_orthogonal_in_the_actual_weighted_pairing(self):
        I, J, masses = z.partition((1, 2, 3, 4), ((0, 1), (2, 3)))
        P = z.a.multiply(I, J)
        W = z.a.diag((1, 2, 3, 4))
        self.assertEqual(z.a.multiply(P, P), P)
        self.assertEqual(z.a.multiply(W, P), z.a.multiply(z.a.transpose(P), W))
        self.assertEqual(P[0][:2], [F(1, 3), F(2, 3)])
        self.assertNotEqual(P, z.a.transpose(P))

    def test_nested_marginals_and_constants_are_preserved(self):
        I, J, weights = z.partition((1, 2, 3, 4), ((0, 1), (2,), (3,)))
        U, Z, twice = z.partition(weights, ((0,), (1, 2)))
        direct_I, direct_J, direct_mass = z.partition((1, 2, 3, 4), ((0, 1), (2, 3)))
        self.assertEqual(z.a.multiply(I, U), direct_I)
        self.assertEqual(z.a.multiply(Z, J), direct_J)
        self.assertEqual(twice, direct_mass)
        self.assertTrue(all(sum(row) == 1 for row in z.a.multiply(I, J)))

    def test_overlapping_missing_empty_and_inexact_partitions_are_refused(self):
        for groups in (((0, 1), (1, 2)), ((0,), (1,)), ((0, 1, 2), ()), ((False,), (1, 2))):
            with self.assertRaises(ValueError):
                z.partition((1, 2, 3), groups)
        for weights in ((1, 0), (1, -1), (1, 0.5)):
            with self.assertRaises(ValueError):
                z.partition(weights, ((0,), (1,)))

    def test_free_cut_overlap_degrades_while_adapted_hidden_gap_stays_one(self):
        rows = self.cert['obstructions']['free_overlap_can_vanish']
        values = [F(row['free_cut_hidden_gap']) for row in rows]
        self.assertTrue(all(x > y > 0 for x, y in zip(values, values[1:])))
        self.assertLess(values[-1], F(1, 1000))
        self.assertEqual({row['adapted_hidden_gap'] for row in rows}, {1})

    def test_operator_fixture_has_full_gap_and_exact_retained_vacuum(self):
        self.assertTrue(z.a.psd(z.a.add(self.H, z.a.diag((0, 1, 1, 1, 1, 1)), -1)))
        self.assertEqual(self.H[0], [F(0)]*6)
        self.assertTrue(z.a.psd(z.a.add(self.D, z.a.identity(3), -2)))

    def test_schur_subtraction_equals_independent_compressed_inverse(self):
        s = F(4, 7)
        shifted = z.a.add(self.H, z.a.identity(6), s)
        compressed_inverse = z.sub(z.a.inverse(shifted), (0, 1, 2))
        correct = z.schur(shifted, (0, 1, 2))
        self.assertEqual(z.a.inverse(compressed_inverse), correct)
        wrong = z.a.add(z.a.add(self.C, z.a.identity(3), s), z.self_energy(self.B, self.D, s))
        self.assertNotEqual(correct, wrong)

    def test_three_successive_spectral_blocks_equal_one_direct_block(self):
        for s in (F(1, 11), F(1), F(7, 3)):
            first = z.effective(self.H, (0, 1, 2, 3), s)
            second = z.effective(first, (0, 1, 2), s)
            third = z.effective(second, (0, 1), s)
            self.assertEqual(third, z.effective(self.H, (0, 1), s))

    def test_static_reset_cannot_replace_the_spectral_family(self):
        rows = self.cert['operators']['nested_blocks_and_reset_counterexamples']
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(row['direct'] != row['static_reset'] for row in rows))

    def test_weighted_conditional_resolvents_compose_without_independence(self):
        data = self.cert['conditional']
        self.assertTrue(data['conditional_tower'])
        self.assertTrue(data['unweighted_averaging_differs'])
        for row in data['exact_nested_resolvents']:
            self.assertTrue(all(sum(F(x) for x in line) == 0 for line in row['coarse_operator']))

    def test_every_inverse_moment_obeys_the_hidden_gap_bound(self):
        base = z.moment(self.B, self.D, 0)
        for n in range(9):
            value = z.moment(self.B, self.D, n)
            self.assertTrue(z.a.psd(value))
            self.assertTrue(z.a.psd(z.a.add(z.a.scale(base, F(1, 2)**n), value, -1)))

    def test_positive_moment_operators_need_not_commute(self):
        first, second = z.moment(self.B, self.D, 0), z.moment(self.B, self.D, 1)
        self.assertNotEqual(z.a.multiply(first, second), z.a.multiply(second, first))

    def test_hankel_shift_and_gap_localizer_are_positive(self):
        moments = [z.moment(self.B, self.D, n) for n in range(6)]
        H0, H1 = z.hankel(moments, 3), z.hankel(moments, 3, 1)
        self.assertTrue(z.a.psd(H0))
        self.assertTrue(z.a.psd(H1))
        self.assertTrue(z.a.psd(z.a.add(z.a.scale(H0, F(1, 2)), H1, -1)))
        self.assertEqual(z.p.y37.inertia(H0), (3, 0, 6))

    def test_positive_individual_coefficients_do_not_suffice_for_a_moment_tower(self):
        impossible = [[[F(1)]], [[F(2)]], [[F(1)]]]
        self.assertTrue(all(z.a.psd(value) for value in impossible))
        self.assertFalse(z.a.psd(z.hankel(impossible, 2)))

    def test_form_taylor_error_is_bounded_on_both_sides_of_zero(self):
        moments = [z.moment(self.B, self.D, n) for n in range(6)]
        for s in (F(-1, 3), F(1, 3)):
            trunc = [[F(0)]*3 for _ in range(3)]
            for n, value in enumerate(moments):
                trunc = z.a.add(trunc, value, (-s)**n)
            error = z.a.add(z.self_energy(self.B, self.D, s), trunc, -1)
            bound = z.a.scale(moments[0], z.tower_budget(1, s, 5)['remainder_factor'])
            self.assertTrue(z.a.psd(z.a.add(bound, error)))
            self.assertTrue(z.a.psd(z.a.add(bound, error, -1)))

    def test_positive_remainder_bound_does_not_require_series_convergence(self):
        s = F(3)
        base = z.moment(self.B, self.D, 0)
        error = z.a.add(base, z.self_energy(self.B, self.D, s), -1)
        self.assertTrue(z.a.psd(error))
        self.assertTrue(z.a.psd(z.a.add(z.a.scale(base, F(3, 4)), error, -1)))
        with self.assertRaises(ValueError):
            z.tower_budget(1, s, 0)

    def test_gap_endpoint_zero_gap_and_inexact_tower_inputs_are_refused(self):
        for gap, s, order in ((1, 1, 2), (1, -1, 2), (0, 0, 2), (-1, 0, 2),
                              (1.0, 0, 2), (1, 0.1, 2), (1, 0, True), (1, 0, -1)):
            with self.assertRaises(ValueError):
                z.tower_budget(gap, s, order)

    def test_integrated_memory_tails_decay_in_energy_relative_order(self):
        base = z.moment(self.B, self.D, 0)
        for row in self.cert['operators']['integrated_memory_tails']:
            tail = [[F(v) for v in line] for line in row['tail']]
            self.assertTrue(z.a.psd(tail))
            self.assertTrue(z.a.psd(z.a.add(z.a.scale(base, F(row['relative_bound'])), tail, -1)))

    def test_clock_rescaling_changes_each_order_by_its_own_dimension(self):
        kappa = F(5, 2)
        for n in range(5):
            changed = z.moment(z.a.scale(self.B, kappa), z.a.scale(self.D, kappa), n)
            self.assertEqual(changed, z.a.scale(z.moment(self.B, self.D, n), kappa**(1-n)))

    def test_zero_coupling_closes_the_retained_resolvent_exactly(self):
        free = z.a.diag((0, 2, 3, 5))
        for s in (F(1, 7), F(2)):
            self.assertEqual(z.effective(free, (0, 1), s), z.a.diag((0, 2)))
        zero = [[F(0)]*3 for _ in range(2)]
        self.assertEqual(z.moment(zero, self.D, 4), [[F(0)]*2 for _ in range(2)])

    def test_native_pointwise_gradient_pullbacks_retain_the_linear_clock(self):
        rows = self.cert['native']['pointwise_carre_du_champ_pullbacks']
        self.assertEqual([row['factor'] for row in rows], [2, 2, 3])
        self.assertTrue(all(row['terms'] > 0 for row in rows))

    def test_native_positive_and_inverse_word_derivatives_match_polynomials(self):
        self.assertEqual(self.cert['native']['word_derivative_sign_checks'], 12)
        for axis in (True, 1.0, 0, 4):
            with self.assertRaises(ValueError):
                z.word_derivative({}, z.model.ORIENTED_WORD, 0, axis)

    def test_physical_coarse_loop_has_nonzero_adapted_first_order_leakage(self):
        data = self.cert['native']
        self.assertEqual(data['fine_boundary_links'], 12)
        self.assertEqual(data['adjacent_plaquettes'], 44)
        self.assertEqual(F(data['gamma_V_probe_at_rational_configuration']), F(16, 25))
        self.assertEqual(F(data['first_order_hidden_leakage']), F(-32, 75))

    def test_local_gauge_change_preserves_the_physical_response(self):
        self.assertTrue(self.cert['native']['local_gauge_transformation_preserves_response'])

    def test_distant_native_readouts_have_disjoint_leading_response_supports(self):
        self.assertTrue(self.cert['native']['distant_leading_response_supports_disjoint'])

    def test_uniform_gap_does_not_imply_spatial_locality(self):
        rows = self.cert['obstructions']['gap_alone_does_not_imply_spatial_locality']
        self.assertEqual([row['separation'] for row in rows], [4, 8, 16])
        self.assertTrue(all(abs(F(row['off_diagonal'])) > 1 for row in rows))
        self.assertEqual({row['gap_lower'] for row in rows}, {1})

    def test_singular_hidden_inverse_and_bad_retained_sets_are_refused(self):
        with self.assertRaises(ValueError):
            z.schur(z.a.diag((1, 0)), (0,))
        for indices in ((), (0, 0), (False,), (6,)):
            with self.assertRaises(ValueError):
                z.schur(self.H, indices)
        with self.assertRaises(ValueError):
            z.effective(self.H, (0, 1), 0)

    def test_scope_keeps_relative_bounds_separate_from_continuum_and_locality(self):
        model = self.cert['model']
        for key in ('YM62_operator_retained', 'actual_vacuum_weighted_cut',
                    'full_constant_vacuum_retained', 'uniform_hidden_gap_from_YM62',
                    'nested_resolvent_composition', 'positive_inverse_moment_tower',
                    'uniform_energy_relative_memory_bounds', 'local_first_order_physical_leakage'):
            self.assertTrue(model[key])
        for key in ('full_spatial_locality_bound', 'uniform_absolute_memory_norm',
                    'explicit_interacting_vacuum_computed', 'exact_scalar_RG',
                    'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
