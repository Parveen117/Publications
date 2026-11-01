import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym55_anisotropic_joint_limit as z


class AnisotropicJointLimitTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_digest(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())
        self.assertEqual(self.cert['evidence_scope']['written_results'], 7)

    def test_full_open_window_has_spatial_parameters(self):
        rows = self.cert['parameter_gates']['entire_open_window_controls']
        self.assertEqual(len(rows), 5)
        for row in rows:
            g = z.gate(1, F(row['lam']), row['J'], F(row['R']), F(row['rho']))
            self.assertGreater(g['margin'], 0)
        for lam in [F(0), F(1, 3000), F(9999, 24000000)]:
            J, R, rho, u = z.parameters(1, -lam)
            self.assertLess(u*u, rho)
            self.assertEqual(z.gate(1, lam, J, R, rho)['margin'],
                             z.gate(1, -lam, J, R, rho)['margin'])

    def test_endpoint_abstains_without_claiming_gap_loss(self):
        for J in range(1, 25):
            self.assertIsNone(z.gate(1, F(1, 2400), J, F(10001, 10000), F(10001, 10000)))
        for lam in [F(1, 2400), F(1, 1000)]:
            with self.assertRaises(ValueError):
                z.parameters(1, lam)

    def test_spatial_weight_cannot_be_omitted(self):
        beta, theta, J, R = F(3, 2), F(1, 4096), 8, F(4, 3)
        self.assertIsNotNone(z.old.block_gate(beta, theta, J, R))
        self.assertIsNotNone(z.gate(beta, theta, J, R, F(3, 2)))
        self.assertIsNone(z.gate(beta, theta, J, R, 4))

    def test_rank_two_cell_retains_prior_rate(self):
        c = z.rank_two_cell()
        self.assertEqual(c['margin'], F(7, 72))
        self.assertEqual(c['trace_ratio'], 2)
        self.assertEqual(c['gamma'].lo, z.old.block_gate(F(3, 2), F(1, 4096), 8, F(4, 3))['gamma'].lo)
        self.assertGreater(c['gamma'].lo, F(539403, 10**8))

    def test_thermo_curvature_supplies_both_uniform_budgets(self):
        c, f = z.thermo_cell(), z.native.fixture()
        self.assertEqual(c['beta'], 4*f['marker']**2/f['tau'])
        self.assertEqual(c['trace_max'], f['tau']/4)
        self.assertEqual(c['trace_ratio'], F(169, 36))
        self.assertEqual(c['margin'], F(10249, 98304))
        self.assertGreater(c['gamma'].lo, F(79467076, 10**12))

    def test_actual_profile_bounds_allow_rotation_and_reject_rank_one(self):
        self.assertEqual(self.cert['native_profiles']['fixed_periodic_profile_sites'], 3)
        self.assertTrue(self.cert['native_profiles']['common_diagonalization_not_used'])
        for C in [z.a.diag((0, 0, 3)), z.a.diag((0, F(1, 10), 3)),
                  z.a.diag((0, 2, 5)), z.a.diag((-1, 2, 2))]:
            with self.assertRaises(ValueError):
                z.trace_profile([C], 1, 4)
        self.assertEqual(z.trace_profile([z.a.diag((0, 1, 1))], 1, 2),
                         [z.a.diag((0, 1, 1))])

    def test_fine_step_boundary_budget_has_checked_upper_constant(self):
        self.assertEqual(self.cert['parameter_gates']['fine_step_comparisons'], 8)
        c = z.rank_two_cell()
        for step in [F(1, 3), F(1, 19), F(1, 100)]:
            d = z.step_budget(step, c)
            self.assertTrue(9 <= step*d['r'] < 10)
            self.assertLessEqual(d['side'], c['side']/step)
            self.assertGreater(d['half_time'], d['time'])

    def test_duplicate_block_labels_and_short_strips_are_counted(self):
        d = self.cert['weighted_incidence']
        self.assertEqual(d['layouts'], 18)
        self.assertEqual(d['interior_checks'], 96)
        self.assertEqual(d['boundary_checks'], 168)
        self.assertTrue(d['short_strips_and_duplicate_labels'])

    def test_zero_midpoint_pairs_do_not_require_an_invented_inverse(self):
        d = self.cert['singular_midpoints']
        self.assertEqual(d['admissible_midpoint_pairs'], 14)
        self.assertEqual(d['null_pairs'], 2)
        self.assertEqual(d['actual_ordered_path_identities'], 64)
        self.assertTrue(d['arbitrary_null_fillers_immaterial'])
        H = z.old.path_kernel()
        K = z.a.multiply(H, H)
        self.assertEqual(K[0][3], 0)
        self.assertTrue(all(H[0][j]*H[j][3] == 0 for j in range(4)))

    def test_grouped_endpoint_comparison_avoids_invalid_intermediate(self):
        d = self.cert['singular_midpoints']['old_grouped_endpoint_control']
        self.assertTrue(d['zero_normalizer_rejected'])
        self.assertEqual(d['invalid_intermediate'], [0, 3])
        self.assertLessEqual(F(d['grouped_endpoint_hamming_cost']), 1)
        with self.assertRaises(ValueError):
            z.old.short_bridge(0, 3)

    def test_independent_normalized_histories_obey_trace_sensitive_budget(self):
        rows = self.cert['normalized_refinement']['independent_scaled_two_state_histories']
        self.assertEqual(len(rows), 6)
        for row in rows:
            self.assertLess(F(row['operator_error_upper']), F(row['operator_budget_upper']))
            self.assertLess(F(row['history_error_upper']), F(row['history_budget_upper']))
        c = z.rank_two_cell()
        larger = dict(c, cstar=c['cstar']*2)
        self.assertEqual(z.normalized_budget(4, F(1, 16), 2, larger),
                         2*z.normalized_budget(4, F(1, 16), 2, c))

    def test_refinement_coefficient_reduces_to_previous_isotropic_formula(self):
        c = dict(z.rank_two_cell(), cstar=F(6))
        old = dict(d=c['d'])
        for width, step, T in [(2, F(1, 16), F(2)), (17, F(1, 256), F(3))]:
            self.assertEqual(z.history_budget(width, step, T, 2, 1, c),
                             z.joint.history_budget(width, c['lam'], step, T, 2, 1, old))

    def test_tail_formula_has_independent_finite_difference_controls(self):
        self.assertEqual(self.cert['geometric_tails']['exact_tail_difference_identities'], 54)
        self.assertEqual(self.cert['geometric_tails']['asymmetric_extension_paths'], 144)
        c = z.rank_two_cell()
        N, error = z.tail_radius(3, 6, 2, 2, 1, F(1, 10**6), c)
        self.assertEqual(N, 256)
        self.assertLess(error, F(1, 10**6))
        self.assertLess(z.tail(3, 6, 2, 2, 1, 512, c), z.tail(3, 6, 2, 2, 1, 256, c))

    def test_unrestricted_joint_moduli_and_asymmetric_crop(self):
        rows = self.cert['joint_cutoffs']['joint_cutoff_cells']
        self.assertEqual([r['crop_padding'] for r in rows], [110, 219, 438])
        errors = [F(r['total_error_upper']) for r in rows]
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[-1], F(1, 10**6))
        self.assertTrue(self.cert['joint_cutoffs']['slow_exhaustion_not_silently_dropped'])

    def test_close_boundary_is_not_erased_by_a_tiny_heat_step(self):
        c = z.rank_two_cell()
        b = z.joint_budget(F(1, 2**128), 16, 10**9, 3, 6, 2, 2, 1, c)
        self.assertEqual(b['padding_used'], 16)
        self.assertGreater(b['volume_error'], 1)
        with self.assertRaises(ValueError):
            z.joint_budget(F(1, 2**128), 0, 10**9, 3, 6, 2, 2, 1, c)

    def test_cutoff_varying_profile_can_have_distinct_subsequential_histories(self):
        lo, hi = self.cert['native_profiles']['alternating_profile_failure']
        self.assertGreater(F(lo['lo']), F(hi['hi']))
        self.assertEqual(self.cert['native_profiles']['native_energy_controls'], 9)

    def test_positive_reflection_forms_on_a_singular_kernel_control(self):
        H = z.old.path_kernel()
        S = z.a.multiply(z.a.multiply(H, z.a.diag((1, 2, 3, 4))), H)
        vectors = [[F(1), F(-1), F(0), F(2)], [F(0), F(1), F(-2), F(1)]]
        columns = z.a.transpose(vectors)
        Q = z.a.identity(4)
        for _ in range(6):
            form = z.a.multiply(z.a.multiply(vectors, Q), columns)
            self.assertTrue(z.a.psd(form))
            Q = z.a.multiply(Q, S)

    def test_clock_covariance_does_not_calibrate_physical_time(self):
        self.assertEqual(self.cert['scaling']['clock_rescaling_checks'], 3)
        self.assertFalse(self.cert['scaling']['physical_clock_selected'])

    def test_exact_domains_and_missing_hypotheses_are_rejected(self):
        c = z.rank_two_cell()
        invalid = [
            lambda: z.gate(0, 0, 5, 2, 2),
            lambda: z.gate(1, 0.5, 5, 2, 2),
            lambda: z.gate(1, 0, True, 2, 2),
            lambda: z.gate(1, 0, 5, 1, 2),
            lambda: z.gate(1, 0, 5, 2, 1),
            lambda: z.constants(1, 0, 1, 5, F(3, 2), 2, F(9, 8), 5),
            lambda: z.constants(1, 0, 2, 5, F(3, 2), 2, F(9, 8), 4),
            lambda: z.constants(1, 0, 2, 5, F(3, 2), 2, 2, 5),
            lambda: z.trace_profile([], 1, 2),
            lambda: z.normalized_budget(True, F(1, 2), 1, c),
            lambda: z.normalized_budget(2, 0, 1, c),
            lambda: z.history_budget(2, F(1, 2), 2, 2, 1, c),
            lambda: z.history_budget(2, F(1, 8), 1, 2, 1, c),
            lambda: z.tail_radius(3, 6, 2, 2, 1, 0, c),
            lambda: z.crop_radius(1, F(3, 2)),
            lambda: z.step_budget(2, c)]
        for operation in invalid:
            with self.assertRaises(ValueError):
                operation()

    def test_zero_interaction_has_no_spatial_or_operator_error(self):
        J, R, rho, u = z.parameters(1, 0)
        c = z.constants(1, 0, 2, J, R, rho, u, 5)
        self.assertEqual(c['side'], 0)
        self.assertEqual(z.normalized_budget(20, F(1, 16), 3, c), 0)
        self.assertEqual(z.history_budget(20, F(1, 16), 0, 0, 1, c), 0)

    def test_tampered_certificate_and_upstream_pin_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, pin, sources = [Path(tmp)/name for name in ['result', 'pin', 'sources']]
            changed = copy.deepcopy(self.cert)
            changed['claim_status'] = 'CLAY_SOLVED'
            result.write_text(json.dumps(changed))
            pin.write_text(z.old.old.canonical_sha(changed))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)
            p = json.loads(z.SOURCES.read_text())
            key = next(iter(p['upstream_sha256']))
            p['upstream_sha256'][key] = '0'*64
            sources.write_text(json.dumps(p))
            original = z.SOURCES
            try:
                z.SOURCES = sources
                with self.assertRaises(ValueError):
                    z.source_checks()
            finally:
                z.SOURCES = original


if __name__ == '__main__':
    unittest.main()
