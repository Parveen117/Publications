"""YM93: native identities, complete tails, continuation and actual budgets."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym93_bochner_continuation as y


class NativeBochnerTests(unittest.TestCase):
    def test_covariant_vector_and_scalar_identities_on_three_carriers(self):
        rows = [y.vector_control(n) for n in (1, 2, 3)]
        self.assertEqual(sum(r['covariant_one_form_components'] for r in rows), 18)
        self.assertEqual(sum(r['scalar_partial_Bochner_identities'] for r in rows), 6)
        self.assertEqual(len({r['potential_sha256'] for r in rows}), 3)

    def test_scalar_laplacian_and_missing_ricci_or_drift_are_detected(self):
        for n in (1, 2, 3):
            row = y.vector_control(n)
            for key in ('components_changed_by_scalar_Laplacian',
                        'components_changed_by_dropping_Ricci',
                        'components_changed_by_dropping_drift'):
                self.assertEqual(row[key], 3*n)
            self.assertEqual(row['native_Ricci_damping'], F(1, 2))

    def test_second_derivative_commutators_keep_the_two_product_terms(self):
        row = y.derivative_control()
        self.assertEqual(row['unit_direction_cases'], 8)
        self.assertEqual(row['second_plaquette_derivative_factor'], -F(1, 4))
        self.assertTrue(row['electric_commutators_and_full_product_rules'])
        self.assertFalse(row['test_wave_is_actual_Wilson_ground_state'])

    def test_all_coupling_gradient_helper_does_not_assert_a_gap(self):
        for t in (0, F(1, 10000), F(1, 100), 1, 10, 1000):
            self.assertEqual(y.gradient_upper(t), 4*t)
        self.assertFalse(y.MODEL['all_coupling_gradient_bound_implies_all_coupling_gap'])
        with self.assertRaises(ValueError):
            y.actual_budget(1)

    def test_native_comparisons_keep_their_carrier_and_scope(self):
        for n in (0, 4, True, F(2)):
            with self.assertRaises(ValueError):
                y.comparison_chi(n)
        self.assertIn('not_Wilson_vacuum', y.vector_control(3)['kind'])
        with self.assertRaises(ValueError):
            y.gradient_upper(-1)


class HeatAndDiagonalTests(unittest.TestCase):
    def test_time_four_kernel_bound_includes_all_representation_channels(self):
        row = y.heat_control()
        self.assertGreater(row['exp3_lower'], 20)
        self.assertGreater(row['exp4_lower'], 50)
        self.assertEqual(row['channels_n_at_least_2_full_tail'], F(10977, 2941225))
        self.assertEqual(row['complete_kernel_deviation_upper'], F(599222, 2941225))
        self.assertLess(row['complete_kernel_deviation_upper'], F(21, 100))

    def test_closed_geometric_tail_agrees_with_independent_finite_sums(self):
        for q in (F(0), F(1, 50), F(1, 3)):
            for start in (0, 1, 2, 5):
                tail = y.squared_geometric_tail(q, start)
                finite = sum(((n+1)**2*q**n for n in range(start, start+12)), F(0))
                self.assertEqual(tail, finite+y.squared_geometric_tail(q, start+12))
                self.assertGreaterEqual(tail, finite)

    def test_positive_exponential_enclosure_retains_a_nonzero_tail(self):
        for x in (F(1, 100), F(4, 25), F(88, 175), F(2)):
            self.assertGreater(y.exp_upper(x), y.exp_lower(x, 30))
            self.assertGreater(y.exp_upper(x), 1+x+x*x/2)
        self.assertEqual(y.exp_upper(0), 1)

    def test_exact_endpoint_duhamel_denominator_is_positive(self):
        row = y.diagonal_budget(F(1, 100))
        self.assertEqual(row['heat_exp_upper'], F(2083, 1775))
        self.assertEqual(row['square_Harnack_ratio_upper'], F(126547, 76475))
        self.assertEqual(row['Duhamel_denominator'], F(102957, 177500))
        self.assertEqual(row['adaptive_diagonal_Hessian_upper'], F(4081841509, 28120130625))
        self.assertLess(row['adaptive_diagonal_Hessian_upper'], F(3, 20))

    def test_ratio_monotonicity_supports_closed_interval_diagonal_bound(self):
        last = F(0)
        for n in range(1, 101):
            t = F(n, 10000)
            row = y.diagonal_budget(t)
            ratio = row['adaptive_diagonal_Hessian_upper']/t
            self.assertGreater(ratio, last)
            self.assertLess(ratio, 15)
            last = ratio
        self.assertEqual(y.diagonal_budget(0)['adaptive_diagonal_Hessian_upper'], 0)

    def test_relative_second_wave_derivative_is_not_log_hessian_alone(self):
        t = F(1, 100)
        row = y.diagonal_budget(t)
        self.assertEqual(row['adaptive_diagonal_Hessian_upper'],
                         row['relative_second_wave_derivative_upper']+(4*t)**2)
        expected = row['square_Harnack_ratio_upper']*(row['heat_exp_upper']-1)*(F(1, 4)+4*t)/row['Duhamel_denominator']
        self.assertEqual(row['relative_second_wave_derivative_upper'], expected)
        self.assertGreater(row['relative_second_wave_derivative_upper'],
                           row['square_Harnack_ratio_upper']*(row['heat_exp_upper']-1)/4/row['Duhamel_denominator'])

    def test_heat_enclosures_reject_divergent_or_inexact_inputs(self):
        for value in (-1, 3, 4):
            with self.assertRaises(ValueError):
                y.exp_upper(value)
        for q in (-1, 1, 2):
            with self.assertRaises(ValueError):
                y.squared_geometric_tail(q, 2)
        for value in (True, 0.2):
            with self.assertRaises(TypeError):
                y.exp_upper(value)
        with self.assertRaises(ValueError):
            y.exp_lower(1, -1)


class ContinuationTests(unittest.TestCase):
    def test_constant_barrier_is_strictly_forbidden_on_the_entire_window(self):
        for n in range(101):
            t = F(n, 10000)
            self.assertLessEqual(y.barrier_polynomial(t, y.BARRIER), -F(13, 62500))
        self.assertEqual(y.barrier_polynomial(F(1, 100), y.BARRIER), -F(13, 62500))

    def test_large_branch_remains_admissible_without_continuity(self):
        t = F(1, 100)
        self.assertGreater(y.barrier_polynomial(t, 1), 0)
        self.assertGreater(F(1), y.remainder_budget(t)['mixed_remainder_upper'])
        self.assertFalse(y.remainder_budget(t)['quadratic_inequality_alone_selects_small_branch'])

    def test_refinement_uses_the_previously_proved_constant_barrier(self):
        for n in range(101):
            t = F(n, 10000)
            row = y.remainder_budget(t)
            self.assertGreaterEqual(row['refinement_denominator'], F(77, 250))
            self.assertLessEqual(row['refined_remainder_upper'], F(19500, 77)*t*t)
            self.assertLessEqual(row['refined_remainder_upper'], 260*t*t)
        self.assertEqual(y.remainder_budget(0)['refined_remainder_upper'], 0)

    def test_new_estimates_never_call_old_window_numeric_helpers(self):
        with patch.object(y.previous, 'actual_budget', side_effect=AssertionError('old budget called')), \
             patch.object(y.previous, 'remainder_budget', side_effect=AssertionError('old remainder called')), \
             patch.object(y.geometry, 'parameters', side_effect=AssertionError('old local window called')):
            for t in (F(1, 200), F(1, 125), F(1, 100)):
                self.assertGreater(y.diagonal_budget(t)['Duhamel_denominator'], 0)
                self.assertGreater(y.remainder_budget(t)['refinement_denominator'], 0)
                self.assertGreater(y.actual_budget(t)['bare_full_conditional_and_coarse_gap_lower'], 0)


class ActualGapAndFluxTests(unittest.TestCase):
    def test_new_endpoint_has_fourfold_window_and_positive_actual_gap(self):
        row = y.endpoint_control()
        self.assertEqual(row['endpoint_window_factor'], 4)
        end = row['new_endpoint']
        self.assertEqual(end['bare_full_conditional_and_coarse_gap_lower'], F(16, 125))
        self.assertEqual(end['relative_memory_upper'], F(418609, 2073600))
        self.assertLess(end['relative_memory_upper'], F(1, 4))

    def test_same_old_endpoint_improves_both_certificates(self):
        row = y.endpoint_control()
        new, old = row['old_endpoint_with_new_bounds'], row['old_endpoint_YM92']
        self.assertGreater(new['bare_full_conditional_and_coarse_gap_lower'], old['bare_full_and_conditional_gap_lower'])
        self.assertLess(new['relative_memory_upper'], F(1, 8800))
        self.assertGreater(row['same_theta_memory_upper_improvement'], F(29, 5))

    def test_complete_flux_keeps_all_three_gradient_residual_terms(self):
        t = F(1, 100)
        row = y.actual_budget(t)
        expected = 520*t*t+2*(16*t+260*t*t)*(8*t/9)+F(8, 3)*t*(4*t)
        self.assertEqual(row['residual_gradient_over_sqrt_b'], expected)
        self.assertEqual(row['completed_flux_norm_over_sqrt_b'],
                         row['trial_flux_norm_over_sqrt_b']+expected/row['bare_full_conditional_and_coarse_gap_lower'])
        self.assertGreater(row['relative_memory_upper'], row['trial_flux_norm_over_sqrt_b']**2)

    def test_block_factor_cancels_in_relative_memory_on_new_window(self):
        rows = [y.actual_budget(F(1, 100), b) for b in (2, 3, 8, 64)]
        self.assertEqual(len({row['relative_memory_upper'] for row in rows}), 1)
        for row in rows:
            self.assertEqual(row['conditional_tensor_norm_upper'], row['factor']*row['relative_memory_upper'])

    def test_gap_and_memory_endpoint_directions_are_monotone(self):
        old_gap, old_eta = F(1, 2), F(0)
        for n in range(101):
            row = y.actual_budget(F(n, 10000))
            gap, eta = row['bare_full_conditional_and_coarse_gap_lower'], row['relative_memory_upper']
            self.assertLessEqual(gap, old_gap)
            self.assertGreaterEqual(eta, old_eta)
            self.assertLessEqual(eta, 4*row['mixed_Hessian_row_upper']**2/gap**2)
            old_gap, old_eta = gap, eta

    def test_frequency_bound_uses_two_separately_justified_bounds(self):
        for s in (0, 1, 10, 100):
            row = y.actual_budget(F(1, 100), 2, s)
            expected = min(row['relative_memory_upper'], row['source_energy_factor'] /
                           (row['bare_full_conditional_and_coarse_gap_lower']+s))
            self.assertEqual(row['frequency_memory_upper'], expected)
            self.assertEqual(row['relative_resolvent_error'], expected/(1-expected))
        self.assertLess(y.actual_budget(F(1, 100), 2, 100)['frequency_memory_upper'],
                        y.actual_budget(F(1, 100))['relative_memory_upper'])

    def test_zero_frequency_does_not_grant_unproved_extra_decay(self):
        self.assertTrue(y.previous.previous.previous.insufficient_zero_frequency_control())
        self.assertGreater(F(1, 101), F(1, 200))

    def test_free_case_has_exactly_zero_residual_and_memory(self):
        row = y.actual_budget(0, 8)
        for key in ('log_gradient_upper', 'mixed_remainder_upper', 'relative_memory_upper',
                    'residual_gradient_over_sqrt_b', 'relative_resolvent_error'):
            self.assertEqual(row[key], 0)
        self.assertEqual(row['bare_full_conditional_and_coarse_gap_lower'], F(1, 2))

    def test_scope_and_parameters_are_enforced(self):
        for t in (-1, F(1, 99)):
            with self.assertRaises(ValueError):
                y.actual_budget(t)
            with self.assertRaises(ValueError):
                y.diagonal_budget(t)
        for value in (True, 0.01):
            with self.assertRaises(TypeError):
                y.actual_budget(value)
        for b in (1, True, F(2)):
            with self.assertRaises(ValueError):
                y.actual_budget(0, b)
        with self.assertRaises(ValueError):
            y.actual_budget(0, 2, -1)


class EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fresh = y.run()

    def test_fresh_result_and_pin(self):
        self.assertEqual(y.check(self.fresh), y.PIN.read_text().strip())
        self.assertEqual(len(self.fresh['input_sha256']), 439)

    def test_modified_scope_fails_even_with_new_matching_hash(self):
        changed = copy.deepcopy(self.fresh)
        changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(changed))
            pin.write_text(y.p.canonical_sha(changed))
            with self.assertRaisesRegex(ValueError, 'scope'):
                y.check(changed, result, pin)

    def test_changed_result_or_pin_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(self.fresh))
            pin.write_text('0'*64)
            with self.assertRaisesRegex(ValueError, 'result/pin'):
                y.check(self.fresh, result, pin)
            changed = copy.deepcopy(self.fresh)
            changed['controls']['endpoint']['endpoint_window_factor'] = '5'
            result.write_text(json.dumps(changed))
            pin.write_text(y.p.canonical_sha(self.fresh))
            with self.assertRaisesRegex(ValueError, 'result/pin'):
                y.check(self.fresh, result, pin)

    def test_changed_parent_or_manifest_input_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'manifest'
            wrong_base = copy.deepcopy(manifest)
            wrong_base['base_commit'] = '0'*40
            path.write_text(json.dumps(wrong_base))
            with patch.object(y, 'SOURCES', path), self.assertRaisesRegex(ValueError, 'base'):
                y.source_checks()
            key = next(iter(manifest['upstream_sha256']))
            manifest['upstream_sha256'][key] = '0'*64
            path.write_text(json.dumps(manifest))
            with patch.object(y, 'SOURCES', path), self.assertRaisesRegex(ValueError, 'Changed upstream'):
                y.source_checks()

    def test_runtime_and_continuum_scope_remain_explicit(self):
        self.assertFalse(y.MODEL['older_numeric_helpers_used_outside_their_windows'])
        self.assertFalse(y.MODEL['weak_coupling_continuum_trajectory_controlled'])
        self.assertFalse(y.MODEL['positive_four_dimensional_continuum_mass_gap'])
        with patch.object(y.sys, 'version_info', (3, 11, 0)), self.assertRaises(RuntimeError):
            y.run()


if __name__ == '__main__':
    unittest.main()
