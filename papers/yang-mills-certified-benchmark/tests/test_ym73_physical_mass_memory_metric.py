import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym73_physical_mass_memory_metric as z


class PhysicalMassMemoryMetricTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_246_upstream_pins_and_four_local_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 246)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], '37d20649978324dd9dd269ccfdb425f189e7931a')
        for name in ('YM72_GROWING_BLOCK_MATRIX_MEMORY.md', 'YM64_VACUUM_ADAPTED_TOWER.md'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+name, manifest['upstream_sha256'])

    def test_upstream_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_continuum_overclaim_rejected(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['continuum_R4_QFT'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_no_binary_floats_or_unintended_control_characters(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_native_free_physical_gap_normalization(self):
        control = z.native_controls()
        self.assertEqual(control['native_plaquette']['free_electric_rate'], 3)
        self.assertEqual(control['native_plaquette']['variance'], F(1, 4))
        self.assertEqual(control['native_plaquette']['energy'], F(3, 4))
        self.assertEqual([row['shortest_cycle'] for row in control['geometry']], [4, 4])

    def test_exact_matrix_and_shape_domains(self):
        for bad in ([], [[1, 2]], [[1, 2], [3, 4]], [[1.0]], [[True]]):
            with self.assertRaises(ValueError): z.symmetric(bad)
        with self.assertRaises(ValueError): z.block_data([[1]], [[0, 0]], [[1]])
        for hidden in ([[0]], [[-1]], [[1, 0], [0, 0]]):
            with self.assertRaises(ValueError): z.positive_definite(hidden)
        with self.assertRaises(ValueError): z.block_data([[1]], [[2]], [[1]])

    def test_pencil_interval_independent_known_generalized_eigenvalues(self):
        S, Z = z.a.diag((6, 10)), z.a.diag((2, 5))
        row = z.pencil_interval(S, Z, F(1, 10**10))
        self.assertLessEqual(row['lower'], 2)
        self.assertGreaterEqual(row['upper'], 2)
        self.assertLessEqual(row['upper']-row['lower'], F(1, 10**10))
        # An invertible nonorthogonal congruence preserves this pencil spectrum.
        T = z.matrix(((1, 2), (0, 1)))
        congruence = lambda A: z.a.multiply(z.a.multiply(z.a.transpose(T), A), T)
        changed = z.pencil_interval(congruence(S), congruence(Z))
        self.assertLessEqual(changed['lower'], 2)
        self.assertGreaterEqual(changed['upper'], 2)

    def test_pencil_and_bracket_refuse_invalid_claims(self):
        for S, Z, tolerance in (([[1]], [[0]], 1), ([[-1]], [[1]], 1),
                                 ([[1]], [[1]], 0), ([[1]], [[1]], 0.1)):
            with self.assertRaises(ValueError): z.pencil_interval(S, Z, tolerance)
        for lo, hi, hidden in ((2, 1, 1), (-1, 1, 1), (0, 1, 0), (0.0, 1, 1)):
            with self.assertRaises(ValueError): z.gap_bracket(lo, hi, hidden)

    def test_harmonic_extension_changes_norm_and_preserves_static_energy(self):
        cell = z.fixture(4)
        for vector in ([F(1), F(0)], [F(2), -F(3, 7)]):
            extended = z.a.apply(cell['harmonic_extension'], vector)
            self.assertEqual(z.a.dot(extended, extended), z.a.dot(vector, z.a.apply(cell['Z'], vector)))
            self.assertEqual(z.a.dot(extended, z.a.apply(cell['H'], extended)),
                             z.a.dot(vector, z.a.apply(cell['S0'], vector)))
            self.assertGreater(z.a.dot(extended, extended), z.a.dot(vector, vector))

    def test_gap_bound_agrees_with_independent_full_matrix_inertia(self):
        for row in z.matrix_controls()['gap_convergence']:
            cell = z.fixture(row['clock'])
            bound = row['full_gap_bracket']
            physical = z.a.scale(cell['H'], row['clock'])
            self.assertTrue(z.a.psd(z.a.add(physical, z.a.identity(4), -bound['lower'])))
            self.assertFalse(z.a.psd(z.a.add(physical, z.a.identity(4), -(bound['upper']+F(1, 100)))))
            self.assertGreaterEqual(row['independently_bisected_full_gap']['lower'], bound['lower'])
            self.assertLessEqual(row['independently_bisected_full_gap']['upper'], bound['upper'])

    def test_hidden_floor_must_be_independently_true(self):
        with self.assertRaises(ValueError): z.gap_certificate([[3]], [[0]], [[F(1, 100)]], 1)
        with self.assertRaises(ValueError): z.resolvent_control([[3]], [[0]], [[F(1, 100)]], 1, 1, 1)
        honest = z.gap_certificate([[3]], [[0]], [[F(1, 100)]], F(1, 100))
        self.assertLessEqual(honest['full_gap']['lower'], F(1, 100))

    def test_sharp_lower_boundary_is_approached_from_above(self):
        rows = z.sharpness_controls()
        self.assertEqual(rows[0]['lower_bound'], F(6, 7))
        self.assertLess(rows[-1]['distance_upper'], F(1, 10000))
        self.assertGreater(rows[0]['distance_upper'], rows[-1]['distance_upper'])
        # With no coupling, the upper endpoint can be exact.
        cell = z.gap_certificate([[F(3, 2)]], [[0]], [[2]], 2)
        self.assertEqual(cell['full_gap']['upper'], F(3, 2))

    def test_zero_static_form_does_not_certify_positive_gap(self):
        cell = z.block_data([[1]], [[1]], [[1]])
        self.assertEqual(cell['S0'], [[0]])
        row = z.gap_certificate([[1]], [[1]], [[1]], 1)
        self.assertEqual(row['full_gap']['lower'], 0)
        self.assertEqual(row['full_gap']['upper'], 0)

    def test_noncommuting_fixture_does_not_use_simultaneous_diagonalization(self):
        for C in z.matrix_controls()['nonzero_static_metric_commutators']:
            self.assertNotEqual(C[0][1], 0)
            self.assertEqual(C[0][1], -C[1][0])

    def test_resolvent_difference_is_positive_with_exact_congruence_bound(self):
        for row in z.matrix_controls()['exact_resolvent_checks']:
            difference = row['positive_difference']
            self.assertTrue(z.a.psd(difference))
            self.assertTrue(z.a.psd(z.a.add(row['unweighted_difference_upper'], difference, -1)))
            self.assertEqual(row['normalized_norm_error_upper'], 1/(2*row['clock']))

    def test_resolvent_domains_and_high_spectral_parameter(self):
        cell = z.fixture(4)
        for hidden, clock, spectral in ((0, 1, 1), (2, 0, 1), (2, 1, 0), (2, 1, 0.1)):
            with self.assertRaises(ValueError):
                z.resolvent_control(cell['C'], cell['B'], cell['D'], hidden, clock, spectral)
        row = z.resolvent_control(cell['C'], cell['B'], cell['D'], 2, 4, F(10**12))
        self.assertTrue(z.a.psd(row['positive_difference']))

    def test_interval_pencil_uses_opposite_endpoints_for_memory(self):
        cert = z.interval_pencil([[5]], [[F(51, 10)]], [[F(19, 10)]], [[2]],
                                [[F(9, 10)]], [[1]], 1)
        self.assertLessEqual(cert['effective_lower'], F(3, 2))
        self.assertGreater(cert['effective_lower'], F(3, 2)-F(1, 10**10))
        self.assertEqual(cert['effective_upper'], F(32, 19))
        self.assertGreaterEqual(cert['full_gap']['lower'], F(3, 5)-F(1, 10**10))

    def test_interval_pencil_rejects_reversal_inconclusiveness_and_zero_probe(self):
        valid = ([[5]], [[6]], [[1]], [[2]], [[0]], [[1]])
        with self.assertRaises(ValueError): z.interval_pencil(*valid, 1, probe=[0])
        with self.assertRaises(ValueError): z.interval_pencil([[6]], [[5]], *valid[2:], 1)
        with self.assertRaises(ValueError): z.interval_pencil([[1]], [[2]], [[2]], [[3]], [[0]], [[1]], 1)
        with self.assertRaises(ValueError): z.interval_pencil(*valid, 1, probe=[1, 0])

    def test_full_uncertain_pencils_enclose_known_effective_energy(self):
        for row in z.interval_controls():
            self.assertLessEqual(row['effective_lower'], F(3, 2))
            self.assertGreaterEqual(row['effective_upper'], F(3, 2))
            self.assertLess(row['effective_upper']-row['effective_lower'], F(1, 1000))

    def test_clock_scaling_leaves_first_moment_and_metric_invariant(self):
        cell = z.fixture(4)
        clock = F(7, 3)
        changed = z.block_data(*(z.a.scale(cell[key], clock) for key in ('C', 'B', 'D')))
        self.assertEqual(changed['M0'], z.a.scale(cell['M0'], clock))
        self.assertEqual(changed['M1'], cell['M1'])
        self.assertEqual(changed['Z'], cell['Z'])
        self.assertEqual(changed['S0'], z.a.scale(cell['S0'], clock))

    def test_missing_first_moment_gives_wrong_limiting_mass(self):
        rows = z.scope_counterexamples()['missing_first_moment']
        for row in rows:
            self.assertEqual(row['unnormalized_static_energy']/row['correct_effective_energy'], 5)
            self.assertLess(row['actual_gap_interval']['upper'], row['correct_effective_energy'])
        self.assertGreater(rows[-1]['actual_gap_interval']['lower'], F(149, 100))

    def test_entire_observed_memory_tower_can_miss_an_uncoupled_sector(self):
        rows = z.scope_counterexamples()['invisible_sector']
        for row in rows:
            self.assertTrue(row['all_memory_moments_zero'])
            self.assertEqual(row['retained_response_at_one'], F(1, 4))
        self.assertNotEqual(rows[0]['hidden_gap'], rows[-1]['hidden_gap'])

    def test_actual_slow_probe_weight_and_enlarged_gram_are_nontrivial(self):
        row = z.slow_probe_budget(F(1, 10**8*16), 16)
        self.assertEqual(row['microscopic_low_weight_lower'], F(89, 114))
        self.assertGreater(row['enlarged_Gram_lower'], F(53, 100))
        self.assertLess(row['enlarged_Gram_upper'], F(147, 100))
        self.assertTrue(row['independent_slow_direction_certified'])
        small = z.slow_probe_budget(0, 3)
        self.assertEqual(small['enlarged_Gram_lower'], 0)
        self.assertFalse(small['independent_slow_direction_certified'])

    def test_slow_probe_count_budget_and_limit_improve(self):
        rows = z.slow_and_scale_controls()['actual_model_slow_probe_enclosures']
        for i in range(len(rows)-1):
            self.assertLess(rows[i]['microscopic_low_weight_lower'], rows[i+1]['microscopic_low_weight_lower'])
            self.assertGreater(rows[i]['micro_large_overlap_upper'], rows[i+1]['micro_large_overlap_upper'])
        with self.assertRaises(ValueError): z.slow_probe_budget(F(1, 1000), 16)

    def test_physical_calibration_forces_spectrum_escape_in_current_window(self):
        row = z.physical_escape_budget(F(1, 16), 4)
        self.assertEqual(row['physical_gap_lower'], 2560)
        self.assertEqual(row['physical_gap_lower_squared'], 2560**2)
        self.assertEqual(row['calibrated_clock_lower_squared'], (125*2560)**2)
        scaled = z.physical_escape_budget(F(1, 16), 4, 1, F(9), F(4))
        self.assertEqual(scaled['physical_gap_lower'], 6*row['physical_gap_lower'])

    def test_square_root_rounds_down_and_small_floor_stays_safe(self):
        for value in (0, F(1, 4), 2, F(1, 10**100)):
            low = z.sqrt_lower(value)
            self.assertLessEqual(low*low, value)
            self.assertGreater((low+F(1, z.GRID))**2, value)
        self.assertEqual(z.escape_heat_upper(F(10**100), 4, 1, 1), 1)
        for value in (-1, True, 0.5):
            with self.assertRaises(ValueError): z.sqrt_lower(value)

    def test_heat_escape_bound_includes_norm_and_tends_to_zero(self):
        rows = z.slow_and_scale_controls()['current_trajectory_spectrum_escape']
        self.assertLess(rows[-1]['time_one_correlation_upper'], F(1, 10**20))
        self.assertEqual(z.escape_heat_upper(F(1, 4), 4, 1, 1, 3),
                         3*z.escape_heat_upper(F(1, 4), 4, 1, 1))
        for spacing, factor, count in ((0, 4, 1), (1, 2, 1), (1, 4, 0), (0.1, 4, 1)):
            with self.assertRaises(ValueError): z.physical_escape_budget(spacing, factor, count)
        with self.assertRaises(ValueError): z.escape_heat_upper(1, 4, 1, 0)

    def test_finite_gap_does_not_make_a_selected_state_survive(self):
        rows = z.scope_counterexamples()['finite_gap_without_local_visibility']
        self.assertTrue(all(row['full_gap'] == 1 and row['selected_state_norm_squared'] == 1 for row in rows))
        self.assertGreater(rows[0]['correlation_at_log_two'], rows[-1]['correlation_at_log_two'])
        self.assertLess(rows[0]['selected_state_energy'], rows[-1]['selected_state_energy'])

    def test_scope_does_not_discharge_the_continuum_inputs(self):
        claim = self.cert['model']
        for field in ('actual_microscopic_slow_probe_enrichment', 'memory_normalization_Z_equals_I_plus_M1',
                      'noncircular_conditional_physical_gap_bracket',
                      'uniform_normalized_resolvent_error_one_over_hidden_physical_floor',
                      'actual_YM72_calibrated_trajectory_has_escaping_spectrum',
                      'continuum_gap_passage_on_a_dense_core_is_conditional'):
            self.assertTrue(claim[field])
        for field in ('independently_certified_hidden_floor_on_weak_coupling_trajectory',
                      'actual_weak_coupling_static_pencil_computed',
                      'full_low_energy_space_retained_by_one_extra_probe',
                      'finite_fixture_is_Yang_Mills_vacuum',
                      'lattice_gap_convergence_alone_constructs_continuum_fields',
                      'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(claim[field])


if __name__ == '__main__':
    unittest.main()
