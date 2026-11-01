import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym75_rh_source_resolved_transfer as z


class SourceResolvedMathematicsTests(unittest.TestCase):
    def test_exact_inertia_includes_zero_diagonal_two_by_two_pivots(self):
        self.assertEqual(z.inertia([[0, 2], [2, 0]]), dict(negative=1, zero=0, positive=1))
        self.assertEqual(z.inertia([[0, 0], [0, 0]]), dict(negative=0, zero=2, positive=0))
        self.assertEqual(z.inertia([[2, 2], [2, 2]]), dict(negative=0, zero=1, positive=1))

    def test_inertia_is_invariant_under_independent_rational_congruences(self):
        diagonal = z.a.diag((-3, 0, 2, 7))
        for U in (z.matrix(((1, 2, 3, 1), (0, 1, -2, 3), (0, 0, 1, 5), (0, 0, 0, 1))),
                  z.matrix(((1, 0, 0, 1), (1, 1, 0, 0), (0, 1, 1, 0), (0, 0, 1, 2)))):
            z.a.inverse(U)
            changed = z.a.multiply(z.a.multiply(z.a.transpose(U), diagonal), U)
            self.assertEqual(z.inertia(changed), dict(negative=1, zero=1, positive=2))

    def test_exact_domains_reject_floats_boolean_and_nonsymmetric_inputs(self):
        for bad in ([], [[1, 2]], [[1, 2], [0, 1]], [[1.0]], [[True]]):
            with self.assertRaises(ValueError): z.inertia(bad)
        for theta in (True, 0.1):
            with self.assertRaises(ValueError): z.cube_parameters(theta)

    def test_full_matrix_spectra_validate_source_resolved_inertia(self):
        rows = z.matrix_controls()['full_matrix_inertia_controls']
        self.assertEqual(len(rows), 18)
        for row in rows:
            self.assertLessEqual(row['exact_full_inertia']['negative'], row['inertia']['negative'])
            self.assertTrue(z.a.psd(z.a.add(row['exact_Schur'], row['lower_Schur_matrix'], -1)))
            if row['lower_eigenvalue_certified']:
                self.assertLessEqual(row['exact_full_inertia']['negative'], row['eigenvalue_index'])

    def test_false_hidden_floor_and_source_envelope_are_rejected(self):
        with self.assertRaises(ValueError): z.finite_block_control([[3]], [[1]], [[2]], 3, 0)
        with self.assertRaises(ValueError): z.finite_block_control([[3]], [[1]], [[2]], 2, 0, envelope=[[F(1, 2)]])
        with self.assertRaises(ValueError): z.finite_block_control([[3]], [[1, 2]], [[2]], 2, 0)
        honest = z.finite_block_control([[3]], [[1]], [[2]], 2, 0, envelope=[[2]])
        self.assertTrue(honest['lower_eigenvalue_certified'])

    def test_equal_source_norm_does_not_determine_orientation_or_gap_count(self):
        left, right = z.matrix_controls()['equal_norm_different_source_orientation']
        for row in (left, right):
            gram = row['actual_source_Gram']
            self.assertEqual(max(gram[i][i] for i in range(2)), F(1, 25))
        self.assertEqual(left['inertia']['negative'], 1)
        self.assertEqual(right['inertia']['negative'], 2)
        self.assertTrue(left['lower_eigenvalue_certified'])
        self.assertFalse(right['lower_eigenvalue_certified'])
        self.assertNotEqual(left['actual_source_Gram'], right['actual_source_Gram'])

    def test_zero_coupling_does_not_allow_dropping_low_hidden_modes(self):
        row = z.matrix_controls()['invisible_uncoupled_hidden_sector']
        self.assertEqual(row['exact_full_inertia'], dict(negative=0, zero=0, positive=3))
        with self.assertRaises(ValueError):
            z.source_resolved_schur([[3, 0], [0, 5]], [[0, 0], [0, 0]], 1, F(3, 2), 0)
        # The retained matrix alone would incorrectly certify the full ground >=3.
        self.assertEqual(z.inertia([[F(3, 2), 0], [0, F(7, 2)]])['negative'], 0)

    def test_Schur_threshold_requires_a_valid_index_and_subhidden_energy(self):
        for hidden, energy, index in ((2, 2, 0), (2, 3, 0), (2, 0, 1), (2, 0, -1), (2, 0, True)):
            with self.assertRaises(ValueError): z.source_resolved_schur([[3]], [[1]], hidden, energy, index)
        with self.assertRaises(ValueError): z.source_resolved_schur([[3]], [[-1]], 2, 0)
        with self.assertRaises(ValueError): z.source_resolved_schur([[3]], [[1, 0], [0, 1]], 2, 0)

    def test_scalar_residual_identity_matches_direct_independent_arithmetic(self):
        d, f, x, shift = F(7), F(3), F(1, 5), F(2)
        row = z.finite_residual_memory([[d]], [[f]], [[x]], 7, shift=shift)
        residual = f-(d-shift)*x
        expected_lower = 2*f*x-(d-shift)*x*x
        self.assertEqual(row['lower'], [[expected_lower]])
        self.assertEqual(row['exact_memory'], [[f*f/(d-shift)]])
        self.assertEqual(row['exact_error'], [[residual*residual/(d-shift)]])
        self.assertEqual(row['exact_memory'][0][0]-row['lower'][0][0], row['exact_error'][0][0])

    def test_residual_upper_bound_is_sharp_when_hidden_operator_equals_floor(self):
        D, Fsource, X = z.a.scale(z.a.identity(2), 3), [[1, 2], [3, -1]], [[0, 1], [F(1, 2), 0]]
        row = z.finite_residual_memory(D, Fsource, X, 3, shift=F(1, 2))
        self.assertEqual(row['upper'], row['exact_memory'])
        self.assertTrue(z.a.psd(z.a.add(row['exact_memory'], row['lower'], -1)))

    def test_memory_error_is_quadratic_in_solve_error_and_zero_for_exact_solve(self):
        rows = z.residual_controls()['quadratic_residual_examples']
        self.assertEqual(len(rows), 12)
        for start in (0, 4, 8):
            group = rows[start:start+4]
            self.assertEqual(group[1]['exact_error'], z.a.scale(group[0]['exact_error'], F(1, 4)))
            self.assertEqual(group[2]['exact_error'], z.a.scale(group[1]['exact_error'], F(1, 4)))
            self.assertEqual(group[3]['lower'], group[3]['upper'])
            self.assertEqual(group[3]['lower'], group[3]['exact_memory'])

    def test_symmetrized_raw_source_Gram_can_overestimate_exact_memory(self):
        controls = z.residual_controls()
        row = controls['symmetrized_raw_Gram_not_a_lower_bound']
        self.assertEqual(controls['raw_Gram_value_for_counterexample'], [[1]])
        self.assertEqual(row['exact_memory'], [[F(1, 2)]])
        self.assertEqual(row['lower'], [[0]])
        self.assertGreater(controls['raw_Gram_value_for_counterexample'][0][0], row['exact_memory'][0][0])

    def test_residual_premises_and_positive_residual_Gram_are_enforced(self):
        with self.assertRaises(ValueError): z.finite_residual_memory([[2]], [[1]], [[0]], 3)
        with self.assertRaises(ValueError): z.finite_residual_memory([[2]], [[1]], [[0]], 2, shift=2)
        with self.assertRaises(ValueError): z.finite_residual_memory([[2]], [[1], [2]], [[0]], 2)
        with self.assertRaises(ValueError): z.residual_enclosure([[1]], [[-1]], 2)
        with self.assertRaises(ValueError): z.residual_enclosure([[1]], [[1]], 0)
        with self.assertRaises(ValueError): z.residual_enclosure([[1]], [[1]], 2, clock=0)

    def test_nondiagonal_memory_intervals_are_valid_in_Loewner_order(self):
        rows = z.residual_controls()['quadratic_residual_examples']
        nonzero_off_diagonal = False
        for row in rows:
            self.assertTrue(z.a.psd(z.a.add(row['exact_memory'], row['lower'], -1)))
            self.assertTrue(z.a.psd(z.a.add(row['upper'], row['exact_memory'], -1)))
            self.assertEqual(row['physical_norm_error_upper'],
                             3*sum(row['error_Gram_upper'][i][i] for i in range(2)))
            nonzero_off_diagonal |= bool(row['exact_error'][0][1])
        self.assertTrue(nonzero_off_diagonal)

    def test_band_decoder_bound_is_attained_by_independent_orthonormal_states(self):
        row = z.decoder_controls()
        self.assertEqual(z.a.multiply(z.a.transpose(row['band_isometry']), row['band_isometry']), z.a.identity(2))
        self.assertEqual(z.a.add(row['retained_Gram'], row['hidden_Gram']), z.a.identity(2))
        self.assertEqual(row['exact_band_decoder'], [[F(5, 3), 0], [0, F(5, 4)]])
        self.assertEqual(row['sharp_budget']['inverse_decoder_norm_squared_upper'], F(25, 9))

    def test_injectivity_does_not_supply_a_uniform_decoder_bound(self):
        rows = z.decoder_controls()['injective_without_uniform_stability']
        self.assertTrue(all(row['retained_scalar'] > 0 for row in rows))
        self.assertEqual([row['inverse_decoder_norm_squared_upper'] for row in rows], [4, 16, 64, 256])
        self.assertLess(rows[-1]['retained_Gram_floor'], rows[0]['retained_Gram_floor'])

    def test_decoder_requires_strict_leakage_margin_and_no_state_existence_claim(self):
        for bad in (-1, 1, 2, 0.5, True):
            with self.assertRaises(ValueError): z.decoder_budget(bad)
        row = z.decoder_budget(0)
        self.assertEqual(row['inverse_decoder_norm_squared_upper'], 1)
        self.assertFalse(row['nonvacuum_state_existence_asserted'])

    def test_actual_cube_structured_interval_improves_both_YM74_endpoints(self):
        row = z.actual_cube_controls()
        self.assertEqual(row['certified_bare_gap']['lower'], F(49069036, 10**8))
        self.assertEqual(row['certified_bare_gap']['upper'], F(300712209, 10**8))
        self.assertGreater(row['certified_bare_gap']['lower'], row['previous_scalarized_gap']['lower'])
        self.assertLess(row['certified_bare_gap']['upper'], row['previous_scalarized_gap']['upper'])
        self.assertEqual(row['ground_threshold']['inertia'], dict(negative=0, zero=0, positive=241))
        self.assertEqual(row['first_excited_threshold']['inertia'], dict(negative=1, zero=0, positive=240))

    def test_actual_cube_inertia_counts_transverse_and_winding_sectors(self):
        row = z.cube_schur_inertia(F(1, 48), F(499, 1000))
        self.assertEqual(row['inertia']['negative'], 241)
        self.assertEqual(row['uniform_inertia']['negative'], 2)
        self.assertLess(row['transverse_pivot'], 0)
        self.assertLess(row['winding_pivot'], 0)
        self.assertEqual(row['transverse_multiplicity'], 191)
        self.assertEqual(row['winding_multiplicity'], 48)

    def test_actual_cube_zero_threshold_free_limit_and_sign_symmetry(self):
        self.assertEqual(z.cube_schur_inertia(0, 0)['inertia'], dict(negative=0, zero=1, positive=240))
        self.assertEqual(z.cube_schur_inertia(0, 3)['inertia'], dict(negative=1, zero=240, positive=0))
        self.assertEqual(z.cube_schur_inertia(F(1, 48), 0)['inertia'], dict(negative=1, zero=0, positive=240))
        positive = z.cube_schur_inertia(F(1, 48), F(12, 25))
        negative = z.cube_schur_inertia(F(-1, 48), F(12, 25))
        self.assertEqual(positive['inertia'], negative['inertia'])
        self.assertEqual(positive['uniform_pivots'], negative['uniform_pivots'])
        with self.assertRaises(ValueError): z.cube_schur_inertia(F(1, 48), F(1, 2))

    def test_comparison_root_intervals_agree_with_independent_exact_cubic(self):
        roots = z.actual_cube_controls()['comparison_root_brackets']
        polynomial = lambda x: 9216*x**3-32256*x**2+13249*x+96
        self.assertLessEqual(polynomial(roots[0]['certified_full_eigenvalue_lower']), 0)
        self.assertGreaterEqual(polynomial(roots[0]['comparison_threshold_upper']), 0)
        self.assertGreaterEqual(polynomial(roots[1]['certified_full_eigenvalue_lower']), 0)
        self.assertLessEqual(polynomial(roots[1]['comparison_threshold_upper']), 0)
        self.assertTrue(all(row['width'] <= F(1, 10**12) for row in roots))
        with self.assertRaises(ValueError): z.cube_threshold_bisection(F(1, 48), 0, 0, F(1, 10))
        with self.assertRaises(ValueError): z.cube_threshold_bisection(F(1, 48), 1, 0, F(1, 10))

    def test_actual_YM74_coverage_supplies_a_small_decoder_price(self):
        row = z.decoder_controls()['actual_YM74_cutoff_decoder_bound']
        self.assertLess(row['leakage_squared'], F(3, 128))
        self.assertLess(row['inverse_decoder_norm_squared_upper'], F(128, 125))
        self.assertFalse(row['nonvacuum_state_existence_asserted'])


class SourceResolvedCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_260_upstream_pins_and_four_local_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 260)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertIn('papers/yang-mills-certified-benchmark/YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES.md', manifest['upstream_sha256'])

    def test_upstream_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_uniform_or_continuum_overclaims_are_rejected(self):
        for field in ('uniform_weak_coupling_positive_gap', 'continuum_R4_QFT'):
            wrong = copy.deepcopy(self.cert)
            wrong['model'][field] = True
            with tempfile.TemporaryDirectory() as directory:
                result, pin = Path(directory)/'result', Path(directory)/'pin'
                result.write_text(json.dumps(wrong))
                pin.write_text(z.p.canonical_sha(wrong))
                with self.assertRaises(ValueError): z.check(wrong, result, pin)
                with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_certificate_has_no_binary_floats_or_unintended_control_characters(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_scope_retains_all_actual_model_and_continuum_boundaries(self):
        claims = self.cert['model']
        for field in ('source_resolved_matrix_Schur_lower_bound', 'exact_full_retained_inertia_not_scalarized',
                      'quadratic_residual_Loewner_memory_enclosures',
                      'actual_four_cube_theta_one_over_48_improvement', 'all_241_physical_retained_directions_counted'):
            self.assertTrue(claims[field])
        for field in ('hidden_auxiliary_modes_may_be_discarded', 'source_injectivity_alone_implies_uniform_gap',
                      'low_energy_coverage_proves_nonvacuum_state_exists', 'uniform_weak_coupling_positive_gap',
                      'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(claims[field])


if __name__ == '__main__':
    unittest.main()
