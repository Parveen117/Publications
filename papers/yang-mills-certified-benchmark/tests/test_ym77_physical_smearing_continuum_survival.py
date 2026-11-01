import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym77_physical_smearing_continuum_survival as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class PhysicalSurvivalMathematicsTests(unittest.TestCase):
    def test_seven_channels_are_needed_even_for_exact_unit_quaternions(self):
        row = z.trace_channels(z.model.UNITS[1:])
        self.assertEqual(row['singleton_squares'], 0)
        self.assertEqual(row['pair_squares'], 0)
        self.assertEqual(row['triple_square'], 1)
        self.assertEqual(row['imaginary_Gram'], z.a.identity(3))
        self.assertEqual(row['imaginary_Gram_determinant'], 1)

    def test_quaternion_Gram_identity_includes_nonzero_real_and_triple_parts(self):
        palette = ((F(1, 3), F(2, 3), F(2, 3), 0),
                   (F(3, 5), 0, F(4, 5), 0), (F(5, 13), 0, 0, F(12, 13)))
        for qs in product(palette, repeat=3):
            row = z.trace_channels(qs)
            self.assertEqual(row['imaginary_Gram_determinant'], row['signed_triple_product']**2)
            self.assertGreater(row['total_trace_squares'], F(1, 5))
            self.assertEqual(sum(row['imaginary_Gram'][i][i] for i in range(3)), 3-row['singleton_squares'])

    def test_written_trace_contradiction_constants_are_strictly_separated(self):
        self.assertLess(F(2, 5), F(4, 9))
        self.assertEqual((F(14, 5)-F(2, 3))/9, F(32, 135))
        self.assertGreater(F(32, 135), F(1, 5))
        # M=aa*-2diag(a_i^2) is bounded in norm by A, including singular cases.
        for vector in ((1, 0, 0), (1, 1, 1), (F(1, 3), F(2, 5), F(-4, 7))):
            A = sum(x*x for x in vector)
            M = [[vector[i]*vector[j]-(2*vector[i]*vector[i] if i == j else 0)
                  for j in range(3)] for i in range(3)]
            self.assertTrue(z.a.psd(z.a.add(z.a.scale(z.a.identity(3), A), M)))
            self.assertTrue(z.a.psd(z.a.add(z.a.scale(z.a.identity(3), A), M, -1)))

    def test_winding_graph_has_exact_center_symmetry_and_distinct_characters(self):
        for side in (4, 5):
            row = z.center_geometry(side)
            self.assertEqual(row['rooted_axial_loop_lengths'], [side]*3)
            self.assertEqual(row['center_seam_sizes'], [side*side]*3)
            self.assertEqual(row['sum_channel_lengths'], 12*side)
            self.assertEqual(row['seven_character_Gram'], z.a.scale(z.a.identity(7), 8))
            self.assertEqual(len({tuple(c['center_parities']) for c in row['channel_controls']}), 7)
            self.assertTrue(row['all_plaquettes_have_even_seam_parity'])
            self.assertTrue(row['edge_disjoint_rooted_cycles'])

    def test_native_long_word_Dirichlet_density_has_perimeter_factor(self):
        for side in (4, 5):
            for row in z.center_geometry(side)['channel_controls']:
                self.assertEqual(row['length'], side*len(row['subset']))
                self.assertEqual(row['exact_Dirichlet_density'], F(row['length'], 4)*(1-row['trace']**2))

    def test_word_derivatives_match_native_inverse_link_polynomial(self):
        qs = ((F(3, 5), F(4, 5), 0, 0), (F(5, 13), 0, F(12, 13), 0),
              (0, 0, 0, 1), (0, 1, 0, 0))
        word = ((3, -1), (0, 1), (2, -1), (1, 1))
        polynomial = z.model.loop_polynomial(word)
        value, gradients = z.word_gradient(qs, word)
        self.assertEqual(value, z.model.evaluate(polynomial, qs))
        for edge in range(4):
            for axis in (1, 2, 3):
                self.assertEqual(gradients[edge][axis-1],
                                 z.model.evaluate(z.native.derivative(polynomial, edge, axis), qs))

    def test_center_energy_budget_gives_a_true_scope_limited_global_upper_bound(self):
        row = z.center_state_budget(8, F(1, 100), threshold=2, time=3)
        self.assertEqual(row['state_norm_squared'], dict(lower=F(1, 5), upper=F(7)))
        self.assertEqual(row['physical_form_energy_upper'], F(6, 25))
        self.assertEqual(row['full_physical_gap_upper'], F(6, 5))
        self.assertEqual(row['normalized_mass_above_threshold_upper'], F(3, 5))
        self.assertEqual(row['unnormalized_heat_loss_upper'], F(18, 25))
        self.assertTrue(row['includes_nontrivial_global_center_sectors'])
        self.assertFalse(row['bound_for_center_neutral_or_local_QFT_gap'])

    def test_actual_bad_trajectory_closes_global_gap_while_local_UV_band_escapes(self):
        for n in (2, 4, 32, 64, 1000):
            row = z.closing_trajectory(n)
            self.assertEqual(row['physical_extent'], n)
            self.assertEqual(row['theta'], n**8)
            self.assertEqual(row['physical_clock']*n**4, n)
            self.assertEqual(row['actual_center_state']['full_physical_gap_upper'], F(15, n))
            self.assertEqual(row['actual_local_UV_band']['lower'], F(n, 45))
            self.assertEqual(row['actual_plaquette_deficiency_upper'], F(9, 4*n**4))
            self.assertFalse(row['large_lattice_geometry_allocated'])
            self.assertFalse(row['selected_asymptotically_free_trajectory'])

    def test_geometry_and_center_budgets_reject_invalid_domains(self):
        for bad in (1, 3, True, 4.0):
            with self.assertRaises(ValueError): z.center_state_budget(bad)
        with self.assertRaises(ValueError): z.trace_channels((z.model.y.ONE,)*2)
        with self.assertRaises(ValueError): z.word_gradient((z.model.y.ONE,)*2, ((0, 1), (0, -1)))
        with self.assertRaises(ValueError): z.center_state_budget(4, 0)
        with self.assertRaises(ValueError): z.center_state_budget(4, threshold=0)
        with self.assertRaises(ValueError): z.closing_trajectory(1)

    def test_exact_cell_averages_reproduce_independent_linear_Jensen_error(self):
        for n in (1, 2, 4, 7):
            row = z.cell_jensen({(1, 0, 0): 1}, n)
            self.assertEqual(row['continuum_L2_squared'], F(1, 3))
            self.assertEqual(row['exact_Jensen_defect'], F(1, 12*n*n))
            self.assertEqual(row['cell_count'], n**3)
        constant = z.cell_jensen({(0, 0, 0): 3}, 3)
        self.assertEqual(constant['exact_Jensen_defect'], 0)
        self.assertEqual(constant['continuum_L2_squared'], 9)

    def test_signed_multivariate_cell_averages_obey_Jensen_and_refine(self):
        poly = {(1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): -1}
        coarse, fine = z.cell_jensen(poly, 2), z.cell_jensen(poly, 4)
        self.assertEqual(coarse['continuum_L2_squared'], F(3, 2))
        self.assertLess(fine['exact_Jensen_defect'], coarse['exact_Jensen_defect'])
        quadratic = z.cell_jensen({(2, 0, 0): 1}, 3)
        self.assertEqual(quadratic['continuum_L2_squared'], F(1, 5))
        self.assertGreater(quadratic['exact_Jensen_defect'], 0)
        self.assertEqual(z.cell_jensen({}, 2)['continuum_L2_squared'], 0)

    def test_physical_smearing_cancels_theta_and_one_inverse_spacing_exactly(self):
        cell = z.cell_jensen({(1, 0, 0): 1}, 4)
        for theta in (1, 25, 26, 10**12):
            row = z.physical_smear_budget(theta, cell['spacing'], cell['cell_average_square_sum'],
                                         cell['continuum_L2_squared'], 3)
            self.assertEqual(row['smear_coefficient_fourth_power'], F(theta, 4**8))
            self.assertEqual(row['uniform_physical_carre_du_champ_upper'], 4)
            self.assertEqual(row['exact_physical_discrete_Dirichlet_upper'], 12*cell['volume_weighted_square_sum'])
            self.assertLessEqual(row['exact_physical_discrete_Dirichlet_upper'], 4)
            self.assertFalse(row['nonzero_variance_lower_bound_provided'])

    def test_smear_and_polynomial_domain_errors_do_not_hide_false_Jensen_premises(self):
        with self.assertRaises(ValueError): z.physical_smear_budget(25, 1, 2, 1)
        with self.assertRaises(ValueError): z.physical_smear_budget(25, 0, 0, 1)
        with self.assertRaises(ValueError): z.cell_jensen({(1, 0): 1}, 2)
        with self.assertRaises(ValueError): z.cell_jensen({(-1, 0, 0): 1}, 2)
        with self.assertRaises(ValueError): z.cell_jensen({(1, 0, 0): 1.0}, 2)
        with self.assertRaises(ValueError): z.integrate_box({}, (0, 0, 1), (1, 1, 0))

    def test_anchored_increment_moment_constants_and_time_scaling(self):
        row = z.increment_moment_budget(2, F(1, 4), 3)
        self.assertEqual(row['even_moment_upper'], 96)
        self.assertEqual(row['second_moment_upper'], 1)
        self.assertTrue(row['requires_uniform_pointwise_carre_du_champ_not_just_mean_energy'])
        smaller = z.increment_moment_budget(2, F(1, 8), 3)
        self.assertEqual(smaller['even_moment_upper'], row['even_moment_upper']/8)
        self.assertEqual(z.increment_moment_budget(1, 0, 2)['even_moment_upper'], 0)
        for order in (0, -1, True):
            with self.assertRaises(ValueError): z.increment_moment_budget(1, 1, order)

    def test_directed_filter_measures_enclose_independent_decimal_computation(self):
        atoms, weights, time = (F(0), F(1, 100), F(1), F(100)), (F(5), F(1), F(1, 2), F(1, 1000)), F(3, 2)
        row = z.finite_filter_control(atoms, weights, time, 1)
        with localcontext() as context:
            context.prec = 90
            mass = sum(decimal(w)*(1-(-decimal(time*x)).exp())**2 for x, w in zip(atoms, weights))
            energy = sum(decimal(w*x)*(1-(-decimal(time*x)).exp())**2 for x, w in zip(atoms, weights))
            inverse = sum(decimal(w/x)*(1-(-decimal(time*x)).exp())**2 for x, w in zip(atoms, weights) if x)
            for key, true in (('filtered_norm_squared', mass), ('filtered_energy', energy), ('filtered_inverse_energy', inverse)):
                self.assertLessEqual(decimal(row[key]['lower']), true)
                self.assertGreaterEqual(decimal(row[key]['upper']), true)
        self.assertEqual(row['filtered_atom_masses'][0], dict(lower=0, upper=0))

    def test_filter_bounds_control_norm_energy_inverse_energy_and_ordinary_increments(self):
        for time in (F(1, 10), F(1), F(3)):
            row = z.finite_filter_control((F(1, 1000), 1, 1000), (2, F(1, 2), F(1, 10000)), time, 1)
            self.assertLessEqual(row['filtered_norm_squared']['upper'], time)
            self.assertLessEqual(row['filtered_energy']['upper'], 1)
            self.assertLessEqual(row['filtered_inverse_energy']['upper'], time*time)
            self.assertLessEqual(row['ordinary_increment_variance']['upper'], 2*time)
        self.assertFalse(z.filter_budget(1, 1)['uniform_positive_filtered_norm_inferred'])

    def test_two_sided_survival_budget_has_a_nonempty_finite_physical_band(self):
        row = z.filter_survival_budget(2, 3, F(1, 2))
        self.assertEqual(row['physical_band'], dict(lower=F(1, 96), upper=24))
        self.assertEqual(row['unnormalized_band_mass_lower'], F(1, 4))
        self.assertEqual(row['normalized_band_mass_lower'], F(1, 2))
        self.assertEqual(row['normalized_physical_mean_energy'], dict(lower=F(1, 24), upper=6))
        self.assertTrue(row['survival_premise_is_not_discharged'])

    def test_independent_finite_survival_fixture_passes_both_tail_bounds(self):
        row = z.finite_survival_control((F(1, 100), 1, 100), (1, F(1, 2), F(1, 1000)), 1, 1, F(1, 8))
        self.assertGreater(row['control']['filtered_norm_squared']['lower'], F(1, 8))
        self.assertGreaterEqual(row['independently_certified_band_mass_lower'], F(1, 16))
        self.assertLessEqual(row['low_mass_upper'], F(1, 32))
        self.assertLessEqual(row['high_mass_upper'], F(1, 32))
        row_integer_floor = z.finite_survival_control((1,), (4,), 1, 4, 1)
        self.assertIsInstance(row_integer_floor['independently_certified_band_mass_lower'], F)

    def test_filter_controls_reject_false_energy_and_survival_claims(self):
        with self.assertRaises(ValueError): z.finite_filter_control((2,), (1,), Cf=1)
        with self.assertRaises(ValueError): z.finite_filter_control((), (), Cf=1)
        with self.assertRaises(ValueError): z.finite_filter_control((1,), (-1,), Cf=1)
        with self.assertRaises(ValueError): z.filter_survival_budget(1, 1, 2)
        with self.assertRaises(ValueError): z.filter_survival_budget(0, 1, F(1, 4))
        with self.assertRaises(ValueError): z.finite_survival_control((1000,), (F(1, 1000),), 1, 1, F(1, 4))

    def test_UV_energy_one_fixture_loses_all_filtered_norm(self):
        for n in (2, 4, 32, 1024):
            row = z.finite_filter_control((n,), (F(1, n),))
            self.assertEqual(row['source_energy'], 1)
            self.assertEqual(row['source_norm_squared'], F(1, n))
            self.assertLessEqual(row['filtered_norm_squared']['upper'], F(1, n))
            self.assertLessEqual(row['ordinary_increment_variance']['upper'], F(2, n))

    def test_IR_energy_one_fixture_has_increment_survival_but_filter_collapse(self):
        for n in (2, 4, 32, 1024):
            row = z.finite_filter_control((F(1, n),), (n,))
            self.assertEqual(row['source_energy'], 1)
            self.assertEqual(row['source_norm_squared'], n)
            self.assertLessEqual(row['filtered_norm_squared']['upper'], F(1, n))
            self.assertGreaterEqual(row['ordinary_increment_variance']['lower'], 2-F(1, n))
            self.assertLessEqual(row['ordinary_increment_variance']['upper'], 2)

    def test_spectral_fixtures_do_not_assert_actual_vacuum_or_nonzero_local_limit(self):
        control = z.all_controls()
        self.assertEqual(control['quaternion_trace_controls']['exact_triples_checked'], 216)
        for key in ('ultraviolet_filter_collapse', 'infrared_filter_collapse'):
            self.assertTrue(all(row['algebraic_fixture_not_Yang_Mills_vacuum'] for row in control[key]))
        self.assertFalse(z.MODEL['nonzero_anchored_limit'])
        self.assertFalse(z.MODEL['nontrivial_temporal_state'])


class PhysicalSurvivalCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_manifest_and_parent_result_bind_the_declared_sources(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 274)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertEqual(self.cert['parent_result_sha256'], z.PARENT_HASH)
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_upstream_and_parent_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_nonzero_continuum_neutral_gap_or_covariance_overclaims_rejected(self):
        for field in ('nonzero_anchored_limit', 'continuum_Yang_Mills', 'uniform_positive_gap',
                      'center_sector_bound_applies_to_center_neutral_gap', 'full_time_Euclidean_covariance'):
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

    def test_model_separates_constructed_temporal_framework_and_missing_YM_inputs(self):
        scope = self.cert['model']
        self.assertTrue(scope['actual_anchored_physical_space_time_compactness'])
        self.assertTrue(scope['anchored_reanchored_temporal_semigroup_constructed'])
        self.assertTrue(scope['actual_explicit_weak_coupling_gap_closing_trajectory'])
        self.assertFalse(scope['local_field_identification'])
        self.assertFalse(scope['heat_filter_survival_lower_bound_proved'])
        self.assertFalse(scope['ordinary_increment_survival_implies_filter_survival'])
        self.assertFalse(scope['finite_controls_are_universal_analytic_proof'])
        self.assertFalse(scope['winding_states_are_local_continuum_fields'])


if __name__ == '__main__':
    unittest.main()
