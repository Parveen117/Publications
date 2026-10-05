import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym78_stationary_derivative_energy_survival as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class DerivativeEnergyMathematicsTests(unittest.TestCase):
    def test_sparse_native_geometry_has_disjoint_radial_gradients_without_independence(self):
        for side in (4, 6):
            row = z.sparse_geometry(side)
            self.assertEqual(row['selected_plaquettes'], side**3//4)
            self.assertEqual(row['occupied_links'], side**3)
            self.assertEqual(row['rectangular_cell_sides_in_lattice_units'], [2, 2, 1])
            self.assertEqual(row['density_among_vertices'], F(1, 4))
            self.assertGreater(row['nonpole_native_radial_diagonals_checked'], 0)
            self.assertTrue(row['distinct_angle_gradients_have_disjoint_link_support'])
            self.assertTrue(row['radial_Dirichlet_Gram_is_identity_almost_everywhere'])
            self.assertFalse(row['interacting_vacuum_independence_assumed'])

    def test_rectangular_projection_preserves_constants_and_cell_volume(self):
        for n in (1, 2, 4):
            row = z.sparse_cell_projection({(0, 0, 0): -3}, n)
            self.assertEqual(row['spacing'], F(1, 2*n))
            self.assertEqual(row['cell_count'], 2*n**3)
            self.assertEqual(row['cell_volume']*row['cell_count'], 1)
            self.assertEqual(row['projected_L2_squared'], 9)
            self.assertEqual(row['exact_projection_error'], 0)
            self.assertEqual(row['conservative_Poincare_lower'], 9)
            self.assertEqual(row['continuum_gradient_L2_squared'], 0)

    def test_anisotropic_linear_projection_has_independent_exact_errors(self):
        for n in (1, 2, 4):
            x = z.sparse_cell_projection({(1, 0, 0): 1}, n)
            vertical = z.sparse_cell_projection({(0, 0, 1): 1}, n)
            mixed = z.sparse_cell_projection({(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}, n)
            self.assertEqual(x['exact_projection_error'], F(1, 12*n*n))
            self.assertEqual(vertical['exact_projection_error'], F(1, 48*n*n))
            self.assertEqual(mixed['exact_projection_error'], F(3, 16*n*n))
            self.assertEqual(mixed['continuum_L2_squared'], F(5, 2))
            self.assertEqual(mixed['continuum_gradient_L2_squared'], 3)

    def test_signed_quadratic_projection_obeys_Jensen_Poincare_and_refinement(self):
        poly = {(2, 0, 0): 1, (0, 1, 1): -1, (0, 0, 0): F(1, 3)}
        rows = [z.sparse_cell_projection(poly, n) for n in (1, 2, 4)]
        for row in rows:
            self.assertEqual(row['continuum_L2_squared'], F(14, 45))
            self.assertEqual(row['continuum_gradient_L2_squared'], 2)
            self.assertLessEqual(row['conservative_Poincare_lower'], row['projected_L2_squared'])
            self.assertLessEqual(row['projected_L2_squared'], F(14, 45))
        self.assertGreater(rows[0]['exact_projection_error'], rows[1]['exact_projection_error'])
        self.assertGreater(rows[1]['exact_projection_error'], rows[2]['exact_projection_error'])
        self.assertEqual(z.sparse_cell_projection({}, 2)['projected_L2_squared'], 0)

    def test_sparse_geometry_and_projection_reject_invalid_domains(self):
        for side in (2, 5, True, 4.0):
            with self.assertRaises(ValueError): z.sparse_geometry(side)
        for n in (0, -1, True):
            with self.assertRaises(ValueError): z.sparse_cell_projection({}, n)
        with self.assertRaises(ValueError): z.sparse_cell_projection({(1, 0): 1}, 2)
        with self.assertRaises(ValueError): z.sparse_cell_projection({(-1, 0, 0): 1}, 2)
        with self.assertRaises(ValueError): z.sparse_cell_projection({(1, 0, 0): 1.0}, 2)

    def test_exact_sparse_physical_energy_has_no_remaining_theta_factor(self):
        cell = z.sparse_cell_projection({(1, 0, 0): 1}, 3)
        for theta in (25, 256, 10**12):
            row = z.smear_moment_budget(theta, cell['spacing'], cell['continuum_L1_squared_upper'],
                cell['continuum_L2_squared'], cell['continuum_gradient_L2_squared'], 2,
                exact_projected_L2_squared=cell['projected_L2_squared'])
            self.assertEqual(row['physical_first_moment'],
                             dict(lower=2*cell['projected_L2_squared'], upper=2*cell['projected_L2_squared']))
            self.assertEqual(row['smear_coefficient_fourth_power'], 16*theta*cell['spacing']**8)
            self.assertTrue(row['actual_sparse_exact_energy_identity'])
            self.assertTrue(row['bounds_are_uniform_in_volume'])
            self.assertFalse(row['positive_continuum_reflected_norm_inferred'])

    def test_sparse_moment_constants_and_clock_scaling_are_exact(self):
        row = z.smear_moment_budget(25, F(1, 4), 1, 1, 0, 2)
        self.assertEqual(row['physical_first_moment'], dict(lower=2, upper=2))
        self.assertEqual(row['physical_source_norm_squared_upper'], 45)
        self.assertEqual(row['physical_second_moment_upper'], 10240)
        rough = z.smear_moment_budget(25, 1, 1, 1, 100)
        self.assertEqual(rough['physical_first_moment']['lower'], 0)
        self.assertEqual(rough['physical_first_moment']['upper'], 1)

    def test_moment_budgets_reject_false_projection_and_spectral_premises(self):
        with self.assertRaises(ValueError): z.smear_moment_budget(24, F(1, 4), 1, 1, 0)
        with self.assertRaises(ValueError): z.smear_moment_budget(25, 0, 1, 1, 0)
        with self.assertRaises(ValueError): z.smear_moment_budget(25, 1, 1, 1, 0, clock_scale=0)
        with self.assertRaises(ValueError): z.smear_moment_budget(25, 1, 0, 1, 0)
        with self.assertRaises(ValueError):
            z.smear_moment_budget(25, F(1, 4), 1, 1, 0, exact_projected_L2_squared=F(1, 2))
        with self.assertRaises(ValueError):
            z.smear_moment_budget(25, F(1, 4), 1, 1, 0, exact_projected_L2_squared=2)

    def test_reflected_moment_lower_bound_uses_distinct_first_moment_enclosures(self):
        row = z.reflected_norm_lower(F(3, 2), 2, F(25, 8), 4, F(5, 2))
        self.assertEqual(row['exact_denominator'], F(127, 8))
        self.assertEqual(row['certified_reflected_norm_lower'], F(72, 127))
        with localcontext() as context:
            context.prec = 90
            true = 2*(1-Decimal('-1.5').exp())**2+(1-Decimal(-6).exp())**2/8
            self.assertLess(decimal(row['certified_reflected_norm_lower']), true)
        self.assertTrue(row['fixed_regulator_lower_bound_only'])

    def test_positive_regulator_lower_bounds_can_still_collapse_as_spacing_to_fourth(self):
        for spacing in (F(1, 4), F(1, 8), F(1, 32)):
            moments = z.smear_moment_budget(256, spacing, 1, 1, 0)
            row = z.reflected_norm_lower(1, 1, moments['physical_source_norm_squared_upper'],
                                        moments['physical_second_moment_upper'])
            exact = spacing**4/(10+F(45, 16)*spacing**2+spacing**4)
            self.assertEqual(row['certified_reflected_norm_lower'], exact)
            self.assertGreater(exact, 0)
            self.assertLess(exact, spacing**4/10)

    def test_reflected_moment_bound_handles_zero_energy_and_rejects_inconsistency(self):
        self.assertEqual(z.reflected_norm_lower(1, 0, 0, 0)['certified_reflected_norm_lower'], 0)
        self.assertEqual(z.reflected_norm_lower(1, 0, 10, 0)['certified_reflected_norm_lower'], 0)
        with self.assertRaises(ValueError): z.reflected_norm_lower(0, 1, 1, 1)
        with self.assertRaises(ValueError): z.reflected_norm_lower(1, 2, 1, 1)
        with self.assertRaises(ValueError): z.reflected_norm_lower(1, 1, 1, 1, F(1, 2))
        with self.assertRaises(ValueError): z.reflected_norm_lower(1, -1, 1, 1)

    def test_scalar_filter_inequalities_match_independent_exponentials_and_remainders(self):
        with localcontext() as context:
            context.prec = 90
            for x in (0, F(1, 1000), F(1, 2), 1, 2, 3, 100):
                row = z.scalar_filter_inequalities(x)
                true = (1-(-decimal(x)).exp())**2
                self.assertLessEqual(decimal(row['reflected_norm_scalar_lower']), true)
                self.assertLessEqual(decimal(row['contrast_reflected_norm_scalar_lower']), true)
                self.assertEqual(row['quadratic_comparison_remainder'], F(x*x*(1+x), 4))
                self.assertEqual(row['contrast_quartic'], row['nonnegative_quartic_decomposition'])
                self.assertGreaterEqual(row['contrast_quartic'], 0)
                self.assertLessEqual(row['frequency_contrast_kernel'], F(2, 3))

    def test_stationary_derivative_moments_have_odd_reflection_and_required_pointwise_control(self):
        row = z.derivative_moment_budget(2, F(1, 3), 3)
        self.assertEqual(row['mgf_quadratic_coefficient'], F(2, 3))
        self.assertEqual(row['variance_upper'], F(4, 3))
        self.assertEqual(row['even_moment_upper'], F(2048, 9))
        smaller = z.derivative_moment_budget(2, F(1, 12), 3)
        self.assertEqual(smaller['even_moment_upper'], row['even_moment_upper']/64)
        self.assertTrue(row['derivative_reflection_is_odd'])
        self.assertTrue(row['uniform_pointwise_Dirichlet_bound_required'])
        self.assertEqual(z.derivative_moment_budget(0, 3, 2)['even_moment_upper'], 0)

    def test_derivative_kernel_ignores_zero_atoms_and_has_independent_rational_values(self):
        atoms, weights = (0, 1, 3), (100, F(1, 2), F(1, 3))
        self.assertEqual(z.derivative_kernel(atoms, weights, 2), F(4, 5)+F(8, 13))
        self.assertEqual(z.derivative_kernel(atoms, weights, 2),
                         z.derivative_kernel(atoms[1:], weights[1:], 2))
        self.assertLess(z.derivative_kernel(atoms, weights, 1), z.derivative_kernel(atoms, weights, 2))
        self.assertLess(z.derivative_kernel(atoms, weights, 2), 3)

    def test_reflected_spectral_interval_encloses_an_independent_decimal_sum(self):
        atoms, weights, frequency = (0, F(1, 100), 1, 100), (5, 2, F(1, 3), F(1, 1000)), F(3, 2)
        row = z.derivative_spectral_control(atoms, weights, frequency)
        with localcontext() as context:
            context.prec = 90
            true = sum(decimal(w)*(1-(-decimal(E/frequency)).exp())**2 for E, w in zip(atoms, weights))
            self.assertLessEqual(decimal(row['reflected_norm_interval']['lower']), true)
            self.assertGreaterEqual(decimal(row['reflected_norm_interval']['upper']), true)
        self.assertEqual(row['frequency_contrast'],
                         row['derivative_kernel_at_double_frequency']-row['derivative_kernel'])
        self.assertLessEqual(row['frequency_contrast'], F(2, 3)*row['source_energy'])
        self.assertTrue(row['algebraic_fixture_not_interacting_vacuum'])

    def test_fixed_energy_has_positive_frequency_contrast_and_reflected_survival(self):
        row = z.derivative_spectral_control((1,), (1,))
        self.assertEqual(row['source_energy'], 1)
        self.assertEqual(row['derivative_kernel'], 1)
        self.assertEqual(row['derivative_kernel_at_double_frequency'], F(8, 5))
        self.assertEqual(row['frequency_contrast'], F(3, 5))
        self.assertEqual(row['contrast_reflected_norm_lower'], F(1, 5))
        self.assertEqual(row['moment_reflected_norm_lower'], F(1, 3))
        self.assertGreaterEqual(row['reflected_norm_interval']['lower'], F(1, 3))

    def test_UV_energy_one_escapes_both_fixed_frequency_kernel_and_reflected_norm(self):
        for n in (2, 4, 32, 256):
            row = z.derivative_spectral_control((n,), (F(1, n),))
            self.assertEqual(row['source_energy'], 1)
            self.assertEqual(row['derivative_kernel'], F(2, n*n+1))
            self.assertEqual(row['frequency_contrast'], F(6*n*n, (n*n+1)*(n*n+4)))
            self.assertLessEqual(row['frequency_contrast'], F(6, n*n))
            self.assertLessEqual(row['reflected_norm_interval']['upper'], F(1, n))

    def test_IR_energy_one_keeps_white_noise_kernel_but_loses_contrast_and_reflected_norm(self):
        for n in (2, 4, 32, 256):
            row = z.derivative_spectral_control((F(1, n),), (n,))
            self.assertEqual(row['source_energy'], 1)
            self.assertEqual(row['derivative_kernel'], F(2*n*n, 1+n*n))
            self.assertGreaterEqual(row['derivative_kernel'], 2-F(2, n*n))
            self.assertEqual(row['frequency_contrast'], F(6*n*n, (1+n*n)*(1+4*n*n)))
            self.assertLessEqual(row['frequency_contrast'], F(3, 2*n*n))
            self.assertLessEqual(row['reflected_norm_interval']['upper'], F(1, n))

    def test_dyadic_telescoping_recovers_energy_and_separates_missing_frequency_tails(self):
        rows = [z.dyadic_contrast_control((0, F(1, 8), 1, 32), (100, 2, F(1, 3), F(1, 64)),
                                         base_frequency=F(3, 2), levels=n) for n in (1, 3, 6)]
        for levels, row in zip((1, 3, 6), rows):
            self.assertEqual(len(row['terms']), 2*levels)
            self.assertEqual(row['exact_sum'], row['endpoint_difference'])
            self.assertEqual(row['all_frequency_sum'], F(13, 6))
            self.assertEqual(row['exact_sum']+row['missing_frequency_tails'], F(13, 6))
            self.assertTrue(all(term['contrast'] >= 0 for term in row['terms']))
            self.assertFalse(row['fixed_frequency_window_survival_inferred'])
        self.assertGreater(rows[0]['missing_frequency_tails'], rows[1]['missing_frequency_tails'])
        self.assertGreater(rows[1]['missing_frequency_tails'], rows[2]['missing_frequency_tails'])
        for E in (F(2)**12, F(2)**-12):
            escaped = z.dyadic_contrast_control((E,), (1/E,), levels=3)
            self.assertEqual(escaped['all_frequency_sum'], 2)
            self.assertLess(escaped['exact_sum'], F(1, 100))

    def test_spectral_and_derivative_controls_reject_invalid_domains(self):
        for atoms, weights in (((), ()), ((1, 2), (1,)), ((-1,), (1,)), ((1,), (-1,)), ((1.0,), (1,))):
            with self.assertRaises(ValueError): z.derivative_spectral_control(atoms, weights)
        with self.assertRaises(ValueError): z.derivative_kernel((1,), (1,), 0)
        with self.assertRaises(ValueError): z.dyadic_contrast_control((1,), (1,), levels=0)
        with self.assertRaises(ValueError): z.derivative_moment_budget(1, 1, True)
        with self.assertRaises(ValueError): z.derivative_moment_budget(-1, 1, 2)
        with self.assertRaises(ValueError): z.scalar_filter_inequalities(-1)

    def test_neutral_cube_removes_only_charged_modes_and_keeps_the_previous_gap_interval(self):
        row = z.neutral_cube_controls()
        self.assertEqual(row['full_R3_dimension'], 241)
        self.assertEqual(row['neutral_R3_dimension'], 193)
        self.assertEqual(row['removed_charged_winding_dimension'], 48)
        charges = {tuple(item['charge']): item['count'] for item in row['winding_character_multiplicities']}
        self.assertEqual(charges, {(1, 0, 0): 16, (0, 1, 0): 16, (0, 0, 1): 16})
        self.assertEqual(row['actual_neutral_hidden_free_floor'], F(9, 2))
        self.assertEqual(row['interacting_hidden_floor'], F(1, 2))
        self.assertEqual(row['exact_neutral_thresholds'][0]['neutral_inertia'], dict(negative=0, zero=0, positive=193))
        self.assertEqual(row['exact_neutral_thresholds'][1]['neutral_inertia'], dict(negative=1, zero=0, positive=192))
        self.assertEqual(row['certified_neutral_gap'], z.stage75.actual_cube_controls()['certified_bare_gap'])
        self.assertTrue(row['previous_full_gap_same_interval'])
        self.assertTrue(row['inverse_must_be_taken_on_neutral_hidden_complement'])
        self.assertFalse(row['neutral_sector_identified_with_local_continuum_core'])


class DerivativeEnergyCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_manifest_and_parent_result_bind_declared_sources(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 281)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertEqual(self.cert['parent_result_sha256'], z.PARENT_HASH)
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_upstream_and_parent_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_survival_gap_covariance_and_neutral_locality_overclaims_are_rejected(self):
        for field in ('fixed_frequency_survival_window_certified', 'nonzero_local_reflected_continuum_state',
                      'positive_continuum_mass_gap', 'full_four_dimensional_Euclidean_covariance',
                      'center_neutral_sector_equals_local_continuum_core',
                      'charged_resolvent_poles_may_be_ignored_in_full_inverse'):
            wrong = copy.deepcopy(self.cert)
            wrong['model'][field] = True
            with tempfile.TemporaryDirectory() as directory:
                result, pin = Path(directory)/'result', Path(directory)/'pin'
                result.write_text(json.dumps(wrong))
                pin.write_text(z.p.canonical_sha(wrong))
                with self.assertRaises(ValueError): z.check(wrong, result, pin)
                with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_certificate_contains_no_binary_floats_or_unintended_control_characters(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_model_separates_stationary_derivative_construction_from_missing_continuum_inputs(self):
        scope = self.cert['model']
        for field in ('actual_stationary_tempered_space_time_derivative_limit',
                      'derivative_reflection_uses_odd_time_parity',
                      'actual_sparse_local_smear_has_exact_physical_energy',
                      'actual_finite_regulator_reflected_norm_lower_bound',
                      'local_sources_have_exactly_zero_charged_center_Schur_coupling',
                      'actual_neutral_four_cube_cut_dimension_193'):
            self.assertTrue(scope[field])
        for field in ('fixed_frequency_survival_window_certified', 'nonzero_local_reflected_continuum_state',
                      'local_Yang_Mills_field_identification', 'positive_continuum_mass_gap',
                      'neutral_restriction_improves_previous_numeric_gap_interval',
                      'finite_replays_are_universal_proof', 'spectral_fixtures_are_interacting_Yang_Mills_vacua'):
            self.assertFalse(scope[field])


if __name__ == '__main__':
    unittest.main()
