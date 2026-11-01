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
import ym76_local_defect_weak_coupling_dynamics as z


def decimal(value):
    return Decimal(value.numerator)/Decimal(value.denominator)


class LocalDynamicsMathematicsTests(unittest.TestCase):
    def test_local_gradients_match_independent_native_polynomial_derivatives(self):
        qs = ((F(3, 5), F(4, 5), 0, 0), (F(5, 13), 0, F(12, 13), 0),
              (F(8, 17), 0, 0, F(15, 17)), (0, 1, 0, 0))
        for signs in ((1, 1, -1, -1), (-1, 1, -1, 1), (1, -1, 1, -1)):
            word = tuple(enumerate(signs))
            poly = z.model.loop_polynomial(word)
            value, gradients = z.face_gradient(qs, word)
            self.assertEqual(value, z.model.evaluate(poly, qs))
            for e in range(4):
                for axis in (1, 2, 3):
                    self.assertEqual(gradients[e][axis-1],
                                     z.model.evaluate(z.native.derivative(poly, e, axis), qs))

    def test_common_link_cross_gradient_matches_independent_polynomial(self):
        qs = ((F(3, 5), F(4, 5), 0, 0), (F(5, 13), 0, F(12, 13), 0),
              (F(8, 17), 0, 0, F(15, 17)), (0, 1, 0, 0))
        words = (((0, 1), (1, 1), (2, -1), (3, -1)),
                 ((2, 1), (0, -1), (3, 1), (1, -1)))
        polys = [z.model.loop_polynomial(word) for word in words]
        gradients = [z.face_gradient(qs, word)[1] for word in words]
        native_cross = z.p.add(*(z.p.mul(z.native.derivative(polys[0], 0, axis),
                                        z.native.derivative(polys[1], 0, axis))
                                 for axis in (1, 2, 3)))
        dot = sum(x*y for x, y in zip(gradients[0][0], gradients[1][0]))
        self.assertNotEqual(dot, 0)
        self.assertEqual(dot, z.model.evaluate(native_cross, qs))

    def test_inverse_links_and_closed_loop_gauge_invariance(self):
        graph, qs, values, _ = z.geometry_data(4)
        palette = (z.model.y.ONE, (F(3, 5), F(4, 5), 0, 0), (0, 0, 1, 0))
        gauges = {vertex: palette[sum(vertex) % 3] for vertex in graph['vertices']}
        transformed = z.model.gauge_transform(qs, gauges, graph['ends'])
        for face, expected in zip(graph['faces'][:12], values[:12]):
            value, gradient = z.face_gradient(transformed, face)
            self.assertEqual(value, expected)
            self.assertEqual(sum(map(z.norm2, gradient.values())), 1-value*value)

    def test_exact_input_and_geometric_domains(self):
        for bad in ((1, 1, 0, 0), (1, 0, 0), (True, 0, 0, 0), (1.0, 0, 0, 0)):
            with self.assertRaises(ValueError): z.quaternion(bad)
        with self.assertRaises(ValueError): z.geometry_data(3)
        with self.assertRaises(ValueError): z.geometry_data(4, 'unknown')
        with self.assertRaises(ValueError): z.face_gradient((z.model.y.ONE,)*4, ((0, 1),))
        with self.assertRaises(ValueError):
            z.face_gradient((z.model.y.ONE,)*4, ((0, 1), (1, 1), (2, -1), (2, -1)))

    def test_actual_periodic_side_four_and_five_have_the_required_incidence(self):
        for side in (4, 5):
            row = z.geometry_replay(side)
            self.assertEqual(row['plaquettes'], 3*side**3)
            self.assertEqual(row['links'], row['plaquettes'])
            self.assertEqual(row['plaquettes_per_link'], 4)
            self.assertEqual(row['compared_incident_pairs'], 18*side**3)

    def test_local_positive_defect_decomposes_exactly_on_noncommuting_configurations(self):
        for side in (4, 5):
            row = z.geometry_replay(side)
            self.assertEqual(8*row['X']-row['Gamma_V'], row['exact_nonnegative_defect'])
            self.assertEqual(row['exact_nonnegative_defect'],
                             row['local_square_deficiency']+row['local_gradient_mismatch'])
            self.assertGreater(row['local_square_deficiency'], 0)
            self.assertGreater(row['local_gradient_mismatch'], 0)

    def test_flat_and_center_configs_separate_onsite_from_gradient_defects(self):
        for side in (4, 5):
            flat, center = z.geometry_replay(side, 'identity'), z.geometry_replay(side, 'center')
            self.assertEqual(flat['X'], 0)
            self.assertEqual(flat['exact_nonnegative_defect'], 0)
            self.assertEqual(center['Gamma_V'], 0)
            self.assertEqual(center['local_gradient_mismatch'], 0)
            self.assertGreater(center['X'], 0)
            self.assertEqual(center['local_square_deficiency'], 8*center['X'])

    def test_square_completion_matches_the_original_potential(self):
        row = z.geometry_replay(5)
        for control in row['square_completion_replays']:
            t, theta = control['t'], control['theta']
            self.assertEqual(theta, 3*t+8*t*t)
            V = row['plaquettes']-row['X']
            right = -3*t*V+t*t*row['Gamma_V']+3*t*row['plaquettes']+t*t*row['exact_nonnegative_defect']
            self.assertEqual(theta*row['X'], right)

    def test_collective_incidence_is_exact_without_independence(self):
        for side in (4, 5):
            row = z.collective_replay(side)
            self.assertFalse(row['independence_of_vacuum_plaquettes_assumed'])
            self.assertGreater(row['edge_disjoint_plaquette_count'], 1)
            for control in row['exact_congruence_controls']:
                self.assertEqual(4*control['radial_coefficient_norm_squared']-control['collective_Dirichlet_energy'],
                                 control['four_I_minus_K_quadratic_form'])
                self.assertGreaterEqual(control['four_I_minus_K_quadratic_form'], 0)

    def test_weak_budgets_use_actual_YM74_energy_bound_and_positive_square_floor(self):
        for theta in (25, 26, F(101, 4), 256, 10**6):
            row = z.weak_budget(theta)
            old = z.stage74.vacuum_budget(theta)
            self.assertEqual(row['vacuum_excess_density']['upper'], old['actual_shifted_vacuum_deficit_upper'])
            self.assertGreater(row['vacuum_excess_density']['lower'], 0)
            self.assertLessEqual(row['vacuum_excess_density']['lower'], row['vacuum_excess_density']['upper'])
            self.assertGreaterEqual(row['actual_trace_mean_lower'], F(1, 2))

    def test_directed_square_roots_enclose_irrational_weak_scales(self):
        with localcontext() as context:
            context.prec = 90
            for theta in (26, 257, 65537):
                row, true_root = z.weak_budget(theta), Decimal(theta).sqrt()
                self.assertLessEqual(decimal(row['sqrt_theta']['lower']), true_root)
                self.assertGreaterEqual(decimal(row['sqrt_theta']['upper']), true_root)
                tau, delta = row['actual_link_electric_energy'], row['actual_plaquette_deficiency']
                self.assertLessEqual(decimal(tau['lower']), true_root/32)
                self.assertGreaterEqual(decimal(tau['upper']), 9*true_root/4)
                self.assertLessEqual(decimal(delta['lower']), 1/(32*true_root))
                self.assertGreaterEqual(decimal(delta['upper']), 9/(4*true_root))

    def test_Ward_constants_follow_from_the_stated_margins(self):
        for s in (5, 16, 256):
            m_floor, q_upper, tau_upper = F(1, 2), F(9, 2*s), F(9*s, 4)
            self.assertEqual(9*m_floor*m_floor/(16*q_upper), F(s, 32))
            self.assertEqual(9*m_floor*m_floor/(32*tau_upper), F(1, 32*s))
            self.assertLessEqual(16*tau_upper+4, 40*s)

    def test_angle_variance_bounds_are_positive_and_contain_sharper_bounds(self):
        for theta in (25, 26, 256, 65536):
            row = z.weak_budget(theta)
            bound = row['actual_angle_variance']
            self.assertGreater(bound['lower'], 0)
            self.assertLessEqual(bound['lower'], row['sharper_angle_variance_lower'])
            self.assertLessEqual(row['sharper_angle_variance_upper'], bound['upper'])
            self.assertLess(row['pi_interval']['upper']**2, 10)
            self.assertEqual(row['explicit_rescaled_state_first_moment'], 1)
            self.assertLess(row['explicit_rescaled_state_norm_squared']['upper'], F(45, 4))

    def test_weak_simple_constants_have_an_explicit_domain(self):
        for bad in (-1, 0, 1, 24, F(2499, 100), True, 25.0):
            with self.assertRaises(ValueError): z.weak_budget(bad)

    def test_beta_radial_energy_has_four_link_not_one_link_normalization(self):
        for n in (0, 1, 2, 7, 64):
            row = z.radial_beta(n)
            self.assertEqual(row['radial_electric_energy'], F(3*n*n, 4*n+1))
            self.assertEqual(row['radial_electric_energy'], 4*z.stage74.product_trial(n)['link_electric_energy'])
            self.assertEqual(row['cotangent_squared_mean'], F(4*n*n+2*n+1, 4*n+1))
            self.assertEqual(row['cosine_over_one_plus_cosine_mean'], F(2*n-1, 4*n+1))

    def test_singular_radial_terms_cancel_even_when_each_is_large(self):
        for n in (0, 1, 1000):
            row = z.radial_beta(n)
            self.assertEqual(row['J_norm_squared']-row['radial_electric_energy'], 1)
            self.assertEqual(row['cotangent_squared_mean']-2*n*row['cosine_over_one_plus_cosine_mean'], 1)
        self.assertEqual(z.radial_beta(0)['J_norm_squared'], 1)
        for bad in (-1, True, F(1, 2)):
            with self.assertRaises(ValueError): z.radial_beta(bad)

    def test_free_Haar_angle_fixture_respects_uncertainty(self):
        row = z.all_controls()['free_Haar_angle']
        self.assertGreater(row['variance']['lower'], F(1, 4))
        self.assertLess(row['variance']['upper'], F(1, 3))
        self.assertEqual(row['Dirichlet_energy'], 1)
        self.assertEqual(row['second_energy_moment'], 4)

    def test_band_constants_and_independent_finite_measure(self):
        row = z.finite_spectral_control((0, F(1, 100), 1, 100),
                                       (F(1, 4), F(1, 4), F(499, 1000), F(1, 1000)))
        self.assertEqual(row['band_lower'], F(1, 45))
        self.assertEqual(row['band_upper'], 80)
        self.assertEqual(row['band_mass_lower'], F(1, 3600))
        self.assertEqual(row['actual_band_mass'], F(499, 1000))
        self.assertLess(row['second_moment'], 40*row['first_moment'])

    def test_generic_band_theorem_is_valid_with_boundary_atoms(self):
        for a, b in ((1, 1), (F(1, 4), 2), (2, 10)):
            for atoms, weights in (((a,), (1,)), ((0, b), (1-F(a, b), F(a, b)))):
                row = z.finite_spectral_control(atoms, weights, a, b)
                self.assertGreaterEqual(row['actual_band_mass'], F(a, 8*b))
        row = z.finite_spectral_control((F(1, 4), 2), (F(6, 7), F(1, 7)), F(1, 4), 2)
        self.assertEqual(row['actual_band_mass'], 1)

    def test_invalid_spectral_measures_or_missing_moment_bounds_rejected(self):
        for atoms, weights in (((), ()), ((1,), (F(1, 2),)), ((-1,), (1,)), ((1, 2), (1,))):
            with self.assertRaises(ValueError): z.finite_spectral_control(atoms, weights)
        with self.assertRaises(ValueError): z.finite_spectral_control((0,), (1,))
        with self.assertRaises(ValueError): z.finite_spectral_control((0, 1000), (F(999, 1000), F(1, 1000)))
        with self.assertRaises(ValueError): z.spectral_band_budget(2, 1)

    def test_positive_source_band_does_not_prevent_closing_full_gap_or_zero_atom(self):
        rows = z.all_controls()['closing_gap_source_fixtures']
        self.assertTrue(all(row['actual_band_mass'] >= F(1, 2) for row in rows))
        self.assertLess(rows[-1]['atoms'][0], rows[0]['atoms'][0])
        zero = z.all_controls()['limiting_zero_atom_fixture']
        self.assertEqual(zero['atoms'], (0, 1))
        self.assertEqual(zero['weights'][0], F(1, 2))
        self.assertFalse(zero['zero_atom_excluded'])
        self.assertFalse(zero['full_gap_inferred'])

    def test_first_moment_alone_allows_loss_of_every_fixed_finite_band(self):
        rows = z.all_controls()['first_moment_without_second_moment_failure']
        for row in rows:
            self.assertEqual(sum(x*w for x, w in zip(row['atoms'], row['weights'])), 1)
            self.assertEqual(sum(x*x*w for x, w in zip(row['atoms'], row['weights'])), row['parameter'])
            self.assertTrue(all(x == 0 or x > 80 for x in row['atoms']))

    def test_directed_exponentials_enclose_independent_decimal_values(self):
        with localcontext() as context:
            context.prec = 100
            for value in (F(0), F(1, 45), F(1, 2), F(1), F(3), F(40), F(1000)):
                bounds = z.exp_negative_interval(value)
                actual = (-decimal(value)).exp()
                self.assertLessEqual(decimal(bounds['lower']), actual)
                self.assertGreaterEqual(decimal(bounds['upper']), actual)
                self.assertLess(bounds['upper']-bounds['lower'], F(1, 10**32))
        for bad in (-1, 0.1, True):
            with self.assertRaises(ValueError): z.exp_negative_interval(bad)

    def test_heat_budgets_contain_independent_zero_atom_correlations(self):
        with localcontext() as context:
            context.prec = 90
            for time in (F(0), F(1, 100), F(1), F(10)):
                row = z.heat_budget(time)
                actual = (1+(-decimal(time)).exp())/2
                self.assertLessEqual(decimal(row['heat_correlation_lower']), actual)
                self.assertGreaterEqual(decimal(row['heat_correlation_upper']), actual)
        self.assertEqual(z.heat_budget(0)['heat_correlation_lower'], 1)
        self.assertEqual(z.heat_budget(0)['heat_correlation_upper'], 1)
        self.assertLess(z.heat_budget(1)['heat_correlation_upper'], 1)
        with self.assertRaises(ValueError): z.heat_budget(-1)

    def test_physical_band_escapes_to_ultraviolet_and_keeps_a_fixed_mass(self):
        coarse, fine = z.physical_clock_budget(1, 1), z.physical_clock_budget(1, F(1, 100))
        self.assertEqual(fine['physical_energy_band']['lower'], 100*coarse['physical_energy_band']['lower'])
        self.assertEqual(fine['physical_energy_band']['upper'], 100*coarse['physical_energy_band']['upper'])
        self.assertEqual(fine['escaping_spectral_mass_lower'], F(1, 3600))
        self.assertLess(fine['physical_heat_correlation_upper'], coarse['physical_heat_correlation_upper'])
        self.assertEqual(fine['positive_time_lim_sup_upper'], F(3599, 3600))
        self.assertFalse(fine['excludes_smeared_fields_or_remaining_infrared_component'])
        self.assertEqual(z.physical_clock_budget(0, F(1, 100))['positive_time_lim_sup_upper'], 1)
        with self.assertRaises(ValueError): z.physical_clock_budget(1, 0)
        with self.assertRaises(ValueError): z.physical_clock_budget(1, 1, 0)


class LocalDynamicsCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_source_manifest_and_parent_result_are_pinned(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 267)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertEqual(self.cert['parent_result_sha256'], z.PARENT_HASH)
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_upstream_and_parent_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_gap_continuum_or_probe_limit_overclaims_rejected(self):
        for field in ('uniform_weak_coupling_positive_gap', 'continuum_R4_QFT',
                      'unfiltered_probe_has_norm_preserving_strongly_continuous_physical_limit'):
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

    def test_scope_separates_written_theorem_replays_and_unproved_continuum(self):
        scope = self.cert['model']
        self.assertTrue(scope['actual_SU2_local_defect_factorization'])
        self.assertTrue(scope['actual_microscopic_spectral_band_weight'])
        self.assertTrue(scope['physical_ultraviolet_weight_escape'])
        self.assertTrue(scope['zero_atom_allowed_in_limiting_source_spectrum'])
        self.assertFalse(scope['finite_geometry_replay_is_universal_proof'])
        self.assertFalse(scope['radial_beta_fixtures_are_interacting_vacua'])
        self.assertFalse(scope['microscopic_clock_is_physical_clock'])
        self.assertFalse(scope['physical_mass_prediction'])


if __name__ == '__main__':
    unittest.main()
