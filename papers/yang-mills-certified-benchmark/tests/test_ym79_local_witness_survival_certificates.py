import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym79_local_witness_survival_certificates as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class LocalWitnessMathematicsTests(unittest.TestCase):
    def test_actual_vacuum_distance_uses_phase_and_the_sharp_overlap_conversion(self):
        row = z.vacuum_distance_budget(-F(3, 433), -F(712209, 10**8), F(48376192, 10**8))
        self.assertEqual(row['excited_weight_upper'], F(8386497, 433*49088401))
        self.assertLess(row['excited_weight_upper'], F(9999, 25000000))
        self.assertGreater(2*row['excited_weight_upper'], F(1, 2500))
        self.assertLessEqual(row['distance_upper_after_certified_phase_alignment'], F(1, 50))
        self.assertTrue(row['phase_alignment_must_be_verified'])
        with localcontext() as context:
            context.prec = 90
            true = (2*(1-(1-decimal(row['excited_weight_upper'])).sqrt())).sqrt()
            self.assertGreaterEqual(decimal(row['distance_upper_after_certified_phase_alignment']), true)

    def test_vacuum_line_distance_handles_exact_ground_states_and_rejects_bad_spectral_order(self):
        exact = z.vacuum_distance_budget(0, 0, 1)
        self.assertEqual(exact['excited_weight_upper'], 0)
        self.assertEqual(exact['distance_upper_after_certified_phase_alignment'], 0)
        for args in ((-1, 0, 1), (1, 0, 1), (2, 0, 1), (0, 0, 0)):
            with self.assertRaises(ValueError): z.vacuum_distance_budget(*args)

    def test_native_single_holonomy_moments_match_Catalan_integrals(self):
        row = z.native_source_moments()
        for n, key in ((1, 'W2'), (2, 'W4'), (3, 'W6')):
            self.assertEqual(row[key], F(comb(2*n, n), (n+1)*4**n))
        self.assertEqual(row['V_W2'], 0)
        self.assertEqual(row['V_W4'], 0)
        self.assertEqual(row['V2'], 48)
        self.assertEqual(row['V3'], 0)
        self.assertEqual(row['V_LE_V'], 144)

    def test_actual_cube_parity_and_private_links_give_mixed_moments_without_independence(self):
        row = z.native_source_moments()
        self.assertEqual(row['distinct_face_boundaries'], 192)
        self.assertEqual(row['odd_link_pair_checks'], comb(192, 2))
        self.assertTrue(row['no_three_plaquette_parity_cancellation'])
        self.assertEqual(row['other_faces_with_a_private_link'], 191)
        self.assertEqual(row['minimum_private_links'], 3)
        self.assertEqual(row['V2_W2'], F(1, 8)+191*F(1, 16))
        self.assertEqual(row['V2_W4'], F(5, 64)+191*F(1, 32))
        self.assertFalse(row['whole_plaquette_independence_assumed'])

    def test_actual_untruncated_local_witness_has_the_declared_rational_lower_bound(self):
        row = z.actual_cube_witness()
        self.assertEqual(row['theta'], F(1, 48))
        self.assertEqual(row['side'], 4)
        self.assertEqual(row['normalization_squared'], F(433, 432))
        self.assertEqual(row['vacuum_approximation']['rayleigh'], -F(3, 433))
        self.assertTrue(row['phase_verified_by_actual_eigen_equation'])
        self.assertEqual(row['rational_hand_witness']['signed_variational_value'], F(32607, 250000))
        self.assertGreater(row['rational_hand_witness']['signed_variational_value'], F(13, 100))
        self.assertGreater(row['sharper_directed_witness']['signed_variational_value'], F(32607, 250000))
        self.assertTrue(row['hidden_space_untruncated'])
        self.assertFalse(row['trial_identified_with_actual_vacuum'])
        self.assertFalse(row['uniform_volume_or_continuum_survival_asserted'])

    def test_dual_witness_matches_independent_scalar_resolvent_and_heat_inequalities(self):
        with localcontext() as context:
            context.prec = 90
            for energy, time in ((F(1, 10), F(2)), (F(1), F(1)), (F(3), F(1, 2)), (F(100), F(1))):
                coefficient = time*energy/(1+time*energy)**2
                row = z.dual_witness_bound(time, coefficient, energy, (1+time*energy)**2)
                self.assertEqual(row['signed_variational_value'], (time*energy/(1+time*energy))**2)
                heat = (1-(-decimal(time*energy)).exp())**2
                self.assertLessEqual(decimal(row['certified_heat_filtered_norm_squared_lower']), heat)
        self.assertEqual(z.dual_witness_bound(1, 0, 100, 1)['signed_variational_value'], 0)

    def test_dual_witness_does_not_turn_an_unfavorable_pairing_into_survival(self):
        row = z.dual_witness_bound(1, 1, -1, 1)
        self.assertEqual(row['signed_variational_value'], -3)
        self.assertEqual(row['certified_heat_filtered_norm_squared_lower'], 0)
        with self.assertRaises(ValueError): z.dual_witness_bound(0, 1, 1, 1)
        with self.assertRaises(ValueError): z.dual_witness_bound(1, -1, 1, 1)
        with self.assertRaises(ValueError): z.dual_witness_bound(1, 1, 1, -1)

    def test_source_residual_budget_keeps_each_error_and_ground_shift_price(self):
        row = z.source_residual_budget(2, 3, F(1, 10), F(1, 4), 2, 2, 3, F(1, 50), F(1, 25))
        self.assertEqual(row['commuting_ground_shift_operator_error'], F(1, 3))
        self.assertEqual(row['vacuum_source_error'], F(3, 10))
        self.assertEqual(row['ground_shift_source_error'], F(2, 3))
        self.assertEqual(row['total_response_error'], F(77, 75))
        self.assertEqual(row['resolvent_response_norm'], dict(lower=F(73, 75), upper=F(302, 75)))
        self.assertEqual(row['resolvent_norm_squared']['lower'], F(73, 75)**2)
        self.assertTrue(row['lower_A_bound_and_phase_alignment_required'])
        self.assertTrue(row['inverse_uses_full_neutral_operator_not_a_hidden_truncation'])

    def test_commuting_ground_shift_bound_is_sharp_at_zero_energy(self):
        time, width = F(3, 2), F(2, 5)
        bound = z.source_residual_budget(time, 1, 0, width, 1, 0, 1, 0)
        price = bound['commuting_ground_shift_operator_error']
        self.assertEqual(price, F(3, 8))
        for energy in (0, F(1, 10), 1, 100):
            exact_difference = time*(energy+width)/(1+time*(energy+width))-time*energy/(1+time*energy)
            self.assertLessEqual(exact_difference, price)
            if not energy: self.assertEqual(exact_difference, price)

    def test_matrix_residual_certificates_enclose_the_exact_full_response(self):
        rows = z.residual_matrix_controls()
        self.assertEqual(len(rows), 3)
        widths = []
        for row in rows:
            interval = row['budget']['resolvent_norm_squared']
            self.assertEqual(row['exact_source_response_norm_squared'], F(4, 9))
            self.assertLessEqual(interval['lower'], F(4, 9))
            self.assertGreaterEqual(interval['upper'], F(4, 9))
            self.assertTrue(row['matrix_fixture_not_Yang_Mills'])
            widths.append(interval['upper']-interval['lower'])
        self.assertGreater(widths[0], widths[1])
        self.assertGreater(widths[1], widths[2])
        self.assertLess(widths[2], F(1, 10**30))
        self.assertEqual(rows[-1]['budget']['total_response_error'], 0)

    def test_source_residual_budget_rejects_invalid_norms_and_clamps_uncertain_lower_bounds(self):
        row = z.source_residual_budget(1, 1, 2, 0, 1, 1, 1, 0)
        self.assertEqual(row['resolvent_response_norm']['lower'], 0)
        with self.assertRaises(ValueError): z.source_residual_budget(1, 1, 0, 0, 1, 2, 1, 0)
        with self.assertRaises(ValueError): z.source_residual_budget(1, 1, -1, 0, 1, 1, 1, 0)
        with self.assertRaises(ValueError): z.source_residual_budget(1, 1, 0, 0, 1, 1, 1, -1)
        with self.assertRaises(ValueError): z.source_residual_budget(0, 1, 0, 0, 1, 1, 1, 0)

    def test_bounded_Weyl_increment_bridge_uses_centered_mass_and_retains_the_premise(self):
        row = z.weyl_increment_budget(F(1, 2), F(3, 4))
        self.assertEqual(row['reflected_norm_lower'], F(1, 12))
        self.assertEqual(row['weaker_universal_reflected_lower'], F(1, 16))
        self.assertTrue(row['actual_positive_increment_premise_still_required'])
        self.assertEqual(z.weyl_increment_budget(0)['reflected_norm_lower'], 0)
        with self.assertRaises(ValueError): z.weyl_increment_budget(1, F(1, 4))
        with self.assertRaises(ValueError): z.weyl_increment_budget(1, 2)
        with self.assertRaises(ValueError): z.weyl_increment_budget(0, 0)

    def test_raw_to_phase_bridge_has_exact_fourth_moment_constants_and_pointwise_premise(self):
        row = z.raw_phase_survival_budget(2, 3, 3)
        self.assertEqual(row['raw_fourth_moment_upper'], 432)
        self.assertEqual(row['chosen_phase_parameter_squared'], F(1, 24))
        self.assertEqual(row['Weyl_increment_lower'], F(1, 16))
        self.assertEqual(row['Weyl_reflected_norm_lower'], F(1, 1024))
        self.assertEqual(row['Weyl_energy_upper'], F(1, 8))
        self.assertTrue(row['requires_actual_raw_increment_lower_bound'])
        self.assertTrue(row['pointwise_carre_du_champ_bound_required_not_just_mean_energy'])
        self.assertFalse(row['Gaussian_joint_independence_of_martingales_assumed'])

    def test_phase_bridge_bounds_an_independent_OU_choice_of_phase(self):
        for rate in (F(1, 8), F(1), F(8)):
            raw = z.ou_phase_control(rate)
            variance_floor = raw['raw_increment_variance']['lower']
            bridge = z.raw_phase_survival_budget(1, 1, variance_floor)
            phase = z.ou_phase_control(rate, phase_parameter_squared=bridge['chosen_phase_parameter_squared'])
            self.assertGreaterEqual(phase['Weyl_increment']['lower'], bridge['Weyl_increment_lower'])
            self.assertGreaterEqual(phase['Weyl_reflected_norm']['lower'], bridge['Weyl_reflected_norm_lower'])

    def test_OU_Poisson_correlations_and_reflected_interval_enclose_independent_exponentials(self):
        row = z.ou_phase_control(F(1, 3), 2, F(3, 2), F(2, 5))
        self.assertEqual(row['Poisson_mean'], F(12, 5))
        self.assertEqual(row['phase_energy'], F(4, 5))
        self.assertEqual(row['phase_second_moment'], F(68, 75))
        with localcontext() as context:
            context.prec = 90
            first = (-Decimal('2.4')*(1-Decimal('-0.5').exp())).exp()
            second = (-Decimal('2.4')*(1-Decimal(-1).exp())).exp()
            mass = 1-Decimal('-2.4').exp()
            for key, true in (('phase_correlation_at_T', first), ('phase_correlation_at_2T', second),
                              ('phase_centered_mass', mass), ('Weyl_reflected_norm', 1-2*first+second)):
                self.assertLessEqual(decimal(row[key]['lower']), true)
                self.assertGreaterEqual(decimal(row[key]['upper']), true)

    def test_OU_IR_phase_rescue_coexists_with_linear_filter_collapse(self):
        for n in (2, 32, 256):
            row = z.ou_phase_control(F(1, n))
            self.assertLessEqual(row['raw_linear_reflected_norm']['upper'], F(1, n))
            self.assertGreater(row['Weyl_reflected_norm']['lower'], F(1, 4))
            self.assertEqual(row['phase_energy'], 1)
            self.assertEqual(row['phase_second_moment'], 1+F(1, n))
            self.assertTrue(row['OU_fixture_not_Yang_Mills_vacuum'])
        self.assertGreater(F(33, 65)**2, F(1, 4))

    def test_OU_UV_escape_is_not_repaired_by_bounded_phases(self):
        for n in (2, 32, 256):
            row = z.ou_phase_control(n)
            self.assertLessEqual(row['phase_centered_mass']['upper'], F(1, n))
            self.assertLessEqual(row['Weyl_reflected_norm']['upper'], F(1, n))
            self.assertEqual(row['phase_energy'], 1)
            self.assertEqual(row['phase_second_moment'], 1+n)

    def test_phase_and_OU_controls_reject_false_or_nonpositive_premises(self):
        with self.assertRaises(ValueError): z.raw_phase_survival_budget(1, 1, 3)
        with self.assertRaises(ValueError): z.raw_phase_survival_budget(1, 1, 0)
        with self.assertRaises(ValueError): z.raw_phase_survival_budget(1, 0, 1)
        with self.assertRaises(ValueError): z.ou_phase_control(0)
        with self.assertRaises(ValueError): z.ou_phase_control(1, time=0)
        with self.assertRaises(ValueError): z.ou_phase_control(1, phase_parameter_squared=-1)

    def test_clipped_angle_witness_has_exact_mixed_energy_and_optimized_dual_identity(self):
        for theta in (25, 256, 10**12):
            row = z.clipped_source_budget(theta, F(1, 32), 1, 1, time=2, clock_scale=3)
            self.assertEqual(row['regularizer_eta'], F(1, 16))
            self.assertEqual(row['regularizer_epsilon_fourth_power'], F(1, 16**4*theta))
            self.assertEqual(row['exact_raw_physical_energy'], 3)
            self.assertEqual(row['mixed_raw_regularized_energy_lower'], F(27, 10))
            self.assertEqual(row['regularized_second_moment_upper'], 2**18*9*32)
            coefficient = row['optimized_witness_coefficient']
            denominator = row['witness_denominator_upper']
            lower = row['raw_source_resolvent_norm_squared_lower']
            self.assertEqual(2*2*coefficient*F(27, 10)-coefficient**2*denominator, lower)
            self.assertEqual(lower, row['raw_source_heat_filtered_norm_squared_lower'])
            self.assertGreater(lower, 0)
            epsilon = row['regularizer_epsilon_squared']
            self.assertLessEqual(epsilon['lower']**2, row['regularizer_epsilon_fourth_power'])
            self.assertGreaterEqual(epsilon['upper']**2, row['regularizer_epsilon_fourth_power'])

    def test_clipped_bound_crosses_the_old_bound_but_still_vanishes_with_spacing(self):
        rows = [z.clipped_source_budget(256, F(1, n), 1, 1) for n in (16, 32, 256, 2**24)]
        self.assertFalse(rows[0]['current_lower_bound_is_larger'])
        self.assertTrue(all(row['current_lower_bound_is_larger'] for row in rows[1:]))
        coefficients = []
        for row in rows:
            self.assertEqual(row['fixed_data_a2_asymptotic_coefficient'], F(36, 125))
            self.assertLess(row['normalized_a2_lower_coefficient'], F(36, 125))
            self.assertFalse(row['continuum_positive_lower_bound_proved'])
            coefficients.append(row['normalized_a2_lower_coefficient'])
        self.assertTrue(all(left < right for left, right in zip(coefficients, coefficients[1:])))
        self.assertLess(F(36, 125)-coefficients[-1], F(1, 100))
        self.assertLess(rows[-1]['raw_source_heat_filtered_norm_squared_lower'],
                        rows[0]['raw_source_heat_filtered_norm_squared_lower'])

    def test_clipped_source_budget_rejects_false_cell_and_weak_coupling_premises(self):
        empty = z.clipped_source_budget(25, F(1, 4), 0, 0)
        self.assertEqual(empty['optimized_witness_coefficient'], 0)
        self.assertEqual(empty['raw_source_heat_filtered_norm_squared_lower'], 0)
        with self.assertRaises(ValueError): z.clipped_source_budget(24, F(1, 4), 1, 1)
        with self.assertRaises(ValueError): z.clipped_source_budget(25, F(1, 4), 1, F(1, 32))
        with self.assertRaises(ValueError): z.clipped_source_budget(25, 0, 1, 1)
        with self.assertRaises(ValueError): z.clipped_source_budget(25, F(1, 4), 1, 1, time=0)

    def test_clipped_regularizer_Bell_constants_fit_the_analytic_commutator_majorant(self):
        row = z.clipped_regularizer_constants()
        self.assertEqual(row['analytic_composite_derivative_bounds'], {'2': 9, '3': 76, '4': 984, '5': 17484})
        for degree, partition_count, Bell_number in ((2, 2, 2), (3, 3, 5), (4, 5, 15), (5, 7, 52)):
            terms = row['Bell_partition_terms'][str(degree)]
            self.assertEqual(len(terms), partition_count)
            self.assertEqual(sum(term['Bell_coefficient'] for term in terms), Bell_number)
            self.assertEqual(sum(term['term'] for term in terms), row['analytic_composite_derivative_bounds'][str(degree)])
            for term in terms:
                self.assertEqual(sum(j*m for j, m in enumerate(term['multiplicities'], 1)), degree)
                self.assertEqual(sum(term['multiplicities']), term['outer_derivative_order'])
        self.assertEqual(row['commutator_K'], 952)
        self.assertEqual(row['commutator_Z_upper'], 419280)
        self.assertEqual(row['squared_second_moment_constant_upper'], 27480127880)
        self.assertLess(row['squared_second_moment_constant_upper'], row['squared_constant_majorant'])
        self.assertEqual(row['second_moment_constant']**2, row['squared_constant_majorant'])
        self.assertTrue(row['analytic_derivative_and_commutator_estimates_required'])
        self.assertFalse(row['finite_enumeration_proves_universal_derivative_estimates'])


class LocalWitnessCertificateTests(unittest.TestCase):
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
        self.assertEqual(len(manifest['upstream_sha256']), 288)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertEqual(self.cert['parent_result_sha256'], z.PARENT_HASH)
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_upstream_and_parent_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_uniform_survival_actual_OU_and_continuum_overclaims_are_rejected(self):
        for field in ('fixed_regulator_refinement_implies_uniform_continuum_survival',
                      'actual_uniform_positive_increment_lower_bound_proved',
                      'OU_IR_phase_rescue_is_actual_Yang_Mills_vacuum',
                      'one_time_Weyl_field_regularity_or_local_YM_identification_proved',
                      'full_Euclidean_covariance', 'positive_continuum_mass_gap'):
            wrong = copy.deepcopy(self.cert)
            wrong['model'][field] = True
            with tempfile.TemporaryDirectory() as directory:
                result, pin = Path(directory)/'result', Path(directory)/'pin'
                result.write_text(json.dumps(wrong))
                pin.write_text(z.p.canonical_sha(wrong))
                with self.assertRaises(ValueError): z.check(wrong, result, pin)
                with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_certificate_has_JSON_stable_types_and_no_unintended_control_characters(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                self.assertTrue(all(isinstance(key, str) for key in value))
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        self.assertEqual(json.loads(json.dumps(self.cert)), self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_model_separates_actual_finite_witnesses_from_unproved_uniform_survival(self):
        scope = self.cert['model']
        for field in ('actual_four_cube_local_Wilson_resolvent_witness',
                      'actual_vacuum_error_certified_not_set_to_trial',
                      'actual_untruncated_neutral_hidden_space',
                      'all_coupling_fixed_regulator_source_residual_certificate',
                      'actual_bounded_Weyl_sources_exclude_IR_energy_mass',
                      'ordinary_increment_lower_bound_implies_bounded_phase_survival'):
            self.assertTrue(scope[field])
        for field in ('fixed_regulator_refinement_implies_uniform_continuum_survival',
                      'actual_uniform_positive_increment_lower_bound_proved',
                      'OU_IR_phase_rescue_is_actual_Yang_Mills_vacuum',
                      'full_Euclidean_covariance', 'positive_continuum_mass_gap',
                      'finite_controls_are_universal_proof'):
            self.assertFalse(scope[field])


if __name__ == '__main__':
    unittest.main()
