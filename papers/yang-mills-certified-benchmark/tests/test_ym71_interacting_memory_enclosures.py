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
import ym71_interacting_memory_enclosures as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class InteractingMemoryEnclosureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_result_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_232_upstream_files_and_four_local_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 232)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], '71d588f2bde2b572c65a36bf9760fbd89eb3ea8d')
        for name in ('YM70_CERTIFIED_SPECTRAL_MEMORY.md', 'YM67_POINTWISE_COARSE_VACUUM.md'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+name, manifest['upstream_sha256'])

    def test_upstream_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_growing_block_overclaim_rejected(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['uniform_as_block_factor_grows'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_no_binary_floats_in_evidence(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)

    def test_exact_hidden_lifts_and_raw_hiddenness_failure(self):
        rows = z.native_controls()['conditional_lifts']
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row['weighted_lift_has_exact_zero_conditional_mean'] for row in rows))
        self.assertTrue(all(row['raw_free_component_has_nonzero_conditional_mean'] for row in rows))

    def test_lift_must_use_weighted_factor_two(self):
        g = z.model.loop_polynomial()
        h = z.p.add(z.model.source.constant(1, 4), z.p.scale(g, F(1, 7)))
        gamma = z.p.add(*(z.p.mul(z.native.derivative(h, 0, axis), z.native.derivative(g, 0, axis))
                          for axis in (1, 2, 3)))
        wrong_numerator = z.p.add(z.p.mul(h, g), z.p.scale(gamma, -F(4, 3)))
        self.assertNotEqual(z.integrate_link(z.p.mul(h, wrong_numerator), 0), {})

    def test_hidden_link_is_absent_from_all_retained_paths(self):
        rows = z.native_controls()['actual_product_cut_geometry']
        for row in rows:
            r = row['block_factor']
            self.assertTrue(row['every_component_has_an_unretained_fundamental_link'])
            self.assertLessEqual(row['actual_union_links'], 52*r-16)
            self.assertEqual(row['incident_plaquettes'], 16*r-4)
            self.assertEqual(row['absolute_word_weight'], F(16*r, 3))

    def test_square_roots_are_exact_upper_bounds(self):
        for value in (F(0), F(1, 4), F(7, 24), F(2), F(1, 10**40)):
            root = z.sqrt_upper(value)
            self.assertGreaterEqual(root*root, value)
            if root:
                self.assertLess((root-1/F(z.GRID))**2, value)
        for value in (-1, True, 0.5):
            with self.assertRaises(ValueError): z.sqrt_upper(value)

    def test_exponential_rounding_with_independent_decimal(self):
        with localcontext() as context:
            context.prec = 90
            for value in (F(0), F(1, 10000), F(1, 3), F(2), F(17)):
                self.assertGreaterEqual(decimal(z.exponential_upper(value)), decimal(value).exp())

    def test_gap_lower_bound_with_independent_decimal(self):
        with localcontext() as context:
            context.prec = 90
            for theta in (F(0), F(1, 10**8), F(1, 16000)):
                alpha = Decimal(1)/2+4000*decimal(theta)
                exact_gap = ((1+1/alpha)/2).ln()/50
                self.assertLessEqual(decimal(z.parameters(theta)['gap_lower']), exact_gap)
                self.assertGreater(z.parameters(theta)['gap_lower'], 0)

    def test_full_window_domain_is_not_extended(self):
        for theta in (F(1, 8000), -F(1, 8000), F(1), True, 0.00001):
            with self.assertRaises(ValueError): z.parameters(theta)
        for factor in (0, 1, 2, True, 3.0):
            with self.assertRaises(ValueError): z.parameters(F(1, 10**8), factor)

    def test_zero_coupling_closes_actual_response_and_moments(self):
        for s in (0, 1, 1000):
            row = z.response_interval(0, 3, s)
            self.assertEqual(row, dict(lower=F(0), upper=F(0), leading=F(0)))
        for n in range(4):
            row = z.moment_interval(0, 3, n)
            self.assertEqual((row['lower'], row['upper'], row['leading']), (0, 0, 0))

    def test_signed_couplings_have_same_certified_budgets(self):
        theta = F(1, 10**8)
        for s in (0, F(1, 7), 1000):
            self.assertEqual(z.response_interval(theta, 3, s), z.response_interval(-theta, 3, s))
        self.assertEqual(z.moment_interval(theta, 3, 1), z.moment_interval(-theta, 3, 1))

    def test_explicit_nonzero_window(self):
        for theta in (F(1, 10**8), F(1, 10**9), F(1, 10**12)):
            cell = z.parameters(theta)
            self.assertLess(cell['relative_lower_error'], F(16, 1000))
            self.assertLess(cell['relative_upper_error'], F(61, 1000))
            for s in (F(0), F(1, 7), F(1000), F(10**12)):
                row = z.response_interval(theta, 3, s)
                self.assertGreater(row['lower'], 0)
                self.assertGreaterEqual(row['lower'], F(49, 50)*row['leading'])
                self.assertLessEqual(row['upper'], F(107, 100)*row['leading'])

    def test_broad_window_does_not_have_small_relative_error(self):
        cell = z.parameters(F(1, 16000))
        self.assertGreater(cell['relative_lower_error'], 1)
        self.assertGreater(cell['relative_upper_error'], 1000)
        self.assertEqual(z.response_interval(F(1, 16000), 3, 0)['lower'], 0)

    def test_no_hidden_volume_parameter_and_visible_block_dependence(self):
        small, larger = z.parameters(F(1, 10**8), 3), z.parameters(F(1, 10**8), 5)
        self.assertGreater(larger['union_support_upper'], small['union_support_upper'])
        self.assertGreater(larger['lift_error'], small['lift_error'])
        self.assertGreater(larger['eigen_equation_residual'], small['eigen_equation_residual'])
        self.assertNotIn('volume', small)

    def test_residual_completion_and_inverse_power_telescopes(self):
        evidence = z.residual_controls()
        self.assertEqual(evidence['inverse_power_telescope_checks'], 30)
        for row in evidence['energy_residual_identities']:
            self.assertEqual(row['exact_response']-row['variational_value'], row['residual_correction'])
            self.assertGreaterEqual(row['residual_correction'], 0)

    def test_stieltjes_stencil_factorial_and_support_correction(self):
        rate, weight, step = F(3), F(2, 5), F(1, 7)
        for order in range(6):
            values = [weight/(rate+k*step) for k in range(order+1)]
            row = z.stencil_interval([(v, v) for v in values], step, rate)
            true_moment = weight/rate**(order+1)
            self.assertEqual(row['upper'], true_moment)
            self.assertLessEqual(row['lower'], true_moment)

    def test_signed_endpoint_arithmetic_for_uncertain_responses(self):
        row = z.stencil_interval([(F(9, 10), F(11, 10)), (F(2, 5), F(3, 5))], 1, 1)
        self.assertEqual(row['divided_difference_lower'], F(3, 10))
        self.assertEqual(row['divided_difference_upper'], F(7, 10))
        self.assertEqual(row['upper'], F(7, 5))

    def test_stencils_refuse_invalid_or_contradictory_inputs(self):
        for values, step, gap in (([], 1, 1), ([(2, 1)], 1, 1), ([(1, 1)], 0, 1),
                                  ([(1, 1)], 1, 0), ([(0, 0), (1, 1)], 1, 1),
                                  ([(1.0, 1.0)], 1, 1)):
            with self.assertRaises(ValueError): z.stencil_interval(values, step, gap)

    def test_response_bounds_cannot_be_differentiated(self):
        row = z.stencil_controls()['response_bounds_cannot_be_differentiated']
        self.assertEqual(row['relative_value_width'], F(1, 100))
        self.assertEqual(row['baseline_first_moment'], 1)
        self.assertEqual(row['alternative_first_moment'], F(10099, 100))

    def test_displayed_interacting_moment_enclosures_are_outward(self):
        for exponent, order, lower, upper in ((12, 1, F(996, 1000), F(1133, 1000)),
                                            (16, 2, F(996, 1000), F(1226, 1000)),
                                            (20, 3, F(997, 1000), F(1029, 1000))):
            row = z.moment_interval(F(1, 10**exponent), 3, order)
            self.assertGreaterEqual(row['lower'], lower*row['leading'])
            self.assertLessEqual(row['upper'], upper*row['leading'])

    def test_moment_intersection_keeps_direct_radius_and_gap_budget(self):
        theta = F(1, 10**12)
        for order in range(4):
            row = z.moment_interval(theta, 3, order)
            self.assertGreaterEqual(row['lower'], max(0, row['leading']-row['direct_radius']))
            self.assertLessEqual(row['upper'], row['leading']+row['direct_radius'])
            self.assertLessEqual(row['upper'], z.response_interval(theta, 3, 0)['upper']
                                 /z.parameters(theta)['gap_lower']**order)

    def test_clock_covariance_changes_each_moment_dimension(self):
        theta, clock = F(1, 10**12), F(7, 3)
        for order in range(4):
            scaled, base = z.moment_interval(theta, 3, order, clock), z.moment_interval(theta, 3, order)
            for field in ('lower', 'upper', 'leading'):
                self.assertEqual(scaled[field], clock**(1-order)*base[field])
        scaled, base = z.response_interval(theta, 3, 5, clock), z.response_interval(theta, 3, 5/clock)
        for field in ('lower', 'upper', 'leading'):
            self.assertEqual(scaled[field], clock*base[field])

    def test_spectral_parameters_orders_and_clocks_reject_inexact_inputs(self):
        for s in (-1, True, 0.5):
            with self.assertRaises(ValueError): z.response_interval(0, 3, s)
        for order in (-1, True, 1.0):
            with self.assertRaises(ValueError): z.moment_interval(0, 3, order)
        for clock in (0, -1, True, 1.0):
            with self.assertRaises(ValueError): z.moment_interval(0, 3, 0, clock)

    def test_scope_preserves_continuum_and_locality_obligations(self):
        claim = self.cert['model']
        for field in ('actual_scalar_interacting_response_enclosed', 'actual_scalar_inverse_moments_enclosed',
                      'volume_uniform_fixed_block_leading_remainder', 'exact_interacting_hidden_trial_vectors'):
            self.assertTrue(claim[field])
        for field in ('full_vacuum_numerically_computed', 'exactly_four_rates_for_full_memory',
                      'stable_interval_Gram_inversion_proved', 'full_spatial_memory_locality',
                      'uniform_as_block_factor_grows', 'enlarged_bare_coupling_window',
                      'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(claim[field])


if __name__ == '__main__':
    unittest.main()
