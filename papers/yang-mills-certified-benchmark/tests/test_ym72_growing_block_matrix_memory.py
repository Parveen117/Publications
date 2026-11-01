import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from math import comb, factorial
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym72_growing_block_matrix_memory as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


class GrowingBlockMatrixMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_result_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_239_upstream_files_and_four_local_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 239)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], '7029c9b1b81a64815e4b7be88aa40c792bb2db04')
        for name in ('YM71_INTERACTING_MEMORY_ENCLOSURES.md', 'YM70_CERTIFIED_SPECTRAL_MEMORY.md'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+name, manifest['upstream_sha256'])

    def test_upstream_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_fixed_coupling_overclaim_rejected(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['fixed_nonzero_coupling_growing_block_accuracy'] = True
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

    def test_proof_and_inputs_have_no_unintended_control_characters(self):
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            raw = (z.ROOT/path).read_bytes()
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in raw), path)

    def test_product_window_and_exact_size_domains(self):
        for theta, factor, count in ((F(1, 10**5), 3, 1), (True, 3, 1), (0.0, 3, 1),
                                     (0, 2, 1), (0, 3.0, 1), (0, True, 1),
                                     (0, 3, 0), (0, 3, True), (0, 3, 1.0)):
            with self.assertRaises(ValueError): z.parameters(theta, factor, count)
        for rho in (-1, True, 0.0, z.MAX_RHO+F(1, 10**20)):
            with self.assertRaises(ValueError): z.envelope(rho)
        self.assertEqual(z.parameters(F(1, 10**6), 3)['rho'], z.MAX_RHO)

    def test_spectral_order_and_clock_domains(self):
        for s in (-1, True, 0.5):
            with self.assertRaises(ValueError): z.matrix_response_interval(0, 3, 1, s)
        for clock in (0, -1, True, 1.0):
            with self.assertRaises(ValueError): z.matrix_response_interval(0, 3, 1, 0, clock)
        for order in (-1, True, 1.0):
            with self.assertRaises(ValueError): z.normalized_moment_interval(F(1, 10**12), 3, 1, order)

    def test_zero_coupling_and_normalization_gate(self):
        for s in (0, 1, 1000):
            row = z.matrix_response_interval(0, 3, 2, s)
            self.assertEqual((row['lower'], row['upper'], row['leading']), (0, 0, 0))
        with self.assertRaises(ValueError): z.scaled_response_interval(0, 3, 2, 0)
        with self.assertRaises(ValueError): z.normalized_moment_interval(0, 3, 2, 1)
        self.assertEqual(z.dynamics_budget(0, 3, 2)['normalized_generator_residual'], 0)

    def test_signed_couplings_preserve_enclosures(self):
        theta = F(1, 10**12)
        self.assertEqual(z.matrix_response_interval(theta, 4, 2, 7),
                         z.matrix_response_interval(-theta, 4, 2, 7))
        self.assertEqual(z.normalized_moment_interval(theta, 4, 2, 1),
                         z.normalized_moment_interval(-theta, 4, 2, 1))

    def test_uniform_gap_floor_with_independent_decimal(self):
        with localcontext() as context:
            context.prec = 90
            alpha = Decimal('0.504')
            true_gap = ((1+1/alpha)/2).ln()/50
            self.assertLess(decimal(z.GAP), true_gap)
            self.assertLessEqual(decimal(z.spatial.parameters(F(1, 10**6))['gamma'].lo), true_gap)

    def test_universal_rounding_against_independent_decimal(self):
        with localcontext() as context:
            context.prec = 90
            pi = Decimal('3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628')
            rho, gap = Decimal('1e-8'), Decimal(1)/125
            K, t0 = (4992*pi*rho).exp(), (Decimal(24)/7).sqrt()
            d0 = Decimal(512)/3*rho
            d1 = d0+768*rho/gap
            e0 = Decimal(16)/9*((240+Decimal(8)/9+8*Decimal(80).sqrt())*rho+7680*rho*rho)
            d, b, e, h = t0*d0, t0*d1, t0*e0, K.sqrt()+t0*d0
            lower = Decimal(13)/9*(K-1+2*K.sqrt()*d+d*d+(2*b+e)*h)
            upper = lower+Decimal(312)/(7*gap)*(d1+e0)**2
            row = z.envelope(F(1, 10**8))
            self.assertGreaterEqual(decimal(row['relative_lower_error']), lower)
            self.assertGreaterEqual(decimal(row['relative_upper_error']), upper)
            self.assertLess(decimal(row['relative_upper_error']), Decimal('0.010608'))

    def test_product_endpoint_and_direct_multiloop_constants(self):
        rows = z.budget_controls()
        for row in rows['direct_vs_universal']:
            self.assertLessEqual(row['direct_errors']['upper_error'], row['universal_error'])
            for s in (F(0), F(1, 7), F(10**9)):
                bound = z.matrix_response_interval(row['theta'], row['block_factor'], row['readouts'], s)
                self.assertGreaterEqual(bound['lower'], F(99, 100)*bound['leading'])
                self.assertLessEqual(bound['upper'], F(51, 50)*bound['leading'])

    def test_no_ambient_volume_in_growing_size_budget(self):
        row = z.parameters(F(1, 10**30), 10**8, 10**10)
        self.assertNotIn('volume', row)
        self.assertLess(row['relative_upper_error'], F(1, 1000))

    def test_single_readout_reduces_to_previous_constants(self):
        theta = F(1, 10**12)
        old = z.previous.parameters(theta, 5)
        new = z.matrix_from_scalar_budget(theta, 5, 1)
        self.assertEqual(new['lower_error'], old['relative_lower_error'])
        self.assertEqual(new['upper_error'], old['relative_upper_error'])

    def test_native_stars_are_disjoint_and_keep_hidden_links(self):
        evidence = z.geometry_controls()
        self.assertTrue(evidence['native_full_lattice_comparison'])
        self.assertEqual(len(evidence['growing_families']), 9)
        for row in evidence['growing_families']:
            self.assertTrue(row['disjoint'])
            self.assertLessEqual(row['each_star_links'], 52*row['block_factor']-16)

    def test_periodic_geometry_translation_and_coarse_origin_gate(self):
        self.assertEqual(z.star_geometry(3, (0, 0, 0), 12), z.star_geometry(3, (12, -12, 24), 12))
        original = z.star_geometry(3, (0, 0, 0), 12)['support']
        shifted = {(((v[0]+9) % 12, v[1], v[2]), axis) for v, axis in original}
        self.assertEqual(shifted, z.star_geometry(3, (9, 0, 0), 12)['support'])
        for factor, origin, side in ((3, (1, 0, 0), 12), (3, (0, 0), 12),
                                     (3, (False, 0, 0), 12), (3, (0, 0, 0), 9), (3, (0, 0, 0), 13)):
            with self.assertRaises(ValueError): z.star_geometry(factor, origin, side)

    def test_noncommuting_matrix_memory_satisfies_quadratic_commutator_bound(self):
        evidence = z.matrix_controls()
        C = evidence['nonzero_normalized_commutator']
        self.assertEqual(C[0][0], 0)
        self.assertEqual(C[1][0], -C[0][1])
        self.assertNotEqual(C[0][1], 0)
        self.assertLessEqual(abs(C[0][1]), 2*evidence['relative_enclosure_error']**2)

    def test_matrix_moments_need_no_commuting_spectral_weights(self):
        for order in range(5):
            step = F(1, 9)
            total = z.a.scale(z.a.identity(2), 0)
            for k in range(order+1):
                total = z.a.add(total, z.spectral_fixture(k*step), (-1)**k*comb(order, k))
            total = z.a.scale(total, 1/(factorial(order)*step**order))
            moment = z.spectral_fixture(0, order)
            correction = z.previous.product(1+k*step/3 for k in range(1, order+1))
            self.assertTrue(z.a.psd(z.a.add(moment, total, -1)))
            self.assertTrue(z.a.psd(z.a.add(z.a.scale(total, correction), moment, -1)))

    def test_alternating_stencil_requires_opposite_interval_endpoints(self):
        intervals = [(F(9, 10), F(11, 10)), (F(2, 5), F(3, 5))]
        row = z.previous.stencil_interval(intervals, 1, 1)
        self.assertEqual(row['divided_difference_lower'], F(3, 10))
        self.assertEqual(row['divided_difference_upper'], F(7, 10))
        wrong_lower = intervals[0][0]-intervals[1][0]
        self.assertGreater(wrong_lower, row['divided_difference_lower'])

    def test_exact_profile_enclosures_include_all_four_rates(self):
        rows = z.profile_controls()['exact_free_profile_checks']
        self.assertEqual(len(rows), 20)
        for row in rows:
            r, value, phi = row['block_factor'], row['exact_profile'], row['limit']
            self.assertGreaterEqual(value, F(6*r+3, 6*r+8)*phi)
            self.assertLessEqual(value, (1+F(1, 2*r))*phi)

    def test_first_correction_sign_and_second_order_remainder(self):
        self.assertEqual(z.finite_block_correction(0), -F(1, 216))
        self.assertEqual(z.finite_block_correction(1), 0)
        self.assertGreater(z.finite_block_correction(2), 0)
        for row in z.profile_controls()['exact_free_profile_checks']:
            r = row['block_factor']
            self.assertLessEqual(abs(row['remainder']), row['limit']*F(43, 36*r*r))

    def test_moment_rescaling_preserves_each_inverse_power(self):
        for r in (3, 16, 100):
            atoms = z.leading.leading_atoms(r)
            scaled = tuple((rate/r, weight/r) for rate, weight in atoms)
            for order in range(5):
                base = z.leading.atom_moments(atoms, order+1)[order][0][0]
                observed = z.leading.atom_moments(scaled, order+1)[order][0][0]
                self.assertEqual(observed, r**order*base)

    def test_actual_moment_intervals_narrow_along_stronger_budget(self):
        rows = z.budget_controls()['actual_matrix_moment_intervals']
        for order in (1, 2):
            errors = []
            for row in rows:
                if row['order'] != order: continue
                self.assertLessEqual(row['lower'], row['leading'])
                self.assertGreaterEqual(row['upper'], row['leading'])
                self.assertLess((row['upper']-row['lower'])/row['leading'], F(1, 1000))
                errors.append(max(abs(row[k]-row['limiting_moment']) for k in ('lower', 'upper')))
            self.assertEqual(len(errors), 3)
            self.assertGreater(errors[0], errors[1])
            self.assertGreater(errors[1], errors[2])

    def test_response_and_higher_moments_have_different_scaling_budgets(self):
        for r in (4, 16, 64):
            N, theta = r, F(1, 10**8*r**3)
            self.assertEqual(z.parameters(theta, r, N)['rho'], F(1, 10**8*r))
            self.assertEqual(abs(theta)*N*r*r, F(1, 10**8))
        with self.assertRaises(ValueError): z.parameters(F(1, 10**8), 301, 1)

    def test_uniform_response_convergence_can_miss_a_higher_moment(self):
        rows = z.counterexample_controls()['higher_moment_budget_is_necessary_for_an_inference_from_these_bounds']
        for row in rows:
            n, r = row['target_order'], row['block_factor']
            self.assertEqual(row['relative_response_error'], F(1, r**n))
            self.assertEqual(row['first_uncontrolled_moment'],
                             (1-F(1, r**n))*F(1, 12*3**(n+1))+F(1, 36*z.GAP**n))
            for j in range(n):
                excess = F(1, 36*z.GAP**j*r**(n-j))
                self.assertLessEqual(excess, F(1, 36*z.GAP**j*r))

    def test_block_clock_memory_amplitude_is_not_normalized_profile(self):
        theta, r, N, spectral = F(1, 10**12), 4, 2, F(7, 3)
        clock = z.matrix_response_interval(theta, r, N, spectral, F(1, r))
        normalized = z.scaled_response_interval(theta, r, N, spectral)
        for field in ('lower', 'upper', 'leading'):
            self.assertEqual(clock[field], theta*theta/r*normalized[field])

    def test_actual_dynamics_budget_has_positive_gram_and_small_time_error(self):
        row = z.dynamics_budget(F(1, 3*10**8), 3, 1)
        self.assertGreater(row['centered_Gram_lower'], 0)
        self.assertLess(row['centered_Gram_lower'], F(1, 4))
        self.assertGreater(row['centered_Gram_upper'], F(1, 4))
        self.assertLess(row['all_nonnegative_times_correlation_error'], F(2134, 10**10))
        self.assertGreaterEqual(row['spectral_weight_below_three_halves'],
                                F(4, 9)*row['normalized_generator_residual']**2)

    def test_visible_relaxation_can_miss_the_full_small_gap(self):
        rows = z.counterexample_controls()['visible_rate_does_not_determine_full_gap']
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertTrue(row['zero_tilt_gives_exact_rate_three_at_the_same_full_gap'])
            self.assertGreater(row['observed_correlation'], row['rate_three_correlation'])
            self.assertLessEqual(abs(row['observed_correlation']-row['rate_three_correlation']),
                                 row['generator_residual']/3)
        for r in (3, 100, 10**6):
            budget = z.dynamics_budget(F(1, 10**8*r), r)
            self.assertEqual(budget['full_physical_gap_upper_after_block_clock'], F(100, 19*r))

    def test_scope_preserves_physical_continuum_and_slow_mode_obligations(self):
        claim = self.cert['model']
        for field in ('matrix_memory_for_separated_Wilson_readouts', 'uniform_in_admissible_finite_volume',
                      'joint_normalized_one_pole_memory_limit', 'actual_matrix_moments_enclosed',
                      'actual_normalized_loop_time_correlations_enclosed',
                      'full_gap_vanishes_on_the_declared_block_clock'):
            self.assertTrue(claim[field])
        for field in ('fixed_nonzero_coupling_growing_block_accuracy',
                      'unnormalized_memory_has_nonzero_limit_on_block_clock',
                      'full_spatial_memory_distance_decay', 'infinite_order_uniformity',
                      'larger_bare_coupling_window', 'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(claim[field])


if __name__ == '__main__':
    unittest.main()
