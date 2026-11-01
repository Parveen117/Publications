import copy
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb, cos, pi, sin, sqrt
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'certificates'))
import ym84_gauge_resolved_scaling as y


def independent_character(n):
    out = [F(0)] * (n + 1)
    for k in range(n // 2 + 1):
        out[n - 2 * k] = F((-1) ** k * comb(n - k, k) * 2 ** (n - 2 * k))
    return out


def derivative(p):
    return [i * p[i] for i in range(1, len(p))] or [F(0)]


def matrix_rank(matrix):
    rows = [list(row) for row in matrix]
    rank = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[pivot], rows[rank] = rows[rank], rows[pivot]
        value = rows[rank][col]
        rows[rank] = [x / value for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            factor = rows[i][col]
            rows[i] = [a - factor * b for a, b in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


class CharacterResponseTests(unittest.TestCase):
    def test_haar_moments_against_beta_integral_recurrence(self):
        value = F(1)
        for k in range(24):
            self.assertEqual(y.haar_moment(2 * k), value)
            self.assertEqual(y.haar_moment(2 * k + 1), 0)
            value *= F(2 * k + 1, 2 * k + 4)

    def test_character_recurrence_against_closed_binomial_formula(self):
        for n in range(25):
            self.assertEqual(y.character(n), independent_character(n))
            self.assertEqual(y.poly_value(y.character(n), 1), n + 1)

    def test_haar_character_orthogonality_exact(self):
        for n in range(14):
            for m in range(14):
                self.assertEqual(y.haar_integral(y.poly_multiply(y.character(n), y.character(m))), F(n == m))

    def test_four_link_radial_Casimir_normalization(self):
        for n in range(18):
            p = y.character(n)
            radial = y.poly_add(y.poly_multiply([1, 0, -1], derivative(derivative(p))),
                                y.poly_multiply([0, 3], derivative(p)), -1)
            self.assertEqual(radial, y.polynomial([-n * (n + 2) * x for x in p]))
        self.assertEqual(1 * (1 + 2), 4 * F(3, 4))

    def test_character_defect_constant_from_weights_and_derivatives(self):
        for n in range(1, 50):
            weights = range(-n, n + 1, 2)
            expected = F(n * (n + 2), 3)
            self.assertEqual(F(sum(k * k for k in weights), n + 1), expected)
            self.assertEqual(y.poly_value(derivative(y.character(n)), 1) / (n + 1), expected)

    def test_uniform_density_has_no_centered_response(self):
        out = y.central_response([7])
        self.assertEqual(out['normalized_density'], [1])
        self.assertEqual(out['correlation'], 0)
        self.assertEqual(out['delta'], 1)
        self.assertEqual(out['total_hs_squared'], 0)

    def test_exact_linear_density_response_includes_representation_multiplicity(self):
        # q=1+aW: lambda_1=a/4, with four matrix coefficients.
        a = F(3, 4)
        out = y.central_response([1, a], 1)
        self.assertEqual(out['lambda_by_twice_spin'], {1: a / 4})
        self.assertEqual(out['total_hs_squared'], a * a / 4)
        self.assertEqual(out['omitted_hs_squared'], 0)

    def test_response_reconstructs_polynomial_and_exact_omitted_tail(self):
        for k in (1, 2, 5, 9):
            out = y.concentrated_control(k, F(2, 3), 2)
            reconstructed = [F(1)]
            for n, value in out['lambda_by_twice_spin'].items():
                reconstructed = y.poly_add(reconstructed, independent_character(n), (n + 1) * value)
            self.assertEqual(reconstructed, out['normalized_density'])
            residual = list(reconstructed)
            residual = y.poly_add(residual, [1], -1)
            for n in (1, 2):
                residual = y.poly_add(residual, independent_character(n),
                                      -(n + 1) * out['lambda_by_twice_spin'].get(n, F(0)))
            self.assertEqual(y.haar_integral(y.poly_multiply(residual, residual)), out['omitted_hs_squared'])
            self.assertGreater(out['density_lower'], 0)

    def test_fundamental_channel_does_not_bound_higher_response(self):
        for n in (2, 3, 6, 9):
            out = y.hidden_channel_control(n)
            self.assertEqual(out['delta'], 1)
            self.assertEqual(out['lambda_by_twice_spin'][1], 0)
            self.assertGreater(out['correlation'], 0)
            self.assertEqual(out['captured_hs_squared'], 0)
            self.assertEqual(out['total_hs_squared'], out['omitted_hs_squared'])

    def test_comparison_concentration_activates_the_proven_tail_bound(self):
        out = y.concentrated_control(12)
        n = 1
        self.assertLessEqual(n * (n + 2) * out['delta'], F(3, 2))
        lower = y.tail_budget(n, out['delta'], 0)['omitted_hs_squared_lower']
        self.assertGreaterEqual(out['total_hs_squared'], lower)

    def test_rank_dimension_and_actual_fourth_power_constants(self):
        for n in range(100):
            self.assertEqual(y.channel_dimension(n), sum((j + 1) ** 2 for j in range(1, n + 1)))
        for numerator in range(16, 513):
            s = F(numerator, 4)
            out = y.actual_fourth_power_budget(s)
            self.assertGreaterEqual(F(out['dimension'], 4), out['xi_lower'])
            n = out['cutoff']
            self.assertLessEqual(n * (n + 2) * out['delta_upper'], F(3, 2))
            self.assertEqual(out['theta'], s ** 4)

    def test_rank_budget_boundary_and_invalid_inputs(self):
        self.assertEqual(y.tail_budget(2, F(3, 16), 8)['omitted_hs_squared_lower'], F(5, 4))
        self.assertEqual(y.tail_budget(2, F(3, 16), 20)['omitted_hs_squared_lower'], 0)
        for fn, args in ((y.tail_budget, (2, F(1, 5), 0)),
                         (y.actual_fourth_power_budget, (F(7, 2),)),
                         (y.concentrated_control, (2, 1)),
                         (y.concentrated_control, (2, -1)),
                         (y.hidden_channel_control, (1,)),
                         (y.haar_moment, (-1,)),
                         (y.central_response, ([0],))):
            with self.assertRaises(ValueError):
                fn(*args)


class GaugeAndScaleTests(unittest.TestCase):
    def test_gauge_projection_is_weighted_orthogonal_with_two_orbits(self):
        for t in (F(0), F(3, 4), F(-4, 5)):
            data = y.cycle_matrices(t)
            g, w = data['gauge'], data['weights']
            self.assertEqual(matrix_rank(g), 2)
            self.assertEqual(y.prev.multiply(g, g), g)
            for i in range(16):
                self.assertEqual(sum(g[i]), 1)
                for j in range(16):
                    self.assertEqual(w[i] * g[i][j], w[j] * g[j][i])
                    self.assertEqual(g[i][j] > 0,
                                     y.prev.prod(data['states'][i]) == y.prev.prod(data['states'][j]))

    def test_all_sixteen_resampling_eigenvectors_and_exact_gap(self):
        for t in (F(0), F(1, 3), F(-7, 8)):
            data = y.cycle_matrices(t)
            states, k = data['states'], data['k']
            left = [s[0] * s[1] for s in states]
            right = [s[2] * s[3] for s in states]
            modes = [(F(0), [F(1)] * 16), (F(2), [y.prev.prod(s) - t for s in states])]
            modes += [(1 - sign * t, [a + sign * b for a, b in zip(left, right)]) for sign in (-1, 1)]
            modes += [(F(1), [F(s[i]) for s in states]) for i in range(4)]
            modes += [(F(2), [F(s[i] * s[j]) for s in states]) for i in (0, 1) for j in (2, 3)]
            for triple in combinations(range(4), 3):
                missing = next(i for i in range(4) if i not in triple)
                modes.append((F(2), [y.prev.prod(s[i] for i in triple) - t * s[missing] for s in states]))
            self.assertEqual(matrix_rank([vector for _, vector in modes]), 16)
            for value, vector in modes:
                self.assertEqual(y.apply(k, vector), [value * x for x in vector])
            self.assertEqual(min(value for value, _ in modes if value), 1 - abs(t))

    def test_small_full_response_gap_is_removed_by_Gauss_projection(self):
        for t in (F(1, 2), F(99, 100), F(9999, 10000)):
            out = y.cycle_control(t)
            self.assertEqual(out['physical_resampling_gap'], 2)
            self.assertEqual(out['full_resampling_gap'], 1 - t)
            self.assertTrue(out['physical_pair_response_zero'])
            self.assertEqual(out['gauge_projection_rank'], 2)

    def test_restricted_projection_identity_keeps_the_physical_loop(self):
        data = y.cycle_matrices(F(4, 5))
        physical = [3 + 7 * y.prev.prod(s) for s in data['states']]
        self.assertEqual(y.apply(data['gauge'], physical), physical)
        mean = sum(w * f for w, f in zip(data['weights'], physical))
        for key in ('p1', 'p2'):
            self.assertEqual(y.apply(data[key], physical), [mean] * 16)
        self.assertEqual(y.apply(data['k'], physical), [2 * (f - mean) for f in physical])

    def test_radial_oscillator_constants_have_strict_rational_slack(self):
        out = y.oscillator_constants()
        self.assertEqual(out['ground_constant_squared_upper'], F(821, 147))
        self.assertEqual(out['gap_coefficient_lower'], F(6749, 8712))
        self.assertEqual(out['gap_lower_slack'], F(215, 8712))
        self.assertEqual(out['gap_coefficient_upper'], F(440, 49))
        self.assertGreater(7 * sqrt(2) / pi - sqrt(2 * pi * pi / 3 - 1), float(out['gap_coefficient_lower']))

    def test_sine_trial_second_moment_independent_quadrature(self):
        # Numerical diagnostic of the analytic integral, not a rigorous interval claim.
        n = 4000
        integrand = lambda x: 2 * x * x * sin(pi * x) ** 2
        integral = (integrand(0) + integrand(1)
                    + sum((4 if i % 2 else 2) * integrand(i / n) for i in range(1, n))) / (3 * n)
        self.assertAlmostEqual(integral, 1 / 3 - 1 / (2 * pi * pi), places=12)
        for i in range(1, 1000):
            r = pi * i / 1000
            self.assertLessEqual(2 * r * r / (pi * pi), 1 - cos(r) + 1e-14)
            self.assertLessEqual(1 - cos(r), r * r / 2 + 1e-14)

    def test_trial_lengths_are_admissible_for_theta_at_least_one(self):
        a_constant = 1 / 6 - 1 / (4 * pi * pi)
        for theta in (1, 2, 25, 256, 10 ** 8):
            b0_squared = pi / sqrt(theta * a_constant)
            b1_squared = 2 * sqrt(2) * pi / sqrt(theta)
            self.assertLessEqual(b0_squared, pi * pi)
            self.assertLessEqual(b1_squared, pi * pi)
            e0 = pi * pi / b0_squared + theta * a_constant * b0_squared
            e1 = 4 * pi * pi / b1_squared + theta * b1_squared / 2
            self.assertAlmostEqual(e0 / sqrt(theta), sqrt(2 * pi * pi / 3 - 1), places=12)
            self.assertAlmostEqual(e1 / sqrt(theta), 2 * sqrt(2) * pi, places=12)

    def test_scale_matched_margin_and_constant_margin_have_distinct_limits(self):
        for denominator in (4, 16, 256, 4096):
            a = F(1, denominator)
            matched = y.scale_budget(a, 2, F(3, 8), a, 4, 1)
            constant = y.scale_budget(a, 2, F(3, 8), F(1, 2), 4, 1)
            self.assertEqual(matched['conditional_gap_lower'], F(3, 4))
            self.assertTrue(matched['hypotheses_compatible'])
            self.assertEqual(constant['conditional_gap_lower'], F(3 * denominator, 8))
            self.assertEqual(constant['hypotheses_compatible'], denominator <= F(32, 3))

    def test_physical_clock_cancels_gauge_coupling_exactly(self):
        # cE=cB=1, so sqrt(theta)=1/g^2 and kappa=g^2/a.
        for g in (F(1, 2), F(2, 3), F(1, 100)):
            for a in (F(1, 8), F(1, 1000)):
                theta, kappa = 1 / g ** 4, g ** 2 / a
                self.assertEqual(kappa ** 2 * theta, 1 / a ** 2)
                self.assertEqual(kappa * F(3, 4) / g ** 2, F(3, 4) / a)

    def test_scale_inputs_reject_invalid_margins_and_variances(self):
        for values in ((0, 1, 1, 1, 1, 1), (1, 1, 1, 2, 1, 1),
                       (1, 1, 1, 1, 1, 0), (1, 1, 1, -1, 1, 1)):
            with self.assertRaises(ValueError):
                y.scale_budget(*values)


class EvidenceScopeTests(unittest.TestCase):
    def test_source_chain_and_written_result_reproduce(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']), 366)
        self.assertEqual(y.check(cert), y.PIN.read_text().strip())
        self.assertTrue(cert['model']['actual_homogeneous_vacuum_open_path_response_identity'])
        for key in ('full_lattice_physical_local_sqrt_theta_floor_proved',
                    'full_lattice_physical_response_margin_of_order_a_proved',
                    'nonzero_continuum_state_survival_proved',
                    'positive_four_dimensional_continuum_mass_gap'):
            self.assertFalse(cert['model'][key])

    def test_modified_result_pin_and_scope_are_rejected(self):
        cert = y.run()
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory) / 'result.json', Path(directory) / 'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0' * 64)
            with self.assertRaises(ValueError):
                y.check(cert, result, pin)
            pin.write_text(y.prev.canonical_sha(cert))
            changed = copy.deepcopy(cert)
            changed['controls']['actual_rank_budgets'][0]['dimension'] = 0
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert, result, pin)
            changed = copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
            with self.assertRaises(ValueError):
                y.check(changed, result, pin)

    def test_upstream_source_drift_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        path = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][path] = '0' * 64
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'manifest.json'
            source.write_text(json.dumps(manifest))
            with patch.object(y, 'SOURCES', source):
                with self.assertRaisesRegex(ValueError, 'Changed upstream input'):
                    y.source_checks()

    def test_read_only_check_and_optimization_guard(self):
        before = [path.read_bytes() for path in (y.RESULT, y.PIN, y.SOURCES)]
        checked = subprocess.run([sys.executable, '-B', str(y.HERE / 'ym84_gauge_resolved_scaling.py'), '--check'],
                                 capture_output=True, text=True)
        self.assertEqual(checked.returncode, 0, checked.stderr)
        self.assertIn('YM84 PASS', checked.stdout)
        self.assertEqual(before, [path.read_bytes() for path in (y.RESULT, y.PIN, y.SOURCES)])
        optimized = subprocess.run([sys.executable, '-B', '-O', str(y.HERE / 'ym84_gauge_resolved_scaling.py'), '--check'],
                                   capture_output=True, text=True)
        self.assertNotEqual(optimized.returncode, 0)
        self.assertIn('without optimization', optimized.stderr)


if __name__ == '__main__':
    unittest.main()
