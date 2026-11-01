import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations, product
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'certificates'))
import ym83_actual_vacuum_response as y


def decimal(x):
    x = F(x)
    return Decimal(x.numerator) / Decimal(x.denominator)


def apply(a, v):
    return [sum((x * z for x, z in zip(row, v)), F(0)) for row in a]


def determinant(a):
    """Independent rational elimination with row pivoting."""
    a = [list(row) for row in a]
    out = F(1)
    for i in range(len(a)):
        pivot = next((j for j in range(i, len(a)) if a[j][i]), None)
        if pivot is None:
            return F(0)
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            out = -out
        value = a[i][i]
        out *= value
        for j in range(i + 1, len(a)):
            ratio = a[j][i] / value
            for k in range(i + 1, len(a)):
                a[j][k] -= ratio * a[i][k]
    return out


def all_principal_minors_nonnegative(a):
    for length in range(1, len(a) + 1):
        for indices in combinations(range(len(a)), length):
            if determinant([[a[i][j] for j in indices] for i in indices]) < 0:
                return False
    return True


def conditional_binary_correlation_upper(states, weights, i, j):
    outside = [k for k in range(len(states[0])) if k not in (i, j)]
    best = F(0)
    for fixed in product((-1, 1), repeat=len(outside)):
        table = {(a, b): sum(w for s, w in zip(states, weights)
                            if (s[i], s[j]) == (a, b) and tuple(s[k] for k in outside) == fixed)
                 for a, b in product((-1, 1), repeat=2)}
        numerator = table[-1, -1] * table[1, 1] - table[-1, 1] * table[1, -1]
        denominator = y.prod(sum(table[a, b] for b in (-1, 1)) for a in (-1, 1))
        denominator *= y.prod(sum(table[a, b] for a in (-1, 1)) for b in (-1, 1))
        best = max(best, y.sqrt_interval(numerator * numerator / denominator)['upper'])
    return best


class ConditionalGapTests(unittest.TestCase):
    def test_native_heat_time_nine_tail_with_independent_decimal_series(self):
        with localcontext() as ctx:
            ctx.prec = 70
            actual = sum(Decimal((n + 1) ** 2) * (-Decimal(9 * n * (n + 2)) / 4).exp()
                         for n in range(1, 70))
            upper = 8 * y.exp_negative_interval(F(27, 4))['upper']
            self.assertLess(actual, decimal(upper))
            self.assertLess(upper, F(1, 20))
            self.assertLess(F(9, 4) * y.exp_negative_interval(F(45, 4))['upper'], F(1, 2))

    def test_rational_sqrt_brackets_include_exact_and_nonsquare_values(self):
        for numerator in range(31):
            for denominator in (1, 3, 7, 11):
                x = F(numerator, denominator)
                result = y.sqrt_interval(x)
                self.assertLessEqual(result['lower'] ** 2, x)
                self.assertGreaterEqual(result['upper'] ** 2, x)
                self.assertLessEqual(result['upper'] - result['lower'], F(1, 10 ** 12))
        self.assertEqual(y.sqrt_interval(F(9, 16)), dict(lower=F(3, 4), upper=F(3, 4)))

    def test_heat_time_optimization_by_exact_completed_square(self):
        # Set omega=3n: the exact optimizer is t=pi. Use any positive P
        # for pi; this checks the symbolic coefficient without floating point.
        for n in (1, 4, 17):
            omega, P = 3 * n, F(22, 7)
            for t in (F(1, 3), F(2), P, F(9)):
                exponent = 2 * t * omega + 6 * P * P * n / t
                minimum = 12 * P * n
                self.assertEqual(exponent - minimum, 2 * omega * (t - P) ** 2 / t)

    def test_original_link_constants_and_free_case(self):
        row = y.original_model_budget(1, 4, 25)['budget']
        self.assertEqual(row['omega'], 200)
        self.assertEqual(row['fixed_time']['exponent'], 3600)
        self.assertEqual(row['optimized']['prefactor'], F(3, 256))
        self.assertEqual(row['optimized']['sqrt_argument'], 600)
        self.assertEqual(y.conditional_budget(9, 0)['exact_free_gap'], F(3, 4))
        self.assertFalse(row['amplitude_uniform'])

    def test_optimized_bound_improves_fixed_time_at_large_amplitude(self):
        row = y.original_model_budget(1, 4, 25)['budget']
        with localcontext() as ctx:
            ctx.prec = 60
            fixed_log = decimal(row['fixed_time']['prefactor']).ln() - decimal(row['fixed_time']['exponent'])
            # Replacing pi and sqrt by upper values only weakens this floor.
            optimized_log = decimal(row['optimized']['prefactor']).ln() - decimal(row['optimized']['exponent_upper'])
            self.assertGreater(optimized_log, fixed_log + 3000)

    def test_symbolic_exponential_floor_stays_positive_without_underflow(self):
        for exponent in (F(0), F(1, 7), F(4), F(308), F(11700)):
            bound = y.symbolic_exp_floor(F(3, 256), exponent)
            self.assertEqual(bound['denominator_power'], exponent.__ceil__())
            with localcontext() as ctx:
                ctx.prec = 60
                self.assertGreaterEqual(Decimal(3).ln() * bound['denominator_power'], decimal(exponent))

    def test_small_original_conditional_bound_exceeds_three_fifths(self):
        x = F(144, 16000)
        floor = F(3, 4) * F(19, 21) ** 2 * (1 - x)
        self.assertGreater(floor, F(3, 5))
        with localcontext() as ctx:
            ctx.prec = 50
            self.assertLess(decimal(1 - x), (-decimal(x)).exp())

    def test_native_incidence_with_reverse_link_enumeration(self):
        for side in (4, 6, 8):
            def plus(v, axis, step):
                return tuple((x + (step if axis == i else 0)) % side for i, x in enumerate(v))
            def faces_at_link(v, axis):
                return {(v, min(axis, other), max(axis, other)) for other in range(3) if other != axis} | {
                    (plus(v, other, -1), min(axis, other), max(axis, other))
                    for other in range(3) if other != axis}
            v = (0, 0, 0)
            boundary = ((v, 0), (plus(v, 0, 1), 1), (plus(v, 1, 1), 0), (v, 1))
            actual = set().union(*(faces_at_link(vertex, axis) for vertex, axis in boundary))
            row = y.native_incidence(side)
            self.assertEqual(row['link_incident_faces'], len(faces_at_link(v, 0)))
            self.assertEqual(row['plaquette_block_incident_faces'], len(actual))
            self.assertEqual(len(actual), 13)
            self.assertEqual(y.original_model_budget(4, len(actual), 25)['budget']['omega'], 650)

    def test_bad_local_parameters_are_rejected(self):
        for value in (True, 0.5, None):
            with self.assertRaises(TypeError):
                y.exact(value)
        for args in ((0, 1), (True, 1), (2, -1)):
            with self.assertRaises(ValueError):
                y.conditional_budget(*args)
        with self.assertRaises(ValueError):
            y.exp_negative_interval(40)


class GlobalResponseTests(unittest.TestCase):
    def test_conditional_projections_are_reversible_orthogonal_projections(self):
        states, weights = y.parity_law(3, F(2, 5))
        for block in ([0], [0, 2], [0, 1, 2]):
            p = y.conditional_projection(states, weights, block)
            self.assertEqual(y.multiply(p, p), p)
            self.assertTrue(all(sum(row) == 1 for row in p))
            self.assertTrue(all(weights[i] * p[i][j] == weights[j] * p[j][i]
                                for i in range(len(states)) for j in range(len(states))))

    def test_two_block_gap_and_anticommutator_are_sharp(self):
        t = F(3, 5)
        states = list(product((-1, 1), repeat=2))
        weights = [(1 + t * x * z) / 4 for x, z in states]
        k = y.resampling_generator(states, weights, [[0], [1]])
        for vector, eigenvalue in (([1] * 4, F(0)), ([x + z for x, z in states], 1 - t),
                                   ([x - z for x, z in states], 1 + t), ([x * z - t for x, z in states], F(2))):
            self.assertEqual(apply(k, vector), [eigenvalue * x for x in vector])
        defect = y.add(y.multiply(k, k), k, -(1 - t))
        weighted = [[weights[i] * x for x in row] for i, row in enumerate(defect)]
        self.assertTrue(all_principal_minors_nonnegative(weighted))
        self.assertEqual(apply(defect, [x + z for x, z in states]), [0] * 4)

    def test_nonparity_three_spin_law_satisfies_full_global_inequality(self):
        states = list(product((-1, 1), repeat=3))
        weights = [(1 + F(1, 10) * a * b + F(1, 12) * b * c + F(1, 15) * a * b * c) / 8
                   for a, b, c in states]
        correlations = [[F(0)] * 3 for _ in range(3)]
        for i, j in combinations(range(3), 2):
            correlations[i][j] = correlations[j][i] = conditional_binary_correlation_upper(states, weights, i, j)
        bound = y.global_budget([1] * 3, correlations)
        self.assertTrue(bound['accepted'])
        k = y.resampling_generator(states, weights, [[0], [1], [2]])
        centered = [[F(i == j) - weights[j] for j in range(8)] for i in range(8)]
        defect = y.add(k, centered, -bound['conditional_bare_gap_lower'])
        weighted = [[weights[i] * x for x in row] for i, row in enumerate(defect)]
        self.assertTrue(all_principal_minors_nonnegative(weighted))

    def test_parity_all_proper_marginals_are_independent_but_conditionals_are_not(self):
        for n in (3, 4):
            t = F(3, 4)
            states, weights = y.parity_law(n, t)
            for length in range(1, n):
                for subset in combinations(range(n), length):
                    for values in product((-1, 1), repeat=length):
                        marginal = sum(w for s, w in zip(states, weights) if tuple(s[i] for i in subset) == values)
                        self.assertEqual(marginal, F(1, 2 ** length))
            for i, j in combinations(range(n), 2):
                self.assertEqual(conditional_binary_correlation_upper(states, weights, i, j), t)

    def test_three_spin_full_matrix_characteristic_polynomial(self):
        for t in (F(0), F(1, 4), F(-3, 4), F(99, 100)):
            states, weights = y.parity_law(3, t)
            k = y.resampling_generator(states, weights, [[0], [1], [2]])
            for value in (F(-1), F(1, 2), F(2), F(4)):
                matrix = [[value * F(i == j) - k[i][j] for j in range(8)] for i in range(8)]
                expected = value * (value - 3) * (value ** 2 - 3 * value + 2 * (1 - t * t)) ** 3
                self.assertEqual(determinant(matrix), expected)

    def test_failed_row_gate_is_inconclusive_even_when_finite_gap_positive(self):
        row = y.parity_control(3, F(3, 4))
        self.assertFalse(row['row_test']['accepted'])
        self.assertIsNone(row['row_test']['conditional_bare_gap_lower'])
        self.assertFalse(row['row_test']['failed_test_proves_gaplessness'])
        self.assertGreater(row['exact_gap_interval']['lower'], 0)

    def test_physical_clock_and_rate_one_block_normalization(self):
        c = [[0, F(1, 4)], [F(1, 4), 0]]
        row = y.global_budget([F(2), F(3)], c, kappa=F(5, 7))
        self.assertEqual(row['conditional_bare_gap_lower'], F(3, 2))
        self.assertEqual(row['conditional_physical_gap_lower'], F(15, 14))
        self.assertFalse(row['actual_vacuum_correlation_hypothesis_discharged'])

    def test_bad_global_and_probability_inputs_are_rejected(self):
        for floors, correlations in (([1, 1], [[0, 1], [0, 0]]), ([0], [[0]]), ([1], [[1]]),
                                      ([1, 1], [[0, 2], [2, 0]]), ([], [])):
            with self.assertRaises(ValueError):
                y.global_budget(floors, correlations)
        with self.assertRaises(ValueError):
            y.conditional_projection([(-1, -1), (1, 1)], [F(1, 2)] * 2, [0])
        states, weights = y.parity_law(3, 0)
        with self.assertRaises(ValueError):
            y.resampling_generator(states, weights, [[0], [0, 2]])
        with self.assertRaises(ValueError):
            y.parity_law(3, 1)


class MarkedTailTests(unittest.TestCase):
    def test_hidden_parity_response_recovered_only_with_tail_or_correct_markers(self):
        rows = y.marker_controls()
        hidden, repaired = rows['invisible'], rows['repaired']
        self.assertEqual(hidden['matrix'], [[0, 0], [0, 0]])
        self.assertEqual(hidden['omitted_squared_norm_upper'], F(9, 16))
        self.assertEqual(hidden['correlation_upper'], F(3, 4))
        self.assertEqual(repaired['matrix'], [[F(3, 4)]])
        self.assertEqual(repaired['omitted_squared_norm_upper'], 0)
        self.assertEqual(repaired['correlation_upper'], F(3, 4))

    def test_rectangular_matrix_norm_and_orthogonal_kernel_tail(self):
        matrix = [[F(3, 10), 0, 0], [0, F(2, 5), 0]]
        row = y.response_budget(matrix, F(2, 5), F(13, 50))
        self.assertEqual(row['captured_squared_norm'], F(1, 4))
        self.assertEqual(row['omitted_squared_norm_upper'], F(1, 100))
        self.assertEqual(row['correlation_upper'], F(1, 2))
        # An independent third singular mode of size 1/10 has exactly
        # this squared tail, while full operator norm is still 2/5.
        self.assertGreaterEqual(row['correlation_upper'], F(2, 5))

    def test_additional_operator_uncertainty_and_correlation_cap(self):
        self.assertEqual(y.response_budget([[F(1, 2)]], F(1, 2), F(1, 4), F(1, 10))['correlation_upper'], F(3, 5))
        self.assertEqual(y.response_budget([[F(3, 4)]], F(3, 4), F(25, 16))['correlation_upper'], 1)

    def test_invalid_finite_norm_or_total_kernel_norm_rejected(self):
        for args in (([[F(1, 2)]], F(1, 4), F(1)), ([[F(1, 2)]], F(1, 2), F(1, 8)),
                     ([[1, 2], [3]], 5, 20), ([], 0, 0)):
            with self.assertRaises(ValueError):
                y.response_budget(*args)

    def test_exact_psd_checker_against_all_principal_minors(self):
        matrices = [
            [[1, 1], [1, 1]], [[0, 1], [1, 0]], [[0, 0], [0, 1]],
            [[2, -1, 0], [-1, 2, -1], [0, -1, 2]],
            [[1, 2, 0], [2, 1, 0], [0, 0, 3]],
        ]
        for matrix in matrices:
            matrix = [[F(x) for x in row] for row in matrix]
            self.assertEqual(y.psd(matrix), all_principal_minors_nonnegative(matrix))


class CertificateTests(unittest.TestCase):
    def test_result_replays_and_has_correct_scope(self):
        fresh = y.run()
        self.assertEqual(y.check(fresh), y.PIN.read_text().strip())
        self.assertEqual(len(fresh['input_sha256']), 359)
        self.assertFalse(fresh['model']['positive_four_dimensional_continuum_mass_gap'])
        self.assertTrue(fresh['model']['actual_conditional_gap_at_every_finite_amplitude'])

    def test_claim_and_source_tampering_are_rejected(self):
        fresh = y.run()
        bad = copy.deepcopy(fresh)
        bad['model']['positive_four_dimensional_continuum_mass_gap'] = True
        with self.assertRaisesRegex(ValueError, 'scope'):
            y.check(bad)
        with tempfile.TemporaryDirectory() as tmp:
            manifest = json.loads(y.SOURCES.read_text())
            key = next(iter(manifest['upstream_sha256']))
            manifest['upstream_sha256'][key] = '0' * 64
            path = Path(tmp) / 'manifest.json'
            path.write_text(json.dumps(manifest))
            with patch.object(y, 'SOURCES', path):
                with self.assertRaisesRegex(ValueError, 'Changed upstream'):
                    y.source_checks()

    def test_cli_check_is_read_only(self):
        before = (y.RESULT.read_bytes(), y.PIN.read_bytes(), y.SOURCES.read_bytes())
        process = subprocess.run([sys.executable, '-B', str(Path(y.__file__)), '--check'], capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertIn('YM83 PASS', process.stdout)
        self.assertEqual(before, (y.RESULT.read_bytes(), y.PIN.read_bytes(), y.SOURCES.read_bytes()))


if __name__ == '__main__':
    unittest.main()
