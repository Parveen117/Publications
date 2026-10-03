"""Finite source, identifiability and tamper controls for GE4."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ge4_observer_recovery as g


def pair_moments(count=8):
    return tuple((F(1, 2)**k+F(1)**k)/2 for k in range(count+1))


class ObserverRecoveryTests(unittest.TestCase):
    def test_native_polynomials_give_memory_counts_one_and_two(self):
        result = g.native_controls()
        self.assertEqual(result['native_minimum_memory_dimensions'], [1, 2])
        self.assertEqual(result['native_flat_Gram_ranks'], [[2, 2], [3, 3]])

    def test_zero_variance_has_no_memory(self):
        data = tuple(F(2, 3)**k for k in range(7))
        result = g.three_jet(data[:4], 1)
        self.assertEqual(result['memory_dimension'], 0)
        self.assertEqual(g.realize(data, 1, 1)['memory_dimension'], 0)

    def test_weighted_two_rate_reconstruction_preserves_moments(self):
        data = tuple(F(2, 7)*F(1, 3)**k+F(5, 7)*F(4, 5)**k for k in range(11))
        result = g.realize(data, 2, 1)
        self.assertEqual(g.model_moments(result['L'], result['G'], 10), data)
        self.assertNotEqual(result['G'], g.a.identity(2))

    def test_three_derivatives_recover_GE3_memory_kernel(self):
        result = g.three_jet(pair_moments(3), 1)
        self.assertTrue(result['ceiling_exact'])
        self.assertEqual((result['q'], result['delta']), (F(1, 16), F(3, 4)))
        self.assertEqual(g.model_moments(result['L'], result['G'], 8), pair_moments())

    def test_two_derivatives_are_insufficient_even_with_the_ceiling(self):
        L = [[F(3, 4), F(1, 4)], [F(1, 4), F(1, 2)]]
        I = g.a.identity(2)
        self.assertTrue(g.a.psd(L) and g.a.psd(g.a.add(I, L, -1)))
        data = g.model_moments(L, I, 3)
        self.assertEqual(data[:3], pair_moments(2))
        self.assertNotEqual(data[3], pair_moments(3)[3])

    def test_ceiling_is_an_independent_source_requirement(self):
        result = g.boundary_controls()
        self.assertNotEqual(F(result['missing_ceiling_counterexample_fourth_moment']), F(17, 32))
        self.assertEqual(tuple(map(F, result['same_first_three_moments'])), pair_moments(3))

    def test_positive_localizer_does_not_prove_extra_memory(self):
        data = tuple((F(1, 4)**k+F(3, 4)**k)/2 for k in range(7))
        self.assertGreater(g.three_jet(data[:4], 1)['Xi'], 0)
        self.assertEqual(g.realize(data, 2, 1)['memory_dimension'], 1)

    def test_approximate_recovery_has_outward_error_enclosures(self):
        rows = g.native_controls()['two_state_approximation_enclosures']
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertLessEqual(F(row['error_upper'])**2, F(row['bound_squared']))

    def test_three_rate_probe_refuses_a_false_two_state_realization(self):
        data = tuple(sum((x**k for x in (F(1, 2), F(2, 3), F(5, 6))), F(0))/3 for k in range(9))
        with self.assertRaises(ValueError):
            g.realize(data, 2, 1)
        self.assertEqual(g.realize(data, 3, 1)['memory_dimension'], 2)

    def test_invisible_ambient_modes_cannot_be_counted_by_this_probe(self):
        self.assertEqual(g.boundary_controls()['invisible_appended_rates'], ['1/8', '1/100', '1/10000'])

    def test_tiny_nonzero_Gram_rank_is_not_rounded_to_zero(self):
        eps = F(1, 10**12)
        data = tuple((1-eps)*(F(1, 2)**k+1)/2+eps*F(3, 4)**k for k in range(9))
        self.assertEqual(g.psd_rank(g.hankel(data, 3)), 3)
        with self.assertRaises(ValueError):
            g.realize(data, 2, 1)
        self.assertEqual(g.realize(data, 3, 1)['memory_dimension'], 2)

    def test_input_flatness_and_source_interval_gates(self):
        self.assertEqual(len(g.refusals()), 9)

    def test_source_pin_mutation_fails(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source = root/'prior.md'; source.write_text('fixed\n')
            pins = {'upstream_sha256': {'prior.md': g.old.digest(source)}, 'local_inputs': []}
            g.source_checks(root, pins); source.write_text('changed\n')
            with self.assertRaises(ValueError):
                g.source_checks(root, pins)

    def test_observable_recovery_cannot_be_promoted_to_global_gap(self):
        cert = {'scope': 'minimal observable cyclic sector', 'verdict': 'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            result, expected = Path(folder)/'result.json', Path(folder)/'expected.sha256'
            result.write_text(json.dumps(cert)); expected.write_text(g.y.canonical_sha(cert))
            g.check(cert, result, expected)
            with self.assertRaises(ValueError):
                g.check(dict(cert, scope='global Yang-Mills mass gap'), result, expected)


if __name__ == '__main__':
    unittest.main()
