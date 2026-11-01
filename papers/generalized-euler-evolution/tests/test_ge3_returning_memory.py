"""Independent and adversarial finite controls for GE3's written proofs."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ge3_returning_memory as g


class ReturningMemoryTests(unittest.TestCase):
    def test_sector_is_from_native_polynomials_and_reference(self):
        result = g.native_sector_controls()
        self.assertEqual(result['native_eigenvalues'], ['1/2', '1'])
        self.assertEqual(result['Phi_pairing_in_visible_hidden_basis'], 'I/6')
        self.assertEqual(len(result['native_substitution_controls']), 4)

    def test_nonuniform_elimination_keeps_all_hidden_initial_sources(self):
        result = g.word_controls()
        self.assertEqual(result['nonuniform_recursion_checks'], 12)
        self.assertEqual(result['hidden_initial_sources'], 3)

    def test_nonstationary_return_need_not_be_positive(self):
        self.assertEqual(g.word_controls()['nonstationary_positive_step_negative_return'], '-1/100')

    def test_stationary_word_kernel_matches_direct_history(self):
        W = g.rational_step(F(1, 2))[0]
        a, b, d = W[0][0], W[0][1], W[1][1]
        visible = [F(1)]
        for n in range(6):
            visible.append(a*visible[n]+sum((b*b*d**(n-1-j)*visible[j] for j in range(n)), F(0)))
        for n, value in enumerate(visible):
            self.assertEqual(value, g.a.power(W, n)[0][0])
        self.assertGreater(visible[-1], a**6)

    def test_unequal_step_retained_and_reset_have_distinct_limits(self):
        result = g.refinement_controls()
        self.assertEqual(len(result['unequal_step_enclosures']), 4)
        self.assertGreater(F(result['different_limit_separation_lower']), 0)
        for row in result['unequal_step_enclosures']:
            self.assertLessEqual(F(row['retained_error_upper']), F(row['bound']))
            self.assertLessEqual(F(row['reset_error_upper']), F(row['bound']))

    def test_first_jet_agrees_while_second_jet_detects_memory(self):
        # F'(0)=-A is shared; F''(0)=A^2+B B^dagger differs from the reset law.
        second = g.a.multiply(g.L, g.L)[0][0]
        self.assertEqual(second-g.L[0][0]**2, F(1, 16))
        self.assertEqual(g.L[0][0], F(3, 4))
        self.assertEqual(second, F(5, 8))

    def test_same_visible_state_different_hidden_state_changes_the_next_step(self):
        W = g.rational_step(F(1, 3))[0]
        plus = g.a.apply(W, [F(1), F(1)])[0]
        minus = g.a.apply(W, [F(1), F(-1)])[0]
        self.assertGreater(plus, minus)
        self.assertEqual(plus-minus, 2*W[0][1])

    def test_continuum_kernel_resolvent_and_strict_defect(self):
        result = g.continuum_controls()
        self.assertEqual(result['closure_marker'], '1/16')
        self.assertEqual(len(result['resolvent_checks']), 3)
        self.assertEqual(len(result['strict_defect_and_tail_enclosures']), 4)

    def test_linear_degree_is_blind_to_this_quadratic_memory(self):
        C = g.a.diag((0, F(1, 2), F(1, 2)))
        self.assertEqual(g.q.operator(1, C), g.a.scale(g.a.identity(4), F(1, 4)))
        self.assertNotEqual(g.L[0][1], 0)

    def test_reset_rate_is_not_the_retained_asymptotic_rate(self):
        result = g.continuum_controls()
        self.assertEqual(result['initial_visible_rate'], '3/4')
        self.assertEqual(result['retained_asymptotic_rate'], '1/2')
        self.assertTrue(g.a.psd(g.a.add(g.L, g.a.scale(g.a.identity(2), F(1, 2)), -1)))
        self.assertFalse(g.a.psd(g.a.add(g.L, g.a.scale(g.a.identity(2), F(3, 4)), -1)))

    def test_type_history_and_tail_gates(self):
        self.assertEqual(len(g.refusals()), 6)

    def test_frozen_source_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); source = root/'prior.md'
            source.write_text('frozen\n')
            pins = {'upstream_sha256': {'prior.md': g.old.digest(source)}, 'local_inputs': []}
            g.source_checks(root, pins)
            source.write_text('edited\n')
            with self.assertRaises(ValueError):
                g.source_checks(root, pins)

    def test_finite_observer_cannot_be_promoted_to_the_interacting_row(self):
        cert = {'verdict': 'PASS', 'scope': 'declared native quadratic observer'}
        with tempfile.TemporaryDirectory() as folder:
            result, expected = Path(folder)/'result.json', Path(folder)/'expected.sha256'
            result.write_text(json.dumps(cert)); expected.write_text(g.y.canonical_sha(cert))
            g.check(cert, result, expected)
            with self.assertRaises(ValueError):
                g.check(dict(cert, scope='actual interacting row closure solved'), result, expected)


if __name__ == '__main__':
    unittest.main()
