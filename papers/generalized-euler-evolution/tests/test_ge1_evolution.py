"""Adversarial and independent finite witnesses for GE1, Python 3.12."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ge1_evolution as g


class GeneralizedEulerTests(unittest.TestCase):
    def test_complex_kernel_sign_and_resonant_line(self):
        self.assertEqual(g.change_residual((1,), (0,), (1,)), (0, 0))
        self.assertEqual(g.change_residual((-1,), (0,), (1,)), (2, 0))
        for k in range(-10, 11):
            self.assertEqual(g.change_residual((k, 1-k), (0, 0), (1, 1)), (0, 0))

    def test_real_change_operator_has_no_character_kernel(self):
        for k in range(-5, 6):
            for l in range(-5, 6):
                self.assertEqual(g.change_residual((k, l), (F(2, 3), F(-4, 7)), (0, 0))[0], 1)

    def test_bilateral_complex_flow_cannot_be_called_contractive(self):
        rows = g.domain_controls()['bilateral_complex_domain_controls']
        self.assertTrue(all(F(r['growth_lower']) > 1 for r in rows))
        self.assertTrue(all(F(r['decay_upper']) < 1 for r in rows))

    def test_one_fixed_smooth_vector_defeats_forward_joint_limit(self):
        for m in (4, 8, 16):
            n = m*m
            self.assertEqual(g.coefficient(n), F(1, 2**m))
            self.assertEqual(g.forward_growth(m), F(2)**(m*m-2*m))
        self.assertGreater(g.forward_growth(16), g.forward_growth(8))

    def test_all_power_domains_are_not_an_analytic_vector_certificate(self):
        rows = g.domain_controls()['smooth_nonanalytic_taylor_failure']
        self.assertTrue(all(F(r['single_coefficient_term']) >= F(r['lower_bound']) > 1 for r in rows))

    def test_cayley_matches_independent_rational_circle_formula(self):
        for alpha in (F(-7, 3), F(0), F(5, 2)):
            for h in (F(-1), F(1, 3), F(9)):
                x = alpha*h/2
                c, s = (1-x*x)/(1+x*x), 2*x/(1+x*x)
                D = [[F(0), -alpha], [alpha, F(0)]]
                self.assertEqual(g.cayley(D, h, g.a.identity(2)), [[c, -s], [s, c]])

    def test_cayley_actual_error_against_factorial_enclosures(self):
        result = g.phase_controls()
        self.assertEqual(result['exact_isometry_and_reversal_checks'], 32)
        self.assertEqual(len(result['independent_scalar_enclosures']), 24)

    def test_phase_step_handles_zero_and_large_frequencies(self):
        for alpha in (F(0), F(10**10)):
            D = [[F(0), -alpha], [alpha, F(0)]]
            I = g.a.identity(2)
            V = g.cayley(D, 1, I)
            self.assertEqual(g.a.multiply(g.a.transpose(V), V), I)

    def test_native_noncommuting_and_rotated_rank_two_controls(self):
        result = g.native_controls()
        counts = result['checks']
        self.assertEqual(counts['native_noncommuting_brackets'], 3)
        self.assertEqual(counts['native_positive_resolvents'], 12)
        self.assertEqual(counts['reference_pairing_contractions'], 12)
        self.assertNotEqual(F(result['rotated_rank_two_tensor'][0][1]), 0)

    def test_heat_factors_are_checked_on_native_modes(self):
        result = g.heat_controls()
        self.assertEqual(result['native_eigenmode_inverse_checks'], 32)
        for r in result['independent_heat_enclosures']:
            self.assertLessEqual(F(r['error_upper']), F(r['bound_per_unit_norm']))

    def test_heat_midpoint_stability_does_not_imply_positive_readout(self):
        L = [[F(2) if i == j else F(-1) for j in range(3)] for i in range(3)]
        I = g.a.identity(3)
        M = g.a.multiply(g.a.inverse(g.a.add(I, g.a.scale(L, 2))),
                         g.a.add(I, g.a.scale(L, -2)))
        self.assertEqual(M[0][0], F(-1, 7))
        B = g.resolvent(L, 4, I)
        self.assertTrue(all(x >= 0 for row in B for x in row))
        self.assertTrue(all(sum(row) == 1 for row in B))

    def test_zero_generator_resolvent_is_identity(self):
        I = g.a.identity(3)
        self.assertEqual(g.resolvent(g.zero(3), 123, I), I)

    def test_pell_exact_initial_pairs_and_strict_positive_rates(self):
        self.assertEqual([g.pell_pair(i) for i in range(1, 5)], [(1, 1), (3, 2), (7, 5), (17, 12)])
        for m in range(1, 25):
            u, v = g.pell_pair(m)
            self.assertEqual(abs(u*u-2*v*v), 1)
            gap = g.pell_gap_bracket(m)
            self.assertGreater(gap.lo, 0)
            self.assertGreater(gap.hi, gap.lo)
        self.assertLess(g.pell_gap_bracket(24).hi, F(1, 10**17))

    def test_near_resonance_is_not_promoted_to_a_uniform_floor(self):
        result = g.resonance_controls()
        for left, right in zip(result['near_resonant_modes'], result['near_resonant_modes'][1:]):
            self.assertLess(F(right['gap_upper']), F(left['gap_lower']))
        self.assertEqual(result['nonconstant_rational_stationary_mode'], [2, -1])

    def test_domain_and_type_gates_reject_incompatible_inputs(self):
        self.assertEqual(len(g.refusal_controls()), 14)

    def test_upstream_pin_tampering_fails(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path = root/'source.py'
            path.write_text('original\n')
            pins = {'upstream_sha256': {'source.py': g.digest(path)}, 'local_inputs': []}
            g.source_checks(root, pins)
            path.write_text('changed\n')
            with self.assertRaises(ValueError):
                g.source_checks(root, pins)

    def test_certificate_mutation_fails(self):
        cert = {'verdict': 'PASS', 'scope': 'finite controls'}
        with tempfile.TemporaryDirectory() as temp:
            result, expected = Path(temp)/'result.json', Path(temp)/'expected.sha256'
            result.write_text(json.dumps(cert))
            expected.write_text(g.y.canonical_sha(cert))
            g.check(cert, result, expected)
            with self.assertRaises(ValueError):
                g.check(dict(cert, scope='Clay proved'), result, expected)


if __name__ == '__main__':
    unittest.main()
