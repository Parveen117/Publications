import copy
from fractions import Fraction as F
import importlib.util
import itertools
import json
from pathlib import Path
import tempfile
import unittest

MODULE = Path(__file__).resolve().parents[1]/'certificates/nm1_modular_observers.py'
spec = importlib.util.spec_from_file_location('nm1', MODULE)
n = importlib.util.module_from_spec(spec)
spec.loader.exec_module(n)


class NativeModularTests(unittest.TestCase):
    def test_01_complete_exact_controls(self):
        result = n.controls()
        self.assertEqual(result['finite_RK_group_order'], 8)
        self.assertGreater(result['integer_matrices_factored'], 100)
        self.assertEqual(result['cocycle_triples'], 1728)

    def test_02_source_basis_is_not_silently_replaced(self):
        self.assertEqual(n.K, n.mat(((0, 1), (1, 0))))
        self.assertNotEqual(n.K, n.adapted(n.K))
        self.assertEqual(n.native_chart(n.adapted(n.K)), n.K)
        self.assertEqual(n.T, n.mat(((F(1, 2), F(1, 2)), (-F(1, 2), F(3, 2)))))

    def test_03_lattice_is_load_bearing(self):
        with self.assertRaises(ValueError):
            n.sl2(n.T)
        self.assertEqual(n.sl2(n.adapted(n.T)), n.U)
        self.assertEqual(n.det(n.S), 2)
        for a, b in itertools.product(range(-3, 4), repeat=2):
            vector = (a-b, a+b)
            output = tuple(sum(n.T[i][j]*vector[j] for j in range(2)) for i in range(2))
            self.assertTrue(all(x.denominator == 1 for x in output))
            self.assertEqual((output[0]-output[1]) % 2, 0)

    def test_04_finite_rotation_does_not_bound_shear_history(self):
        self.assertEqual(n.power(n.R, 8), n.I)
        self.assertEqual(n.power(n.K, 8), n.I)
        self.assertEqual(len({n.power(n.T, k) for k in range(-30, 31)}), 61)

    def test_05_euclidean_recovery_handles_large_and_signed_inputs(self):
        for a, b, c in [(100003, -99991, 77), (-143, 233, -377), (0, -1, 1)]:
            A = n.mm(n.mm(n.power(n.U, a), n.R), n.mm(n.power(n.U, b), n.mm(n.R, n.power(n.U, c))))
            self.assertEqual(n.eval_word(n.factor(A)), A)
            self.assertEqual(n.eval_word(n.factor(n.scale(A, -1))), n.scale(A, -1))

    def test_06_naive_signless_lift_is_wrong(self):
        g = n.section(n.R)
        self.assertEqual(n.star(g, g), n.I)
        self.assertEqual(n.cocycle(g, g), -1)
        self.assertNotEqual(n.mm(g, g), n.star(g, g))
        self.assertEqual(n.unlift(n.compose((1, g), (1, g))), n.scale(n.I, -1))

    def test_07_every_representative_retains_the_order_two_obstruction(self):
        g = n.section(n.R)
        for choice in (-1, 1):
            self.assertEqual(n.power(n.scale(g, choice), 2), n.scale(n.I, -1))
        self.assertEqual(n.cocycle(n.I, g), 1)
        self.assertEqual(n.cocycle(g, n.I), 1)

    def test_08_section_change_obeys_coboundary_law(self):
        # A nonconstant normalized sign choice, not a homomorphic section.
        b = lambda A: -1 if n.section(A)[0][1] > 0 else 1
        self.assertEqual(b(n.I), 1)
        for A, B in itertools.product(n.sample_group()[:10], repeat=2):
            gh = n.star(A, B)
            left = n.mm(n.scale(n.section(A), b(A)), n.scale(n.section(B), b(B)))
            sigma = b(A)*b(B)*b(gh)*n.cocycle(A, B)
            self.assertEqual(left, n.scale(n.section(gh), sigma*b(gh)))

    def test_09_projective_and_conjugation_blindness_persist(self):
        A = n.mm(n.R, n.U)
        for B, C in itertools.product(n.sample_group()[:6], repeat=2):
            self.assertEqual(n.section(n.mm(n.mm(B, A), C)), n.section(n.mm(n.mm(B, n.scale(A, -1)), C)))
        for X in (n.I, n.K, n.R, n.N):
            self.assertEqual(n.mm(n.mm(A, X), n.inv(A)), n.mm(n.mm(n.scale(A, -1), X), n.inv(n.scale(A, -1))))

    def test_10_scalar_blindness_can_return(self):
        tr = lambda A: A[0][0]+A[1][1]
        A, B, continuation = n.U, n.inv(n.U), n.mm(n.R, n.U)
        self.assertEqual(tr(A), tr(B))
        self.assertNotEqual(n.section(A), n.section(B))
        self.assertEqual((tr(n.mm(continuation, A))**2, tr(n.mm(continuation, B))**2), (4, 0))

    def test_11_one_sign_is_not_a_full_history(self):
        self.assertEqual(n.lift(n.power(n.R, 0)), n.lift(n.power(n.R, 4)))
        self.assertNotEqual(n.lift(n.power(n.R, 0)), n.lift(n.power(n.R, 2)))
        self.assertEqual(n.section(n.power(n.R, 0)), n.section(n.power(n.R, 2)))

    def test_12_even_odd_weight_obstruction(self):
        z, g = n.zpair(F(2, 3), F(4, 5)), n.section(n.R)
        left = n.zm(n.j(g, n.act(g, z)), n.j(g, z))
        self.assertEqual(left, n.zpair(-1))
        self.assertEqual(n.zp(left, 2), n.zpair(1))
        self.assertNotEqual(left, n.j(n.star(g, g), z))

    def test_13_completion_equation_is_not_unique(self):
        z = n.zpair(1, 2); f = lambda w: w
        for A in n.sample_group():
            for constant in (0, 1):
                C = lambda w: n.zs(n.zpair(constant), w)
                self.assertEqual(n.za(n.defect(f, A, 0, z), n.defect(C, A, 0, z)), n.zpair(0))

    def test_14_refuses_wrong_domains_and_approximate_arithmetic(self):
        bad = [lambda: n.sl2(n.K), lambda: n.sl2(((2, 0), (0, 1))),
               lambda: n.sl2(((1, 0.0), (0, 1))), lambda: n.power(n.R, True),
               lambda: n.factor(((1, 0), (0, F(1, 2)))), lambda: n.act(n.R, n.zpair(0, 0)),
               lambda: n.zpair(0.1, 1), lambda: n.zp(n.zpair(1, 1), F(1, 2)),
               lambda: n.unlift((0, n.I)), lambda: n.unlift((True, n.I)),
               lambda: n.unlift((1, n.scale(n.I, -1))), lambda: n.eval_word([('K', 1)])]
        for fn in bad:
            with self.assertRaises((ValueError, TypeError)):
                fn()

    def test_15_source_tamper_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); p = root/'input'; p.write_text('native')
            pins = {'upstream_sha256': {'input': n.digest(p)}, 'local_inputs': []}
            n.source_checks(root, pins)
            p.write_text('changed')
            with self.assertRaises(ValueError):
                n.source_checks(root, pins)

    def test_16_evidence_and_claim_tamper_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            result, expected = Path(tmp)/'result.json', Path(tmp)/'digest'
            cert = {'scope': copy.deepcopy(n.SCOPE), 'value': 1}
            result.write_text(json.dumps(cert)); expected.write_text(n.canonical(cert))
            n.check(cert, result, expected)
            cert['value'] = 2
            with self.assertRaises(ValueError):
                n.check(cert, result, expected)
            cert['scope']['open'] = []
            result.write_text(json.dumps(cert)); expected.write_text(n.canonical(cert))
            with self.assertRaises(ValueError):
                n.check(cert, result, expected)


if __name__ == '__main__':
    unittest.main()
