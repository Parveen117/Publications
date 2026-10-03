import unittest
from fractions import Fraction as Q
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import interaction_model as im
import grading_model as gr
import source_response_model as sr


class SourceResponseTests(unittest.TestCase):
    def setUp(self):
        self.net = sr.network(4, [(0, 1, 1), (1, 2, 4), (2, 3, 9)], [3])

    def test_exact_source_solution_and_flux(self):
        field = sr.solve(self.net, [1, 0, 0])
        self.assertEqual(field, (Q(49, 36), Q(13, 36), Q(1, 9), 0))
        self.assertEqual(sr.currents(self.net, field), (1, 1, 1))
        self.assertEqual(sr.cut_flux(self.net, field, [0, 2]), 1)

    def test_energy_completion(self):
        field = sr.solve(self.net, [2, -1, 3])
        delta = (Q(2), Q(-3), Q(1), Q(0))
        shifted = tuple(a+b for a, b in zip(field, delta))
        self.assertEqual(sr.energy(self.net, shifted, [2, -1, 3])
                         - sr.energy(self.net, field, [2, -1, 3]),
                         sr.energy(self.net, delta, [0, 0, 0]))

    def test_unanchored_component_rejected(self):
        with self.assertRaises(ValueError):
            sr.network(4, [(0, 1, 1), (2, 3, 1)], [3])

    def test_nonpositive_stiffness_rejected(self):
        for c in [0, -1]:
            with self.assertRaises(ValueError):
                sr.network(2, [(0, 1, c)], [1])

    def test_symmetric_shells_reduce(self):
        net = sr.network(4, [(0, 1, 1), (0, 2, 1), (1, 3, 2), (2, 3, 2)], [3])
        field = sr.solve(net, [1, 0, 0])
        reduced = sr.radial([2, 4])['field']
        self.assertEqual(field, (reduced[0], reduced[1], reduced[1], 0))

    def test_flux_conservation_does_not_give_equal_edge_response(self):
        net = sr.network(4, [(0, 1, 1), (0, 2, 1), (1, 3, 1), (2, 3, 3)], [3])
        field = sr.solve(net, [1, 0, 0])
        self.assertEqual(sr.cut_flux(net, field, [0]), 1)
        self.assertNotEqual(sr.currents(net, field)[0], sr.currents(net, field)[1])

    def test_native_generator_is_licensed(self):
        g = sr.native_generator(self.net['h'])
        j = im.native_grading(3)
        self.assertEqual(gm.transpose(g), gm.scale(g, -1))
        self.assertEqual(gm.mul(gm.mul(j, g), j), gm.scale(g, -1))

    def test_generator_decomposes_into_native_brackets(self):
        h = self.net['h']
        zero = gm.scale(gm.identity(6), 0)
        result = zero
        for a in range(3):
            result = gm.add(result, gm.scale(gr.edge_generator(3, a, -1), h[a][a]))
        for a, b in [(0, 1), (1, 2)]:
            bracket = gm.commutator(im.mixing_generator(3, a, b), gr.edge_generator(3, a, -1))
            result = gm.add(result, gm.scale(bracket, h[a][b]))
        self.assertEqual(result, sr.native_generator(h))

    def test_affine_fixed_point(self):
        field = sr.solve(self.net, [1, 0, 0])
        state = tuple(v for x in field[:-1] for v in (x, 0))
        self.assertEqual(sr.affine_step(self.net['h'], [1, 0, 0, 0, 0, 0], state), state)

    def test_cost_preserving_flow_is_not_relaxation(self):
        h, b, state = self.net['h'], [1, 0, 0, 0, 0, 0], [0]*6
        after = sr.affine_step(h, b, state)
        self.assertNotEqual(after, tuple(state))
        self.assertEqual(sr.paired_cost(h, b, after), sr.paired_cost(h, b, state))
        self.assertNotEqual(sr.paired_cost(h, b, after), sr.stationary_cost(h, [1, 0, 0]))

    def test_cut_preserves_stationary_source_and_constant(self):
        h, b = self.net['h'], [1, 2, -1]
        cut = sr.eliminate(h, b, [0, 2])
        original = sr.vector(gm.mul(gm.inverse(h), sr.column(b)))
        reduced = sr.vector(gm.mul(gm.inverse(cut['h']), sr.column(cut['b'])))
        self.assertEqual(reduced, (original[0], original[2]))
        self.assertEqual(sr.stationary_cost(cut['h'], cut['b'])+cut['constant'], sr.stationary_cost(h, b))

    def test_naive_erasure_changes_answer(self):
        cut = sr.eliminate(self.net['h'], [1, 0, 0], [0])
        self.assertEqual(cut['h'], ((Q(36, 49),),))
        self.assertNotEqual(cut['h'][0][0], self.net['h'][0][0])

    def test_pair_cross_energy_has_no_probe_self_energy(self):
        h = self.net['h']
        values = [sr.source_only_cross_energy(h, [2, 0, 0], probe)
                  for probe in ([3, 0, 0], [0, 3, 0], [0, 0, 3])]
        self.assertEqual(values[1]-values[0], 6)
        self.assertEqual(values[2]-values[1], Q(3, 2))

    def test_inverse_square_and_source_linearity(self):
        f = sr.source_probe_force([n*n for n in range(1, 9)], 1)
        four = sr.source_probe_force([n*n for n in range(1, 9)], 4)
        self.assertEqual(f[3]/f[1], Q(1, 4))
        self.assertEqual(four[3]/f[3], 4)

    def test_other_growth_profiles_are_admitted(self):
        ratios = []
        for alpha in [0, 1, 2, 3]:
            f = sr.source_probe_force([n**alpha for n in range(1, 5)])
            ratios.append(f[3]/f[1])
        self.assertEqual(ratios, [Q(1), Q(1, 2), Q(1, 4), Q(1, 8)])

    def test_nonlinear_same_distance_different_source_control(self):
        for n in range(1, 9):
            c, d1, d4 = n**4, Q(1, n*n), Q(2, n*n)
            self.assertEqual(c*d1**2, 1)
            self.assertEqual(c*d4**2, 4)
            self.assertEqual(d4/d1, 2)

    def test_onsite_cost_breaks_conserved_flux(self):
        h = gm.add(self.net['h'], gm.identity(3))
        field = sr.vector(gm.mul(gm.inverse(h), sr.column([1, 0, 0])))+(Q(0),)
        currents = sr.currents(self.net, field)
        self.assertTrue(1 > currents[0] > currents[1] > currents[2] > 0)

    def test_radius_must_be_independent(self):
        # Relabelling n=k^2 makes the 1/n profile look exactly like 1/k^2.
        for k in range(1, 8):
            self.assertEqual(Q(1, k*k), Q(1, k)**2)

    def test_scale_gate_rejects_wrong_distance_law(self):
        gate = sr.scale_test(1, Q(1, 100), Q(1, 2), Q(1, 100), Q(1, 4))
        self.assertEqual(gate['verdict'], 'REJECTED')

    def test_scale_gate_keeps_compatible_data_unresolved(self):
        self.assertEqual(sr.scale_test(1, Q(1, 100), Q(1, 4), Q(1, 100), Q(1, 4))['verdict'], 'UNRESOLVED')
        self.assertEqual(sr.scale_test(1, 0, Q(3, 10), 0, Q(1, 4), Q(1, 10))['verdict'], 'UNRESOLVED')
        self.assertEqual(sr.scale_test(1, 0, Q(3, 10), 0, Q(1, 4), Q(1, 100))['verdict'], 'REJECTED')

    def test_scale_gate_rejects_nonlinear_source_response(self):
        self.assertEqual(sr.scale_test(1, Q(1, 100), 2, Q(1, 100), 4)['verdict'], 'REJECTED')

    def test_invalid_error_or_displacement_rejected(self):
        with self.assertRaises(ValueError):
            sr.scale_test(1, 1, 2, 0, 2)
        with self.assertRaises(ValueError):
            sr.source_probe_force([1], spacing=0)


if __name__ == '__main__':
    unittest.main()
