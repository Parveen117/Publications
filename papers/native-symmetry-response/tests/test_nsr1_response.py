"""Independent exact and failure controls for the electromagnetic adapter."""

from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

SOURCE = Path(__file__).resolve().parents[1]/'certificates/nsr1_response.py'
SPEC = importlib.util.spec_from_file_location('nsr1_test_adapter', SOURCE)
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)
ZERO = m.nt.mat([[0, 0], [0, 0]])


class ResponseTests(unittest.TestCase):
    def test_runtime_is_python_312(self):
        m.require_runtime()

    def test_taylor_square_root_replays_radicand(self):
        x, y = m.variables(F(1, 2), F(1, 2))
        radicand = (1+x*x)*(1+y*y)
        root = radicand.power(F(1, 2))
        self.assertEqual((root*root).c, radicand.c)
        self.assertEqual((root/root).c, {(0, 0): F(1)})
        self.assertEqual((x**2*y**2).derivative(0, 0, 1, 1), 4)

    def test_reciprocal_responses_solve_imposed_constraints(self):
        h = m.nt.mat([[3, 1], [1, 2]])
        r = m.responses(h)
        self.assertEqual(m.nt.mul(h, [[r['epsilon_h']], [r['mixed']]]), [[F(1)], [F(0)]])
        self.assertEqual(m.nt.mul(h, [[r['mixed']], [r['mu_e']]]), [[F(0)], [F(1)]])
        self.assertEqual(r['epsilon_h']/r['epsilon_y'], r['mu_e']/r['mu_x'])
        self.assertGreaterEqual(r['gain'], 1)

    def test_born_infeld_exact_half_bias(self):
        info = m.analyze(m.born_infeld(F(1, 2), F(1, 2)))
        self.assertEqual(info['H'], m.nt.mat([[F(4, 5), F(1, 5)], [F(1, 5), F(4, 5)]]))
        self.assertEqual(info['dH'][0], m.nt.mat([[F(-24, 25), F(8, 25)], [F(8, 25), F(8, 25)]]))
        self.assertEqual(info['delta'], F(3, 5))
        self.assertEqual(info['kappa'], F(128, 225))

    def test_point_matched_linear_model_misses_curvature(self):
        energy = m.born_infeld(F(1, 2), F(1, 2))
        linear = m.quadratic_tangent(energy)
        for axes in ((), (0,), (1,), (0, 0), (0, 1), (1, 1)):
            self.assertEqual(energy.derivative(*axes), linear.derivative(*axes))
        self.assertNotEqual(m.analyze(energy)['F'], ZERO)
        self.assertEqual(m.analyze(linear)['F'], ZERO)

    def test_zero_point_coupling_does_not_force_flatness(self):
        info = m.analyze(m.born_infeld(F(3, 4), 0))
        self.assertEqual(info['H'][0][1], 0)
        self.assertEqual(info['kappa'], F(9, 20))

    def test_separable_nonlinear_response_is_flat(self):
        x, y = m.variables(F(2, 3), F(1, 4))
        info = m.analyze(x*x/2+x**4/4+y*y/2+y**4/3)
        self.assertNotEqual(info['dH'][0], ZERO)
        self.assertEqual(info['F'], ZERO)

    def test_wave_reconstruction_and_signed_impedance(self):
        h = m.analyze(m.born_infeld(F(1, 2), F(1, 2)))['H']
        w = m.wave_data(h)
        self.assertEqual((w['v_plus'], w['v_minus'], w['Z']), (F(1), F(-3, 5), F(1)))
        self.assertEqual(m.reconstruct(w['v_plus'], w['v_minus'], w['Z']), h)
        for speed, sign in ((w['v_plus'], 1), (w['v_minus'], -1)):
            v = [[F(1)], [sign*w['Z']]]
            self.assertEqual(m.nt.mul([h[1], h[0]], v), m.nt.scale(v, speed))
            fields = m.nt.mul(h, v)
            self.assertEqual(fields[0][0]/fields[1][0], sign*w['Z'])

    def test_born_infeld_speed_bound_has_exact_square_witness(self):
        for x in (F(-3, 2), F(0), F(1, 2), F(2)):
            for y in (F(-2), F(0), F(1, 3), F(2)):
                q = (1+x*x)*(1+y*y)
                self.assertEqual(q-(x*y+1)**2, (x-y)**2)
                self.assertEqual(q-(x*y-1)**2, (x+y)**2)

    def test_common_variable_stiffness_changes_only_normalized_scalar(self):
        info = m.analyze(m.born_infeld(F(1, 2), F(1, 2)))
        w, dw = F(2), [F(3), F(-1)]
        h = m.nt.scale(info['H'], w)
        dh = [m.nt.add(m.nt.scale(info['dH'][i], w), m.nt.scale(info['H'], dw[i])) for i in range(2)]
        hi = m.nt.inv(h)
        xs = [m.nt.mul(hi, d) for d in dh]
        curvature = m.nt.scale(m.nt.sub(m.nt.mul(xs[0], xs[1]), m.nt.mul(xs[1], xs[0])), F(-1, 4))
        self.assertEqual(curvature, info['F'])
        self.assertEqual(m.shape_kappa(h, dh), info['kappa']/w)
        self.assertNotEqual(m.shape_kappa(h, dh), info['kappa'])

    def test_orientation_and_wrong_sign_control(self):
        h, dh, _ = m.metric_jets(m.born_infeld(F(1, 2), F(1, 2)))
        self.assertEqual(m.shape_kappa(h, list(reversed(dh))), -m.shape_kappa(h, dh))
        self.assertGreater(m.shape_kappa(h, dh), 0)

    def test_even_symmetry_hides_zero_bias_higher_response(self):
        u = m.born_infeld(0, 0)
        self.assertTrue(all(i % 2 == j % 2 == 0 for i, j in u.c))
        self.assertEqual(u.derivative(0, 0, 1, 1), 1)
        self.assertEqual(m.analyze(u)['F'], ZERO)
        self.assertNotEqual(m.analyze(m.born_infeld(F(1, 2), F(1, 2)))['F'], ZERO)

    def test_rotational_ward_identity_and_nonzero_curvature(self):
        self.assertEqual(m.radial_ward(F(2, 3), F(-3, 4), F(2)), ZERO)
        self.assertEqual(m.analyze(m.radial(F(1, 2), 0))['kappa'], F(-128, 1225))
        x, y = m.variables(F(2, 3), F(-3, 4))
        xr, yr = F(3, 5)*x-F(4, 5)*y, F(4, 5)*x+F(3, 5)*y
        self.assertEqual((xr*xr+yr*yr).c, (x*x+y*y).c)

    def test_full_vector_duality_uses_all_six_channels(self):
        d, b = list(map(F, (1, 2, 3))), list(map(F, (-2, 4, 1)))
        dr = [F(3, 5)*x-F(4, 5)*y for x, y in zip(d, b)]
        br = [F(4, 5)*x+F(3, 5)*y for x, y in zip(d, b)]
        self.assertEqual(m.bi_radicand(d, b), m.bi_radicand(dr, br))

    def test_schur_elimination_and_full_stability(self):
        x = y = F(1, 2)
        z = x*y
        hq = m.analyze(m.clamped(x, y, z))['H']
        effective = m.analyze(m.relaxed(x, y))
        self.assertEqual(m.schur_one(hq, [-y, -x], 1), effective['H'])
        full = m.nt.mat([[1, -z, -y], [-z, 1, -x], [-y, -x, 1]])
        self.assertGreater(full[0][0], 0)
        self.assertGreater(m.nt.det(hq), 0)
        self.assertGreater(m.nt.det(full), 0)
        self.assertEqual(effective['kappa'], F(128, 25))
        self.assertEqual(m.analyze(m.clamped(x, y, z))['kappa'], 0)

    def test_stability_and_domain_failures_are_not_silent(self):
        for point in ((1, 1), (2, 2)):
            with self.assertRaises(ValueError):
                m.analyze(m.born_infeld(*point))
        with self.assertRaises(ValueError):
            m.reconstruct(2, 1, 1)
        with self.assertRaises(ValueError):
            m.schur_one(m.nt.mat([[1, 0], [0, 1]]), [1, 1], 0)

    def test_arbitrary_fourth_order_changes_cancel_from_curvature(self):
        energy = m.born_infeld(F(1, 2), F(1, 2))
        changed = energy+m.Jet({(4, 0): 17, (3, 1): -2, (2, 2): 11, (0, 4): 3})
        first, second = m.analyze(energy), m.analyze(changed)
        self.assertNotEqual(first['ddH'], second['ddH'])
        self.assertEqual(first['F'], second['F'])

    def test_independent_existing_christoffel_fixtures(self):
        previous = m.nt.run()
        self.assertEqual(len(previous['fixtures']), 48)
        for fixture in previous['fixtures']:
            h = m.nt.mat(fixture['H'])
            dh = [m.nt.mat(d) for d in fixture['dH']]
            self.assertEqual(m.shape_kappa(h, dh), F(fixture['kappa']))


if __name__ == '__main__':
    unittest.main()
