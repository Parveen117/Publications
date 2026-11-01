"""NSB2 adversarial checks, including independent full Christoffel curvature."""
from fractions import Fraction as Q
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path
import tempfile
import unittest

PATH = Path(__file__).resolve().parents[1]/'certificates/nsb2_nonlinear_response.py'
SPEC = importlib.util.spec_from_file_location('nsb2_under_test', PATH)
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def christoffel_cubic(h, dh):
    """Differentiate the full coordinate formula; fourth derivatives vanish.

    This route does not invoke the response commutator or spectral reconstruction.
    """
    n, hi = len(h), m.inv(h)
    gamma = [[[sum(hi[k][a]*(dh[i][a][j]+dh[j][a][i]-dh[a][i][j])
                   for a in range(n))/2 for j in range(n)]
              for i in range(n)] for k in range(n)]
    dhi = [[[ -sum(hi[a][p]*dh[ell][p][q]*hi[q][b]
                   for p in range(n) for q in range(n))
               for b in range(n)] for a in range(n)] for ell in range(n)]
    dg = [[[[sum(dhi[ell][k][a]*(dh[i][a][j]+dh[j][a][i]-dh[a][i][j])
                    for a in range(n))/2 for j in range(n)]
               for i in range(n)] for k in range(n)] for ell in range(n)]
    return {f'{i},{j}': [[dg[i][k][j][ell]-dg[j][k][i][ell]
                          +sum(gamma[k][i][a]*gamma[a][j][ell]
                               -gamma[k][j][a]*gamma[a][i][ell] for a in range(n))
                          for ell in range(n)] for k in range(n)]
            for i, j in combinations(range(n), 2)}


class NonlinearResponseTests(unittest.TestCase):
    def test_native_basis_adapter(self):
        self.assertEqual(m.K, [[1, 0], [0, -1]])
        self.assertEqual(m.S, m.emk.K2)
        for h in ([[3, 2], [2, 5]], [[1, -3], [4, 2]]):
            out = m.scale(m.eye(2), 0)
            for v, basis in zip(m.emk_coefficients(h), (m.eye(2), m.K, m.R, m.S)):
                out = m.add(out, m.scale(basis, v))
            self.assertEqual(out, h)

    def test_symmetric_powers_have_no_R_channel(self):
        h = m.matrix([[3, 1], [1, 2]])
        for coefficients in ([0, 0, 1], [2, -3, 4, 5], [1]*8):
            out = m.reduce_polynomial(h, coefficients)['matrix']
            self.assertEqual(m.emk_coefficients(out)[2], 0)
            self.assertTrue(m.emk.is_zero(m.comm(h, out)))

    def test_mixed_channel_creates_R(self):
        self.assertEqual(m.comm(m.K, m.S), m.scale(m.R, -2))

    def test_cayley_hamilton_against_direct_powers(self):
        for a, b, c in product((-1, 0, 2), repeat=3):
            h = m.matrix([[a, b], [b, c]])
            coefficients = [Q(1, 3), -2, 4, Q(-3, 2), 7]
            power, direct = m.eye(2), m.scale(h, 0)
            for coefficient in coefficients:
                direct = m.add(direct, m.scale(power, coefficient))
                power = m.mul(power, h)
            self.assertEqual(m.reduce_polynomial(h, coefficients)['matrix'], direct)

    def test_noncommuting_product_derivative_against_polynomial_jets(self):
        X, Y = m.nsr.variables(0, 0)
        u = X*X+3*Y*Y/2+X*X*Y/2+X*Y*Y/3
        h, dh, _ = m.nsr.metric_jets(u)
        # Direct symbolic Hessian entries, before multiplying matrices of jets.
        field = [[2+Y, X+2*Y/3], [X+2*Y/3, 3+2*X/3]]
        power = [[m.nsr.Jet(int(i == j)) for j in range(2)] for i in range(2)]
        direct = [[m.nsr.Jet(0) for _ in range(2)] for _ in range(2)]
        coefficients = [1, -2, 0, 3]
        for coefficient in coefficients:
            direct = [[direct[i][j]+coefficient*power[i][j] for j in range(2)] for i in range(2)]
            power = [[sum((power[i][k]*field[k][j] for k in range(2)), m.nsr.Jet(0))
                      for j in range(2)] for i in range(2)]
        actual = m.polynomial_jet(h, dh, coefficients)
        self.assertEqual(actual['H'], [[e.derivative() for e in row] for row in direct])
        self.assertEqual(actual['dH'], [[[e.derivative(k) for e in row] for row in direct] for k in range(2)])
        self.assertFalse(m.emk.is_zero(m.comm(h, dh[0])))

    def test_square_defect_equivalence_on_dense_cubic_jets(self):
        for n in (2, 3, 4):
            h = [[Q(i+1) if i == j else Q(0) for j in range(n)] for i in range(n)]
            dh = [[[Q((i+j+k) % 3-1) for j in range(n)] for i in range(n)] for k in range(n)]
            square = m.polynomial_jet(h, dh, [0, 0, 1])
            self.assertEqual(square['pointwise_Hessian_jet'], m.zero_tensor([m.comm(h, d) for d in dh]))
            self.assertFalse(square['pointwise_Hessian_jet'])

    def test_scalar_point_condition_does_not_imply_flatness(self):
        row = m.square_controls()[0]
        self.assertTrue(row['square']['pointwise_Hessian_jet'])
        self.assertFalse(row['square']['neighborhood_Hessian_identity_certified'])
        self.assertEqual(row['original_kappa'], Q(1, 4))
        self.assertFalse(m.square_controls()[1]['square']['pointwise_Hessian_jet'])

    def test_flatness_does_not_imply_square_integrability(self):
        row = m.square_controls()[2]
        self.assertEqual(row['original']['kappa'], 0)
        self.assertFalse(row['square']['pointwise_Hessian_jet'])
        defect = row['square']['Codazzi_defect']
        self.assertEqual([defect[1][0][0], defect[1][0][1]], [1, 1])

    def test_separated_positive_extension(self):
        row = m.square_controls()[3]['transformed']
        self.assertTrue(row['pointwise_Hessian_jet'])
        self.assertTrue(all(m.emk.is_zero(f) for f in m.response_curvature(row['H'], row['dH']).values()))

    def test_strictly_convex_cubic_spectral_function(self):
        row = m.geometry.halfplane_response(1, 2, 0)
        self.assertFalse(m.polynomial_jet(row['H'], row['dH'], [0, 0, 0, 1])['pointwise_Hessian_jet'])

    def test_affine_spectral_function_need_not_be_flat(self):
        row = m.geometry.halfplane_response(1, 2, 0)
        out = m.polynomial_jet(row['H'], row['dH'], [1, 2])
        self.assertTrue(out['pointwise_Hessian_jet'])
        self.assertFalse(m.emk.is_zero(m.response_curvature(out['H'], out['dH'])['0,1']))

    def test_square_does_not_automatically_define_LC_response(self):
        row = m.geometry.halfplane_response(1, 2, 0)
        squared = m.polynomial_jet(row['H'], row['dH'], [0, 0, 1])
        with self.assertRaisesRegex(ValueError, 'Non-Hessian'):
            m.response_curvature(squared['H'], squared['dH'])

    def test_negative_lambda_requires_positivity(self):
        h, dh = m.eye(2), [m.scale(m.eye(2), 0)]*2
        out = m.polynomial_jet(h, dh, [0, 1, -2])
        with self.assertRaisesRegex(ValueError, 'positive'):
            m.response_curvature(out['H'], out['dH'])
        accepted = m.polynomial_jet(h, dh, [0, 1, Q(-1, 2)])
        self.assertTrue(m.emk.is_zero(m.response_curvature(accepted['H'], accepted['dH'])['0,1']))

    def test_variable_coefficient_integrability_residual(self):
        row = m.geometry.halfplane_response(1, 2, 0)
        h, dh = row['H'], row['dH']
        gu = [m.emk.trace(d) for d in dh]
        gv = [-m.det(h)*m.emk.trace(m.mul(m.inv(h), d)) for d in dh]
        residual = m.variable_spectral_residual(h, gu, gv)
        defect = m.polynomial_jet(h, dh, [0, 0, 1])['Codazzi_defect']
        self.assertEqual([defect[1][0][0], defect[1][0][1]], [residual[1], -residual[0]])

    def test_affine_covariance_needs_background_metric(self):
        row = m.covariance_control()
        self.assertNotEqual(row['naive_square'], row['covariant_square'])

    def test_potential_composition_radial_control(self):
        row = m.composition_control()
        self.assertEqual(row['transformed']['H'], [[4, 0], [0, 2]])
        self.assertEqual(row['direct_kappa'], Q(-1, 32))

    def test_composition_third_derivative_chain_rule(self):
        X, Y = m.nsr.variables(Q(1, 3), Q(-1, 4))
        u = 2+X*X+2*Y*Y+X*Y/3+X**3/12
        base = m.nsr.analyze(u)
        direct = m.nsr.analyze(u+u*u+u**3)
        value = u.derivative()
        row = m.compose_potential([u.derivative(0), u.derivative(1)], base['H'], base['dH'],
                                  1+2*value+3*value*value, 2+6*value, 6)
        self.assertEqual(row['H'], direct['H'])
        self.assertEqual(row['dH'], direct['dH'])
        self.assertEqual(row['curvature']['0,1'], direct['F'])

    def test_force_outer_product_can_add_noncommuting_information(self):
        h = m.matrix([[2, 0], [0, 3]])
        row = m.compose_potential([1, 1], h, [m.scale(h, 0)]*2, 1, 1, 0)
        self.assertFalse(m.emk.is_zero(m.comm(h, row['H'])))

    def test_invalid_composition_metric_rejected(self):
        with self.assertRaises(ValueError):
            m.compose_potential([0, 0], m.eye(2), [m.scale(m.eye(2), 0)]*2, -1, 0, 0)

    def test_marked_probe_repairs_NSB1_entire_spectrum_pair(self):
        flat, curved = m.marked_pair_control()
        self.assertEqual(flat['response'], [Q(1, 2), 2])
        self.assertEqual(curved['response'], [1, 1])
        self.assertEqual([flat['kappa'], curved['kappa']], [0, Q(1, 16)])

    def test_vector_probe_arbitrary_nonzero_direction(self):
        h = m.matrix([[3, Q(-1, 2)], [Q(-1, 2), 2]])
        for v in ([1, 0], [0, -2], [Q(2, 3), Q(3, 7)]):
            response = [sum(h[i][j]*v[j] for j in range(2)) for i in range(2)]
            self.assertEqual(m.vector_probe(5, v, response), h)
        with self.assertRaises(ValueError):
            m.vector_probe(5, [0, 0], [0, 0])

    def test_one_scalar_probe_leaves_an_orientation_ambiguity(self):
        a, b = m.matrix([[3, 1], [1, 2]]), m.matrix([[3, -1], [-1, 2]])
        self.assertEqual((m.emk.trace(a), m.det(a), a[0][0]), (m.emk.trace(b), m.det(b), b[0][0]))
        self.assertNotEqual([row[0] for row in a], [row[0] for row in b])
        pair = [m.geometry.halfplane_response(1, 2, 0, kind) for kind in ('flat', 'curved')]
        self.assertEqual([sum(sum(row) for row in item['H']) for item in pair], [Q(3, 2)]*2)
        self.assertNotEqual(pair[0]['kappa'], pair[1]['kappa'])

    def test_rotated_spectral_reconstruction_against_potential(self):
        X, Y = m.nsr.variables(0, 0)
        x, y = (3*X+4*Y)/5, (-4*X+3*Y)/5
        u = x*x+3*y*y/2+x**3/6+x*x*y+3*x*y*y/2+2*y**3/3
        row = m.nsr.analyze(u)
        q = [[Q(3, 5), Q(-4, 5)], [Q(4, 5), Q(3, 5)]]
        g = m.spectral_gradients(row['H'], row['dH'], q)
        out = m.reconstruct_cubic([2, 3], g, frame=q)
        self.assertEqual(out['H'], row['H'])
        self.assertEqual(out['dH'], row['dH'])
        self.assertEqual(out['curvature']['0,1'], row['F'])
        self.assertEqual(m.planar_curvature([2, 3], g), row['kappa'])

    def test_eigenvector_sign_choices_leave_physical_curvature_unchanged(self):
        row = m.geometry.halfplane_response(1, 2, 0)
        for q in ([[1, 0], [0, 1]], [[-1, 0], [0, 1]], [[1, 0], [0, -1]]):
            g = m.spectral_gradients(row['H'], row['dH'], q)
            out = m.reconstruct_cubic([1, Q(1, 2)], g, frame=q)
            self.assertEqual(out['dH'], row['dH'])
            self.assertEqual(out['curvature']['0,1'], row['F'])

    def test_hidden_triples_have_real_curvature(self):
        rows = m.mixed_controls()
        self.assertEqual([r['K12'] for r in rows], [0, Q(1, 24), Q(1, 6)]*2)
        self.assertEqual([r['missing_count'] for r in rows], [1]*3+[4]*3)

    def test_full_3D_and_4D_curvature_against_Christoffel(self):
        for n in (3, 4):
            g = [[Q((i+2*j) % 5-2, 3) for j in range(n)] for i in range(n)]
            mixed = {key: Q(sum(key)+1, 5) for key in combinations(range(n), 3)}
            row = m.reconstruct_cubic(list(range(1, n+1)), g, mixed)
            self.assertEqual(m.spectral_gradients(row['H'], row['dH'], m.eye(n)), g)
            self.assertEqual(row['curvature'], christoffel_cubic(row['H'], row['dH']))
            for key, value in mixed.items():
                i, j, k = key
                self.assertEqual(row['dH'][k][i][j], value)

    def test_missing_mixed_components_cannot_be_silently_zero_filled(self):
        with self.assertRaisesRegex(ValueError, 'every missing'):
            m.reconstruct_cubic([1, 2, 3], m.scale(m.eye(3), 0))
        with self.assertRaises(ValueError):
            m.reconstruct_cubic([1, 2, 3], m.scale(m.eye(3), 0), {(0, 0, 1): 0})

    def test_repeated_eigenvalues_rejected(self):
        with self.assertRaisesRegex(ValueError, 'simple'):
            m.reconstruct_cubic([1, 1], m.eye(2))

    def test_nonorthogonal_or_wrong_eigenframe_rejected(self):
        with self.assertRaises(ValueError):
            m.reconstruct_cubic([1, 2], m.eye(2), frame=[[1, 1], [0, 1]])
        with self.assertRaises(ValueError):
            m.spectral_gradients([[2, 1], [1, 3]], [m.scale(m.eye(2), 0)]*2, m.eye(2))

    def test_inexact_and_malformed_data_rejected(self):
        for value in (0.1, True):
            with self.assertRaises(ValueError):
                m.reduce_polynomial(m.eye(2), [value, 1])
        for bad in ([], [[1, 2]], [[1, 0], [0, 1], [0, 0]]):
            with self.assertRaises(ValueError):
                m.matrix(bad)

    def test_evidence_replay_and_scope(self):
        fresh = m.run()
        self.assertEqual(m.check(fresh), m.PIN.read_text().strip())
        self.assertFalse(fresh['model']['positive_four_dimensional_continuum_Yang_Mills_gap_proved'])

    def test_modified_result_and_scope_cannot_pass(self):
        fresh = m.run()
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result.json', Path(directory)/'pin'
            changed = json.loads(json.dumps(fresh))
            changed['controls']['composition']['direct_kappa'] = '0'
            result.write_text(json.dumps(changed))
            pin.write_text(m.canonical_sha(changed))
            with self.assertRaises(ValueError):
                m.check(fresh, result, pin)
            changed['model']['positive_four_dimensional_continuum_Yang_Mills_gap_proved'] = True
            with self.assertRaisesRegex(ValueError, 'scope'):
                m.check(changed, result, pin)


if __name__ == '__main__':
    unittest.main()
