"""NSB2: exact nonlinear response adapters and finite theorem controls.

Python 3.12 stdlib; no floating point, optimization, or inferred field claims.
Existing EMK arithmetic and NSB1/NSR1 geometry are imported unchanged.
Default/--check are read-only; --write updates only this stage's result/pin.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
from itertools import combinations, product
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import nsb1_geometry as geometry

SPEC = importlib.util.spec_from_file_location(
    'nsb2_emk', ROOT/'papers/emk-ugd-algebra/certificates/emk1_determinant_seam_ladder.py')
emk = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(emk)
nsr = geometry.nsr
RESULT = HERE/'NSB2_RESULT.json'
PIN = HERE/'EXPECTED_NSB2.sha256'
SOURCES = HERE/'NSB2_SOURCE_MANIFEST.json'
BASE_COMMIT = 'ff2b35e8c8749da00d212bdb820217e3faf8c5c2'
BASE_TREE = '558ce3f898cc287e1c287910fc9bf5bd47c3dc39'
MODEL = dict(
    EMK_basis_expansion_implies_linear_state_response=False,
    same_matrix_powers_repair_spectral_blindness=False,
    square_Hessian_rigidity_requires_an_open_set_identity=True,
    pointwise_square_Codazzi_test_certifies_neighborhood_flatness=False,
    flat_Hessian_implies_square_is_Hessian=False,
    potential_composition_preserves_Hessian_integrability=True,
    spectral_reconstruction_requires_simple_spectrum_and_calibration=True,
    two_dimensional_repair_uses_one_vector_probe_not_one_scalar=True,
    marked_cubic_deficit_is_binomial_n_3=True,
    cubic_data_count_claimed_minimal_for_curvature_alone=False,
    arbitrary_real_eigenframes_implemented_with_rational_arithmetic=False,
    pointwise_reconstruction_certifies_neighborhood_flatness=False,
    four_state_coordinates_are_four_dimensional_quantum_spacetime=False,
    YM82_plaquette_nonlinearity_is_Hessian_squaring=False,
    actual_Wilson_weak_coupling_window_extended=False,
    positive_four_dimensional_continuum_Yang_Mills_gap_proved=False,
    universal_analytic_proofs_formally_machine_verified=False)

exact = geometry.exact
mul, add, sub, scale = emk.mmul, emk.madd, emk.msub, emk.mscale
eye, inv, det, comm = emk.eye, emk.inv_exact, emk.det_exact, emk.commutator
R = emk.R2
K = scale(emk.RK2, -1)
S = mul(R, K)


def transpose(m):
    return [list(row) for row in zip(*m)]


def matrix(rows, n=None):
    if not isinstance(rows, (list, tuple)) or not rows:
        raise ValueError('A nonempty square matrix is required')
    n = len(rows) if n is None else n
    if len(rows) != n or any(not isinstance(row, (list, tuple)) or len(row) != n for row in rows):
        raise ValueError('Incompatible square matrix dimensions')
    return [[exact(v) for v in row] for row in rows]


def symmetric(rows):
    h = matrix(rows)
    if h != transpose(h):
        raise ValueError('A symmetric matrix is required')
    return h


def positive(rows):
    h = symmetric(rows)
    if any(det([row[:k] for row in h[:k]]) <= 0 for k in range(1, len(h)+1)):
        raise ValueError('A strictly positive metric is required')
    return h


def vector(values, n):
    if not isinstance(values, (list, tuple)) or len(values) != n:
        raise ValueError('Incompatible vector dimension')
    return [exact(v) for v in values]


def derivatives(h, dh):
    n = len(h)
    if not isinstance(dh, (list, tuple)) or len(dh) != n:
        raise ValueError('One matrix derivative per coordinate is required')
    out = [symmetric(m) for m in dh]
    if any(len(m) != n for m in out):
        raise ValueError('Incompatible derivative dimensions')
    return out


def codazzi_defect(dh):
    n = len(dh)
    return [[[dh[k][i][j]-dh[i][k][j] for j in range(n)]
             for i in range(n)] for k in range(n)]


def zero_tensor(t):
    return all(emk.is_zero(m) for m in t)


def hessian_jet(h, dh):
    h = positive(h)
    dh = derivatives(h, dh)
    if not zero_tensor(codazzi_defect(dh)):
        raise ValueError('Non-Hessian derivative data: Levi-Civita shortcut refused')
    return h, dh


def response_curvature(h, dh):
    """Pointwise LC curvature from a positive Hessian third jet, any dimension."""
    h, dh = hessian_jet(h, dh)
    x = [mul(inv(h), m) for m in dh]
    out = {f'{i},{j}': scale(comm(x[i], x[j]), F(-1, 4))
           for i, j in combinations(range(len(h)), 2)}
    return out


def sectional(h, curvature, i, j):
    return mul(h, curvature[f'{i},{j}'])[i][j]/(h[i][i]*h[j][j]-h[i][j]**2)


def emk_coefficients(h):
    h = matrix(h, 2)
    return [(h[0][0]+h[1][1])/2, (h[0][0]-h[1][1])/2,
            (h[1][0]-h[0][1])/2, (h[1][0]+h[0][1])/2]


def reduce_polynomial(h, coefficients):
    """Ascending coefficients; returns p(H)=r I+s H by Cayley-Hamilton."""
    h = matrix(h, 2)
    if not isinstance(coefficients, (list, tuple)) or not coefficients:
        raise ValueError('Supply at least the constant coefficient')
    r, s, t, d = F(0), F(0), emk.trace(h), det(h)
    for value in reversed(coefficients):
        r, s = exact(value)-s*d, r+s*t
    return dict(identity_coefficient=r, H_coefficient=s,
                matrix=add(scale(eye(2), r), scale(h, s)))


def polynomial_jet(h, dh, coefficients):
    """Full noncommutative product rule; does not assume H commutes with dH."""
    h, dh = hessian_jet(h, dh)
    if not isinstance(coefficients, (list, tuple)) or not coefficients:
        raise ValueError('Supply polynomial coefficients')
    coefficients = [exact(v) for v in coefficients]
    n = len(h)
    powers = [eye(n)]
    for _ in coefficients[1:]:
        powers.append(mul(powers[-1], h))
    out, dout = scale(h, 0), [scale(h, 0) for _ in range(n)]
    for k, coefficient in enumerate(coefficients):
        out = add(out, scale(powers[k], coefficient))
        for axis in range(n):
            for r in range(k):
                term = mul(mul(powers[r], dh[axis]), powers[k-1-r])
                dout[axis] = add(dout[axis], scale(term, coefficient))
    defect = codazzi_defect(dout)
    return dict(H=out, dH=dout, Codazzi_defect=defect,
                pointwise_Hessian_jet=zero_tensor(defect),
                neighborhood_Hessian_identity_certified=False,
                H_derivative_commutators=[comm(h, m) for m in dh])


def variable_spectral_residual(h, grad_u, grad_v):
    """For G=u H+v I, integrability requires grad v+adj(H) grad u=0."""
    h = positive(h)
    if len(h) != 2:
        raise ValueError('This scalar criterion is two-dimensional')
    gu, gv = vector(grad_u, 2), vector(grad_v, 2)
    adj = [[h[1][1], -h[0][1]], [-h[1][0], h[0][0]]]
    return [gv[i]+sum(adj[i][j]*gu[j] for j in range(2)) for i in range(2)]


def compose_potential(gradient, h, dh, first, second, third):
    h, dh = hessian_jet(h, dh)
    n = len(h)
    g = vector(gradient, n)
    f1, f2, f3 = map(exact, (first, second, third))
    out = [[f1*h[i][j]+f2*g[i]*g[j] for j in range(n)] for i in range(n)]
    dout = [[[f1*dh[k][i][j]+f2*(h[i][j]*g[k]+h[i][k]*g[j]+h[j][k]*g[i])
              +f3*g[i]*g[j]*g[k] for j in range(n)] for i in range(n)] for k in range(n)]
    out, dout = hessian_jet(out, dout)
    return dict(H=out, dH=dout, curvature=response_curvature(out, dout))


def frame_data(eigenvalues, gradients, frame=None):
    n = len(eigenvalues)
    if n < 2:
        raise ValueError('Use dimension at least two')
    values = vector(eigenvalues, n)
    if min(values) <= 0 or len(set(values)) != n:
        raise ValueError('Positive simple eigenvalues required')
    g = matrix(gradients, n)
    q = eye(n) if frame is None else matrix(frame, n)
    if mul(transpose(q), q) != eye(n):
        raise ValueError('A calibrated orthonormal frame is required')
    return values, g, q


def reconstruct_cubic(eigenvalues, gradients, mixed=None, frame=None):
    """G[i][k]=direction e_k derivative of eigenvalue i; Q columns are e_i.

    Finite adapter accepts rational frames. Every distinct triple must be
    explicitly supplied in n>=3, including a declared zero value.
    """
    values, g, q = frame_data(eigenvalues, gradients, frame)
    n = len(values)
    expected = set(combinations(range(n), 3))
    mixed = {} if mixed is None else mixed
    if not isinstance(mixed, dict) or set(mixed) != expected:
        raise ValueError('Supply exactly every missing marked mixed component')
    mixed = {key: exact(v) for key, v in mixed.items()}
    c = [[[F(0) for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for i, j, k in product(range(n), repeat=3):
        counts = Counter((i, j, k))
        if len(counts) == 3:
            value = mixed[tuple(sorted(counts))]
        else:
            repeated = max(counts, key=counts.get)
            other = next((axis for axis in counts if axis != repeated), repeated)
            value = g[repeated][other]
        c[k][i][j] = value
    diagonal = [[values[i] if i == j else F(0) for j in range(n)] for i in range(n)]
    h = mul(mul(q, diagonal), transpose(q))
    dh = [[[sum(q[a][i]*q[b][j]*q[k][ell]*c[ell][i][j]
                for i, j, ell in product(range(n), repeat=3))
            for b in range(n)] for a in range(n)] for k in range(n)]
    return dict(H=h, dH=dh, curvature=response_curvature(h, dh),
                missing_cubic_components=comb(n, 3),
                reconstructed_cubic_components=comb(n+2, 3))


def spectral_gradients(h, dh, frame):
    h, dh = hessian_jet(h, dh)
    n = len(h)
    q = matrix(frame, n)
    if mul(transpose(q), q) != eye(n):
        raise ValueError('Orthonormal frame required')
    diagonal = mul(mul(transpose(q), h), q)
    if any(diagonal[i][j] for i in range(n) for j in range(n) if i != j):
        raise ValueError('Frame must diagonalize H')
    values = [diagonal[i][i] for i in range(n)]
    if len(set(values)) != n:
        raise ValueError('Simple spectrum required')
    return [[sum(q[a][i]*q[b][i]*q[k][ell]*dh[k][a][b]
                 for a, b, k in product(range(n), repeat=3))
             for ell in range(n)] for i in range(n)]


def vector_probe(trace, v, response):
    v, response = vector(v, 2), vector(response, 2)
    norm2 = sum(x*x for x in v)
    if not norm2:
        raise ValueError('A nonzero probe is required')
    w = [-v[1], v[0]]
    a = sum(v[i]*response[i] for i in range(2))/norm2
    b = sum(w[i]*response[i] for i in range(2))/norm2
    small = [[a, b], [b, exact(trace)-a]]
    basis = [[v[0], w[0]], [v[1], w[1]]]
    return positive(scale(mul(mul(basis, small), transpose(basis)), 1/norm2))


def planar_curvature(eigenvalues, gradients):
    values, g, _ = frame_data(eigenvalues, gradients)
    if len(values) != 2:
        raise ValueError('Use two eigenvalues')
    a, c = values
    p, q, r, s = g[0][0], g[0][1], g[1][0], g[1][1]
    return -(a*(q*s-r*r)+c*(p*r-q*q))/(4*a*a*c*c)


def square_controls():
    out = []
    for label, x, y in [('scalar_point', 0, 0), ('nearby_failure', 1, 0)]:
        X, Y = nsr.variables(x, y)
        row = nsr.analyze((X*X+Y*Y+X*Y*Y)/2)
        squared = polynomial_jet(row['H'], row['dH'], [0, 0, 1])
        assert squared['pointwise_Hessian_jet'] == zero_tensor(squared['H_derivative_commutators'])
        out.append(dict(label=label, original_kappa=row['kappa'], square=squared))
    X, Y = nsr.variables(0, 0)
    Z = X+Y
    flat = nsr.analyze(X*X/2+Z*Z/2+Z**3/6)
    assert flat['kappa'] == 0
    failed = polynomial_jet(flat['H'], flat['dH'], [0, 0, 1])
    assert not failed['pointwise_Hessian_jet']
    out.append(dict(label='flat_does_not_imply_square_integrability', original=flat, square=failed))
    separated = nsr.analyze(X*X+Y*Y/2+X**3/6+Y**3/3)
    accepted = polynomial_jet(separated['H'], separated['dH'], [0, 1, 1])
    assert accepted['pointwise_Hessian_jet']
    assert all(emk.is_zero(m) for m in response_curvature(accepted['H'], accepted['dH']).values())
    out.append(dict(label='separated_positive_control', transformed=accepted))
    return out


def composition_control():
    X, Y = nsr.variables(1, 0)
    u = (X*X+Y*Y)/2
    base, direct = nsr.analyze(u), nsr.analyze(u+u*u)
    built = compose_potential([u.derivative(0), u.derivative(1)],
                              base['H'], base['dH'], 1+2*u.derivative(), 2, 0)
    assert built['H'] == direct['H'] and built['dH'] == direct['dH']
    assert built['curvature']['0,1'] == direct['F'] and direct['kappa'] == F(-1, 32)
    return dict(base=base, transformed=built, direct_kappa=direct['kappa'])


def marked_pair_control():
    rows = []
    for kind, q in [('flat', [[0, 1], [1, 0]]), ('curved', eye(2))]:
        row = geometry.halfplane_response(1, 2, 0, kind)
        g = spectral_gradients(row['H'], row['dH'], q)
        restored = reconstruct_cubic([1, F(1, 2)], g, frame=q)
        assert restored['H'] == row['H'] and restored['dH'] == row['dH']
        assert restored['curvature']['0,1'] == row['F']
        response = [sum(row['H'][i][j]*[1, 2][j] for j in range(2)) for i in range(2)]
        assert vector_probe(F(3, 2), [1, 2], response) == row['H']
        assert planar_curvature([1, F(1, 2)], g) == row['kappa']
        rows.append(dict(kind=kind, probe=[1, 2], response=response,
                         eigenframe=q, eigenvalue_gradients_in_frame=g,
                         kappa=row['kappa'], reconstructed=restored))
    return rows


def mixed_controls():
    out = []
    for n in (3, 4):
        for value in (0, 1, 2):
            mixed = {key: F(0) for key in combinations(range(n), 3)}
            mixed[0, 1, 2] = F(value)
            row = reconstruct_cubic(list(range(1, n+1)), scale(eye(n), 0), mixed)
            assert spectral_gradients(row['H'], row['dH'], eye(n)) == scale(eye(n), 0)
            curvature = sectional(row['H'], row['curvature'], 0, 1)
            assert curvature == F(value*value, 24)
            out.append(dict(dimension=n, tau=value, K12=curvature,
                            missing_count=comb(n, 3), curvature=row['curvature']))
    return out


def covariance_control():
    h = matrix([[2, 1], [1, 3]])
    change = matrix([[1, 1], [0, 1]])
    pull = lambda m: mul(mul(transpose(change), m), change)
    hy, background = pull(h), pull(eye(2))
    wrong, right = mul(hy, hy), mul(mul(hy, inv(background)), hy)
    assert wrong != pull(mul(h, h)) and right == pull(mul(h, h))
    return dict(H=h, coordinate_change=change, transformed_background=background,
                naive_square=wrong, covariant_square=right)


def controls():
    h = matrix([[3, 1], [1, 2]])
    h2 = matrix([[2, -1], [-1, 4]])
    coefficients = emk_coefficients(h)
    assert S == emk.K2 and mul(K, K) == eye(2) and mul(R, R) == scale(eye(2), -1)
    assert comm(K, R) == scale(S, -2)
    assert comm(h, h2) == scale(R, 2*(coefficients[3]*emk_coefficients(h2)[1]
                                                  -coefficients[1]*emk_coefficients(h2)[3]))
    return dict(native_basis=dict(K=K, R=R, RK=S, sample=coefficients),
                polynomial=reduce_polynomial(h, [1, -2, 3, -1, 2]),
                square_rigidity=square_controls(), composition=composition_control(),
                spectral_repair=marked_pair_control(), mixed=mixed_controls(),
                affine_covariance=covariance_control())


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def canonical_sha(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'NSB2' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('NSB2 source base mismatch')
    if not manifest['upstream_sha256'] or len(manifest['local_inputs']) != 4:
        raise ValueError('NSB2 manifest shape mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = sha256((ROOT/path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed upstream source: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return encode(dict(certificate_type='NSB2_NONLINEAR_RESPONSE_CLOSURE', verdict='PASS',
                       runtime_policy='Python 3.12 only', model=MODEL,
                       input_sha256=source_checks(), controls=controls(),
                       claim_status='LOCAL_RESPONSE_REPAIR__HESSIAN_POWER_RIGIDITY__CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'NSB2_NONLINEAR_RESPONSE_CLOSURE':
        raise ValueError('NSB2 scope mismatch')
    digest = canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('NSB2 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(digest+'\n')
    else:
        check(cert)
    print('NSB2 PASS', digest)
    print('Local response closure and potential rigidity; 4D continuum Yang-Mills remains open')


if __name__ == '__main__':
    main()
