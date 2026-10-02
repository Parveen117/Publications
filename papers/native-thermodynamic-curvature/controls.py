#!/usr/bin/env python3
"""Independent rational coordinate controls, not a native-algebra implementation.

The Christoffel route differentiates the full formula; it does not start from
the commutator identity tested by the canonical native compiler.
"""
import json
from fractions import Fraction as F


def mat(rows):
    return [[F(x) for x in row] for row in rows]


def add(a, b):
    return [[x + y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def scale(a, s):
    return [[F(s) * x for x in row] for row in a]


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(len(b)))
             for j in range(len(b[0]))] for i in range(len(a))]


def tr(a):
    return [list(row) for row in zip(*a)]


def det(a):
    if len(a) == 2:
        return a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return sum((-1)**j * a[0][j] * det([r[:j] + r[j+1:] for r in a[1:]])
               for j in range(len(a)))


def inv(a):
    d = det(a)
    return scale([[a[1][1], -a[0][1]], [-a[1][0], a[0][0]]], 1/d)


def christoffel(h, dh, ddh):
    """Gamma_i and its coordinate derivative from the unsimplified formula."""
    hi = inv(h)
    dhi = [scale(mul(mul(hi, dh[p]), hi), -1) for p in range(2)]
    gamma = [[[sum(hi[k][l] * (dh[i][l][j] + dh[j][l][i] - dh[l][i][j])
                   for l in range(2))/2 for j in range(2)]
              for k in range(2)] for i in range(2)]
    dg = [[[[sum(
        dhi[p][k][l] * (dh[i][l][j] + dh[j][l][i] - dh[l][i][j])
        + hi[k][l] * (ddh[p][i][l][j] + ddh[p][j][l][i] - ddh[p][l][i][j])
        for l in range(2))/2 for j in range(2)] for k in range(2)]
        for i in range(2)] for p in range(2)]
    curvature = add(sub(dg[0][1], dg[1][0]),
                    sub(mul(gamma[0], gamma[1]), mul(gamma[1], gamma[0])))
    return gamma, dg, curvature


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, list):
        return [encode(y) for y in x]
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    return x


def run():
    checks = []
    fixtures = []
    z = mat([[0, 0], [0, 0]])
    for n in range(48):
        if n < 2:
            a, b, c = 2, 1, 5
            third = [0, 1 if n == 0 else 0, 0, 0]
        else:
            a, b, c = 2 + n % 5, n % 3 - 1, 3 + n % 7
            third = [(n * (k+2) + k*k) % 7 - 3 for k in range(4)]
        fourth = ([0]*5 if n < 2 else
                  [(n * (k+3) + 2*k) % 9 - 4 for k in range(5)])
        h = mat([[a, b], [b, c]])
        dh = [mat([[third[i], third[i+1]], [third[i+1], third[i+2]]])
              for i in range(2)]
        ddh = [[mat([[fourth[i+j], fourth[i+j+1]],
                     [fourth[i+j+1], fourth[i+j+2]]]) for j in range(2)]
               for i in range(2)]
        assert a > 0 and det(h) > 0
        gamma, dg, curvature = christoffel(h, dh, ddh)
        altered = [[add(ddh[i][j], mat([[i+j+1, i+j+2],
                                        [i+j+2, i+j+3]])) for j in range(2)]
                   for i in range(2)]
        _, _, changed = christoffel(h, dh, altered)
        assert curvature == changed, 'fourth-order cancellation'
        for i in range(2):
            assert add(mul(tr(gamma[i]), h), mul(h, gamma[i])) == dh[i]
            for j in range(2):
                for k in range(2):
                    assert gamma[i][k][j] == gamma[j][k][i], 'torsion'
        delta = det(h)
        kappa = mul(h, curvature)[0][1] / delta
        assert kappa == -det(mat([[a, b, c], third[:3], third[1:]]))/(4*delta**2)
        t, p, s, v = F(10), F(10), F(1), F(1)
        cv, cp, ks, kt = t/a, t*c/delta, v*c, v*delta/a
        response = [-s*delta/c, -s*a, -v*c, -v*delta/a, p/c, p*a/delta]
        assert cv/cp * ks/kt == 1
        assert response[0]*cp/t == response[1]*cv/t == -s
        assert response[2]*response[4] == response[5]*response[3] == -p*v
        fixtures.append(dict(id=n, H=h, dH=dh, ddH=ddh, Gamma=gamma,
                             dGamma=dg, curvature=curvature, kappa=kappa,
                             delta=delta, response=response))
    assert fixtures[0]['kappa'] == F(5, 324)
    assert fixtures[1]['curvature'] == z
    assert fixtures[0]['H'] == fixtures[1]['H']
    assert fixtures[0]['response'] == fixtures[1]['response']
    checks.append(dict(name='independent_full_Christoffel_differentiation', cases=48,
                       fourth_order_cancellation_cases=48, torsion_cases=48,
                       response_identity_cases=48))

    # Each is an affine coordinate pullback, with the tensor transported.
    for h in (fixtures[i]['H'] for i in range(12)):
        a = mat([[2, 1], [0, 3]])
        ai = inv(a)
        hp = mul(mul(tr(ai), h), ai)
        assert det(hp) == det(h)/det(a)**2
        x = mat([[2], [-1]])
        xp = mul(a, x)
        assert mul(mul(tr(xp), hp), xp) == mul(mul(tr(x), h), x)
    checks.append(dict(name='transported_response_metric', cases=12))

    # Variable reciprocal L=diag(1,1+x) has beta = x dx+(1+x)y dy.
    # Compute each oriented edge integral directly on [0,a] x [0,b].
    seam = []
    for a in (F(1, 2), F(1), F(3, 2)):
        for b in (F(1, 3), F(1), F(2)):
            edges = [a*a/2, (1+a)*b*b/2, -a*a/2, -b*b/2]
            flux = a*b*b/2  # integral y dx dy
            assert sum(edges) == flux and flux > 0
            # Piecewise eta=0 left, c dv right, positively oriented rectangle.
            c = F(3, 2)
            seam_edges = [F(0), c*b, F(0), F(0)]
            assert sum(seam_edges) == c*b
            assert sum(-x for x in seam_edges) == -c*b
            seam.append(dict(width=a, height=b, onsager_loop=flux, seam_loop=c*b))
    checks.append(dict(name='Onsager_derivative_and_Seam_Stokes_integrals',
                       variable_reciprocal_loops=9, oriented_seam_loops=18))

    b = mat([[2, 1], [0, 1]])
    x, j = mat([[1], [2]]), mat([[3], [-1]])
    xp, jp = mul(tr(inv(b)), x), mul(b, j)
    assert mul(tr(jp), xp) == mul(tr(j), x)
    assert mul(tr(mul(b, j)), x) != mul(tr(j), x)
    checks.append(dict(name='Onsager_dual_frame_covariance', cases=1,
                       nondual_change_rejected=True))

    negatives = dict(
        same_point_response_determines_curvature=False,
        reciprocal_variable_response_is_always_flat=False,
        scalar_ratio_closure_implies_native_flatness=False,
        flat_bulk_erases_seam=False,
        nondual_force_flux_change_preserves_production=False,
    )
    # Each false alternative has an explicit executed witness above.
    assert fixtures[0]['curvature'] != z and fixtures[0]['response'] == fixtures[1]['response']
    assert seam[0]['onsager_loop'] != 0 and seam[0]['seam_loop'] != 0
    return encode(dict(status='PASS_INDEPENDENT_RATIONAL_CONTROLS', fixtures=fixtures,
                       checks=checks, seam_integrals=seam,
                       rejected_false_alternatives=list(negatives)))


if __name__ == '__main__':
    print(json.dumps(run(), sort_keys=True, separators=(',', ':')))
