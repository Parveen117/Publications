"""Thermo response-plane geometry: exact finite rational implementation.

General statements are proved in THEOREM.md. A rational matrix fixture is not
a selected equation of state, a clock, a spacetime metric, or a particle model.
"""
from fractions import Fraction as Q
from math import isqrt


I = ((Q(1), Q(0)), (Q(0), Q(1)))
R = ((Q(0), Q(-1)), (Q(1), Q(0)))
ZERO = ((Q(0), Q(0)), (Q(0), Q(0)))


def matrix(rows):
    out = tuple(tuple(Q(v) for v in row) for row in rows)
    if len(out) != 2 or any(len(row) != 2 for row in out):
        raise ValueError('This response-plane implementation requires 2 by 2 matrices.')
    return out


def transpose(a):
    return tuple(zip(*a))


def add(a, b):
    return tuple(tuple(x+y for x, y in zip(r, s)) for r, s in zip(a, b))


def scale(a, q):
    return tuple(tuple(Q(q)*v for v in row) for row in a)


def multiply(a, b):
    return tuple(tuple(sum((x*y for x, y in zip(row, col)), Q(0))
                       for col in zip(*b)) for row in a)


def commutator(a, b):
    return add(multiply(a, b), scale(multiply(b, a), -1))


def trace(a):
    return a[0][0]+a[1][1]


def det(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def inverse(a):
    d = det(a)
    if not d:
        raise ValueError('Singular matrix.')
    return scale(((a[1][1], -a[0][1]), (-a[1][0], a[0][0])), 1/d)


def require_stable(h):
    if h != transpose(h) or h[0][0] <= 0 or det(h) <= 0:
        raise ValueError('A symmetric positive response Hessian is required.')


def rational_sqrt(q):
    q = Q(q)
    if q <= 0:
        raise ValueError('Positive value required.')
    n, d = isqrt(q.numerator), isqrt(q.denominator)
    if n*n != q.numerator or d*d != q.denominator:
        raise ValueError('Square root is not rational; use the general written theorem.')
    return Q(n, d)


def responses(h, temperature, volume):
    require_stable(h)
    t, v = Q(temperature), Q(volume)
    if t <= 0 or v <= 0:
        raise ValueError('Positive temperature and volume are required.')
    a, c, d = h[0][0], h[1][1], det(h)
    return dict(C_V=t/a, C_P=t*c/d, K_S=v*c, K_T=v*d/a)


def quarter_turn(h):
    require_stable(h)
    return scale(multiply(R, h), 1/rational_sqrt(det(h)))


def quarter_turn_derivative(h, dh):
    if dh != transpose(dh):
        raise ValueError('A symmetric derivative of the metric is required.')
    j = quarter_turn(h)
    return add(scale(multiply(R, dh), 1/rational_sqrt(det(h))),
               scale(j, -trace(multiply(inverse(h), dh))/2))


def connection(h, dh):
    """Affine Hessian-chart formula; torsion-freedom needs a symmetric cubic jet."""
    require_stable(h)
    if dh != transpose(dh):
        raise ValueError('A symmetric derivative of the metric is required.')
    return scale(multiply(inverse(h), dh), Q(1, 2))


def require_hessian_jet(hs, hv):
    if hs != transpose(hs) or hv != transpose(hv):
        raise ValueError('Each Hessian derivative must be symmetric.')
    if hs[0][1] != hv[0][0] or hs[1][1] != hv[0][1]:
        raise ValueError('Third derivatives must form a symmetric cubic tensor.')


def curvature(h, hs, hv):
    require_stable(h)
    require_hessian_jet(hs, hv)
    h_inv = inverse(h)
    return scale(commutator(multiply(h_inv, hs), multiply(h_inv, hv)), Q(-1, 4))


def gaussian_curvature(h, hs, hv):
    """Convention R_sv = [nabla_s,nabla_v]; K=g(R_sv e_v,e_s)/det(g)."""
    return multiply(h, curvature(h, hs, hv))[0][1]/det(h)


def christoffel_direct(h, dh):
    """Independent full Levi-Civita formula, including non-Hessian metrics.

    dh[i][j][k] denotes partial_i H_jk. Returns A_i^k_j.
    """
    require_stable(h)
    hi = inverse(h)
    return tuple(tuple(tuple(sum((hi[k][ell]*(dh[i][ell][j]
                    + dh[j][ell][i]-dh[ell][i][j]) for ell in range(2)), Q(0))/2
                    for j in range(2)) for k in range(2)) for i in range(2))


def curvature_direct(h, dh, ddh):
    """Curvature from differentiated full Christoffels, independent of -[X,Y]/4.

    ddh[r][i][j][k] denotes partial_r partial_i H_jk.
    """
    hi = inverse(h)
    a = christoffel_direct(h, dh)
    da = []
    for r in range(2):
        dhi = scale(multiply(multiply(hi, dh[r]), hi), -1)
        da.append(tuple(tuple(tuple(sum((
            dhi[k][ell]*(dh[i][ell][j]+dh[j][ell][i]-dh[ell][i][j])
            + hi[k][ell]*(ddh[r][i][ell][j]+ddh[r][j][ell][i]
                         - ddh[r][ell][i][j])
            for ell in range(2)), Q(0))/2
            for j in range(2)) for k in range(2)) for i in range(2)))
    return add(add(da[0][1], scale(da[1][0], -1)), commutator(a[0], a[1]))
