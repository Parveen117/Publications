"""Exact finite phase-field controls; hypotheses and general proofs in THEOREM.md.

Python 3.11/3.12, standard library. Time, cell complex and response weights
are declared model data; this module does not choose a physical vacuum.
"""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import importlib.util
import sys


ROOT = Path(__file__).resolve().parents[2]


def load_source(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


thermo = load_source('phase_upstream_thermo', 'papers/thermo-compass-foundations/model.py')
native = load_source('phase_upstream_scalar',
    'papers/uncut-cut-measurement/certificates/native_sources/native_summability.py')
C = native.NativeCutScalar


def transpose(a):
    return [list(row) for row in zip(*a)]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), 0)


def mv(a, x):
    return [dot(row, x) for row in a]


def mm(a, b):
    return [[dot(row, col) for col in zip(*b)] for row in a]


def ident(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def rank(a):
    a = [[Q(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        scale = a[r][c]
        a[r] = [x/scale for x in a[r]]
        for i in range(r+1, len(a)):
            factor = a[i][c]
            if factor:
                a[i] = [x-factor*y for x, y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def solve(a, b):
    n = len(a)
    if any(len(row) != n for row in a) or len(b) != n:
        raise ValueError('Square linear system required.')
    rows = [[Q(x) for x in row]+[Q(y)] for row, y in zip(a, b)]
    for c in range(n):
        pivot = next((i for i in range(c, n) if rows[i][c]), None)
        if pivot is None:
            raise ValueError('Singular linear system.')
        rows[c], rows[pivot] = rows[pivot], rows[c]
        divisor = rows[c][c]
        rows[c] = [x/divisor for x in rows[c]]
        for i in range(c+1, n):
            factor = rows[i][c]
            if factor:
                rows[i] = [x-factor*y for x, y in zip(rows[i], rows[c])]
    x = [Q(0)]*n
    for i in reversed(range(n)):
        x[i] = rows[i][-1]-dot(rows[i][i+1:n], x[i+1:])
    return x


def scalar_connection(h, dh):
    """Oriented upper-triangular orthonormal coframe; no sqrt(a) in result."""
    thermo.require_stable(h)
    a, b = h[0][0], h[0][1]
    root = thermo.rational_sqrt(thermo.det(h))
    return (dh[0][1]-b*dh[0][0]/a)/(2*root)


def scalar_curvature_direct(h, dh, ddh):
    """Differentiate the scalar connection, retaining the quartic tensor."""
    a, b, d = h[0][0], h[0][1], thermo.det(h)
    root = thermo.rational_sqrt(d)
    da = [v[0][0] for v in dh]
    db = [v[0][1] for v in dh]
    n = [db[i]-b*da[i]/a for i in range(2)]
    def derivative(r, i):
        dd = d*thermo.trace(thermo.multiply(thermo.inverse(h), dh[r]))
        dn = ddh[r][i][0][1]-(db[r]/a-b*da[r]/a**2)*da[i]-b*ddh[r][i][0][0]/a
        return dn/(2*root)-n[i]*dd/(4*d*root)
    return derivative(0, 1)-derivative(1, 0)


def coframe(h):
    a, b, d = h[0][0], h[0][1], thermo.det(h)
    ra, rd = thermo.rational_sqrt(a), thermo.rational_sqrt(d)
    return thermo.matrix([[ra, b/ra], [0, rd/ra]])


def coframe_derivative(h, dh):
    a, b, d = h[0][0], h[0][1], thermo.det(h)
    da, db = dh[0][0], dh[0][1]
    dd = d*thermo.trace(thermo.multiply(thermo.inverse(h), dh))
    ra, rd = thermo.rational_sqrt(a), thermo.rational_sqrt(d)
    return thermo.matrix([[da/(2*ra), db/ra-b*da/(2*a*ra)],
                          [0, dd/(2*rd*ra)-rd*da/(2*a*ra)]])


def wedge(u, v):
    return [[u[i]*v[j]-u[j]*v[i] for j in range(4)] for i in range(4)]


def wedge_square_coefficient(f):
    """Coefficient of dx0 wedge dx1 wedge dx2 wedge dx3 in F wedge F."""
    return 2*(f[0][1]*f[2][3]-f[0][2]*f[1][3]+f[0][3]*f[1][2])


def periodic_complex(n):
    if not isinstance(n, int) or n < 2:
        raise ValueError('Periodic cube side count must be an integer at least two.')
    vertices = list(product(range(n), repeat=3))
    lookup = {v: i for i, v in enumerate(vertices)}
    nv, ne = len(vertices), 3*len(vertices)
    def shift(v, axis):
        w = list(v)
        w[axis] = (w[axis]+1) % n
        return tuple(w)
    def edge(v, axis):
        return 3*lookup[v]+axis
    d0 = [[Q(0)]*nv for _ in range(ne)]
    d1 = [[Q(0)]*ne for _ in range(ne)]
    d2 = [[Q(0)]*ne for _ in range(nv)]
    for v in vertices:
        for axis in range(3):
            e = edge(v, axis)
            d0[e][lookup[v]] -= 1
            d0[e][lookup[shift(v, axis)]] += 1
            b, c = (axis+1) % 3, (axis+2) % 3
            for w, k, sign in ((shift(v, b), c, 1), (v, c, -1),
                               (shift(v, c), b, -1), (v, b, 1)):
                d1[e][edge(w, k)] += sign
            d2[lookup[v]][edge(shift(v, axis), axis)] += 1
            d2[lookup[v]][e] -= 1
    return d0, d1, d2


def phase_product(row, phases):
    result = C.one()
    for power, phase in zip(row, phases):
        if Q(power).denominator != 1:
            raise ValueError('Integral incidence required for compact phases.')
        factor = phase if power >= 0 else phase.dagger()
        for _ in range(abs(int(power))):
            result = result*factor
    return result


def positive_weights(eps, nu):
    if any(Q(x) <= 0 for x in list(eps)+list(nu)):
        raise ValueError('Strictly positive constitutive weights required.')


def field_generator(curl, eps, nu):
    positive_weights(eps, nu)
    ne, nf = len(eps), len(nu)
    if len(curl) != nf or any(len(row) != ne for row in curl):
        raise ValueError('Curl/weight dimensions disagree.')
    k = [[Q(0)]*(ne+nf) for _ in range(ne+nf)]
    for f in range(nf):
        for e in range(ne):
            k[e][ne+f] = curl[f][e]*nu[f]
            k[ne+f][e] = -curl[f][e]/eps[e]
    return k


def energy(state, eps, nu):
    positive_weights(eps, nu)
    ne = len(eps)
    return (sum(x*x/w for x, w in zip(state[:ne], eps))+
            sum(x*x*w for x, w in zip(state[ne:], nu)))/2


def midpoint_step(k, state, step):
    n, step = len(state), Q(step)
    left = [[Q(i == j)-step*k[i][j]/2 for j in range(n)] for i in range(n)]
    right = [x+step*y/2 for x, y in zip(state, mv(k, state))]
    return solve(left, right)


def fourier_curl(z):
    if len(z) != 3 or any(v.norm_square() != 1 for v in z):
        raise ValueError('Three unit native phases required.')
    d = [v-C.one() for v in z]
    zero = C.zero()
    return d, [[zero, -d[2], d[1]], [d[2], zero, -d[0]], [-d[1], d[0], zero]]


def native_adjoint(a):
    return [[x.dagger() for x in row] for row in transpose(a)]


def real_representation(a):
    n = len(a)
    out = [[Q(0)]*(2*n) for _ in range(2*n)]
    for i, row in enumerate(a):
        for j, z in enumerate(row):
            out[2*i][2*j], out[2*i][2*j+1] = z.radial, -z.turn
            out[2*i+1][2*j], out[2*i+1][2*j+1] = z.turn, z.radial
    return out
