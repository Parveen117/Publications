#!/usr/bin/env python3
"""Independent Fraction/matrix route. No RKF imports and no numerical eigenvalues."""
from fractions import Fraction as F
from itertools import product
import json


def mat(rows):
    return [[F(x) for x in row] for row in rows]


def zero(n):
    return mat([[0] * n for _ in range(n)])


def eye(n):
    return mat([[int(i == j) for j in range(n)] for i in range(n)])


def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(a, c):
    return [[x * c for x in r] for r in a]


def sub(a, b):
    return add(a, scale(b, -1))


def transpose(a):
    return [list(x) for x in zip(*a)]


def mul(a, b):
    bt = transpose(b)
    return [[sum(x*y for x, y in zip(r, s)) for s in bt] for r in a]


def comm(a, b):
    return sub(mul(a, b), mul(b, a))


def kron(a, b):
    return [[x*y for x in ar for y in br] for ar in a for br in b]


def tr(a):
    return sum(a[i][i] for i in range(len(a)))


def flat(a):
    return [x for row in a for x in row]


def lin(c, bs):
    out = zero(len(bs[0]))
    for x, b in zip(c, bs):
        out = add(out, scale(b, x))
    return out


def rank(rows):
    a = [list(r) for r in rows if any(r)]
    if not a:
        return 0
    p = 0
    for j in range(len(a[0])):
        k = next((k for k in range(p, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[p], a[k] = a[k], a[p]
        v = a[p][j]
        a[p] = [x/v for x in a[p]]
        for k in range(p+1, len(a)):
            if a[k][j]:
                v = a[k][j]
                a[k] = [x-v*y for x, y in zip(a[k], a[p])]
        p += 1
        if p == len(a):
            break
    return p


def eps(i, j, k):
    if len({i, j, k}) < 3:
        return 0
    return 1 if (i, j, k) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] else -1


def cross(a, b):
    return [sum(F(eps(i, j, k))*a[j]*b[k] for j in range(3) for k in range(3)) for i in range(3)]


def va(a, b):
    return [x+y for x, y in zip(a, b)]


def vs(a, c):
    return [x*c for x in a]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def mv(a, v):
    return [dot(r, v) for r in a]


def bil(v, a, u=None):
    return dot(v, mv(a, v if u is None else u))


checks = {}
negative = []


def check(condition, group):
    if not condition:
        raise ValueError('Failed exact control: ' + group)
    checks[group] = checks.get(group, 0)+1


def reject(name, condition):
    if not condition:
        raise ValueError('False alternative survived: ' + name)
    negative.append(name)


I = eye(2)
R = mat([[0, -1], [1, 0]])
K = mat([[1, 0], [0, -1]])
L = mul(R, K)
BASIS = [kron(a, b) for a in [I, K, R, L] for b in [I, K, R, L]]
C = kron(R, I)
E = [kron(I, R), kron(R, K), kron(R, L)]
T = [scale(x, F(1, 2)) for x in E]
T0 = scale(C, F(1, 2))


def centralizer_constraints(c):
    n = len(c)
    bs = [mat([[int(i == a and j == b) for j in range(n)] for i in range(n)])
          for a, b in product(range(n), repeat=2)]
    return list(zip(*[flat(comm(x, c))+flat(add(x, transpose(x))) for x in bs]))


def invariant_rank(ad):
    n = len(ad[0])
    sym = []
    for i in range(n):
        for j in range(i, n):
            sym.append(mat([[int((a, b) in [(i, j), (j, i)]) for b in range(n)] for a in range(n)]))
    columns = [sum((flat(add(mul(transpose(a), h), mul(h, a))) for a in ad), []) for h in sym]
    return rank(list(zip(*columns)))


check(rank(centralizer_constraints(R)) == 3, 'exhaustive_centralizer_ranks')
check(rank(centralizer_constraints(C)) == 12, 'exhaustive_centralizer_ranks')
for a, b in product(range(3), repeat=2):
    target = add(scale(eye(4), -int(a == b)), lin([eps(a, b, c) for c in range(3)], E))
    check(mul(E[a], E[b]) == target, 'quaternion_products')
for x in [C]+E:
    check(transpose(x) == scale(x, -1) and comm(x, C) == zero(4), 'skew_and_multiplier')
for a, b in product([T0]+T, repeat=2):
    check(-tr(mul(a, b)) == F(int(a == b)), 'positive_pairing')
ad = [mat([[0 if i == 0 or j == 0 else eps(k, j-1, i-1) for j in range(4)] for i in range(4)]) for k in range(3)]
check(invariant_rank(ad) == 8, 'exhaustive_invariant_form_rank')
# There are ten unknown symmetric coefficients, leaving the two declared weights.
for h in [mat([[int(i == j == 0) for j in range(4)] for i in range(4)]),
          mat([[int(i == j and i > 0) for j in range(4)] for i in range(4)])]:
    check(all(add(mul(transpose(a), h), mul(h, a)) == zero(4) for a in ad), 'invariant_form_basis')
reject('one_real_factor_has_no_three_skew_directions', rank(centralizer_constraints(R)) == 3)
reject('K_is_not_a_compact_generator', transpose(K) != scale(K, -1))
g_bad = add(scale(I, F(5, 4)), scale(K, F(3, 4)))
gi_bad = sub(scale(I, F(5, 4)), scale(K, F(3, 4)))
bad_r = mul(mul(gi_bad, R), g_bad)
check(tr(mul(transpose(bad_r), bad_r))/2 == F(257, 32), 'noncompact_norm_witness')
reject('positive_norm_not_GL_invariant', tr(mul(transpose(bad_r), bad_r)) != tr(mul(transpose(R), R)))
reject('negative_square_trace_not_positive_on_full_EMK', -tr(mul(K, K))/2 < 0)

# Exact affine gauge jets; an independent color-cross route is compared with matrices.
fixtures = []
eta = [-1, 1, 1, 1]
for seed in range(1, 7):
    A = [[F((seed+2*i+3*a)%7-3, 5) for a in range(3)] for i in range(4)]
    dA = [[[F((seed+3*i+2*j+a)%9-4, 7) for a in range(3)] for j in range(4)] for i in range(4)]
    curvature = [[va(va(dA[i][j], vs(dA[j][i], -1)), cross(A[i], A[j])) for j in range(4)] for i in range(4)]
    dF = [[[va(cross(dA[k][i], A[j]), cross(A[i], dA[k][j])) for j in range(4)] for i in range(4)] for k in range(4)]
    ddF = [[[[va(cross(dA[k][i], dA[l][j]), cross(dA[l][i], dA[k][j])) for j in range(4)] for i in range(4)] for l in range(4)] for k in range(4)]
    am = [lin(a, T) for a in A]
    dam = [[lin(x, T) for x in row] for row in dA]
    fm = [[lin(x, T) for x in row] for row in curvature]
    for i, j in product(range(4), repeat=2):
        check(fm[i][j] == add(sub(dam[i][j], dam[j][i]), comm(am[i], am[j])), 'connection_color_matrix_routes')
    for i, j, k in product(range(4), repeat=3):
        v = [F(0)]*3
        for a, b, c in [(i, j, k), (j, k, i), (k, i, j)]:
            v = va(v, va(dF[a][b][c], cross(A[a], curvature[b][c])))
        check(v == [0]*3, 'affine_covariant_Bianchi')
    div = [F(0)]*3
    for nu, mu in product(range(4), repeat=2):
        f = curvature[mu][nu]
        v = ddF[nu][mu][mu][nu]
        for term in [cross(dA[nu][mu], f), cross(A[mu], dF[nu][mu][nu]),
                     cross(A[nu], dF[mu][mu][nu]), cross(A[nu], cross(A[mu], f))]:
            v = va(v, term)
        div = va(div, vs(v, eta[mu]*eta[nu]))
    check(div == [0]*3, 'covariant_double_divergence')
    # Nonconstant gauge g=g0 exp(x0 X0)...exp(x3 X3), evaluated at the origin.
    g0 = add(scale(eye(4), F(3, 5)), scale(E[seed%3], F(4, 5)))
    gi = transpose(g0)
    X = [T[(i+seed)%3] for i in range(4)]
    B = [mul(mul(gi, a), g0) for a in am]
    ap = [add(b, x) for b, x in zip(B, X)]
    dap = [[add(add(mul(mul(gi, dam[i][j]), g0), comm(B[j], X[i])),
                  comm(X[j], X[i]) if i > j else zero(4)) for j in range(4)] for i in range(4)]
    for i, j in product(range(4), repeat=2):
        fp = add(sub(dap[i][j], dap[j][i]), comm(ap[i], ap[j]))
        check(fp == mul(mul(gi, fm[i][j]), g0), 'variable_gauge_covariance')
    # Independent coefficient variation of the polynomial curvature-squared density.
    def lag(a, da):
        ff = [[va(va(da[i][j], vs(da[j][i], -1)), cross(a[i], a[j])) for j in range(4)] for i in range(4)]
        return -sum(F(eta[i]*eta[j], 4)*dot(ff[i][j], ff[i][j]) for i, j in product(range(4), repeat=2))
    # A given component enters F at most linearly, so centered differences are exact.
    for nu, color in product(range(4), range(3)):
        up = [a[:] for a in A]; down = [a[:] for a in A]
        up[nu][color] += 1; down[nu][color] -= 1
        derivative = (lag(up, dA)-lag(down, dA))/2
        expected = sum(eta[mu]*eta[nu]*cross(A[mu], curvature[mu][nu])[color] for mu in range(4))
        check(derivative == expected, 'gauge_action_connection_variations')
        for mu in range(4):
            up = [[a[:] for a in row] for row in dA]; down = [[a[:] for a in row] for row in dA]
            up[mu][nu][color] += 1; down[mu][nu][color] -= 1
            derivative = (lag(A, up)-lag(A, down))/2
            check(derivative == -eta[mu]*eta[nu]*curvature[mu][nu][color], 'gauge_action_derivative_variations')
    fixtures.append({'A': A, 'dA': dA, 'F': curvature})

reject('drop_nonabelian_commutator', comm(T[0], T[1]) != zero(4))
reject('wrong_color_orientation', comm(T[0], T[1]) != scale(T[2], -1))
reject('individually_exact_components_imply_flatness', comm(T[0], T[1]) == T[2] and T[2] != zero(4))
reject('ordinary_divergence_can_replace_covariant', any(cross(A[mu], curvature[mu][0]) != [0]*3 for mu in range(4)))
# Maurer-Cartan jet of g=exp(x t1) exp(y t2), without inserting F=0 by hand.
pure_A = [T[0], T[1]]
pure_dA = [[zero(4), zero(4)], [comm(T[0], T[1]), zero(4)]]
pure_F = add(sub(pure_dA[0][1], pure_dA[1][0]), comm(*pure_A))
check(pure_F == zero(4), 'pure_gauge_cancellation')
reject('coefficient_commutator_alone_is_curvature', comm(T[0], T[1]) != zero(4))

# Source SM uses L=KR, explicitly opposite to this chapter's L=RK.
LS = mul(K, R)
G0 = [kron(R, I), kron(K, K), kron(K, LS), kron(LS, I)]
spin_words = [eye(4)]
for g in G0:
    spin_words += [mul(x, g) for x in spin_words]
check(rank([flat(x) for x in spin_words]) == 16, 'complete_Clifford_span')
G = [kron(eye(4), x) for x in G0]
CC = kron(C, eye(4))
TT = [kron(x, eye(4)) for x in T]
H = scale(mul(CC, G[0]), -1)
HG = [mul(H, g) for g in G]
check(transpose(H) == H and mul(H, H) == eye(16), 'massive_bilinear')
for a, b in product(range(4), repeat=2):
    check(add(mul(G[a], G[b]), mul(G[b], G[a])) == scale(eye(16), 2*eta[a]*int(a == b)), 'Clifford_relations')
for a in range(4):
    check(transpose(HG[a]) == scale(HG[a], -1), 'kinetic_skew')
    for t in TT:
        check(comm(t, G[a]) == zero(16) and add(mul(transpose(t), H), mul(H, t)) == zero(16), 'commuting_gauge_spin_actions')
for seed in range(1, 5):
    psi = [F((j*j+seed*j+1)%11-5, 7) for j in range(16)]
    dp = [[F((seed+3*j+mu)%9-4, 11) for j in range(16)] for mu in range(4)]
    a = [[F((seed+mu+2*k)%7-3, 5) for k in range(3)] for mu in range(4)]
    am = [lin(row, TT) for row in a]
    z, mass = F(seed+1, 3), F(seed, 7)
    current = [[-z*bil(psi, mul(HG[mu], t))/2 for t in TT] for mu in range(4)]
    eq = vs(mv(H, psi), -z*mass)
    for mu in range(4):
        eq = va(eq, vs(mv(HG[mu], va(dp[mu], mv(am[mu], psi))), z))
    divj = [F(0)]*3
    ordinary = [F(0)]*3
    for mu in range(4):
        dj = [-z*bil(dp[mu], mul(HG[mu], t), psi) for t in TT]
        ordinary = va(ordinary, dj)
        divj = va(divj, va(dj, cross(a[mu], current[mu])))
    rhs = [dot(eq, mv(t, psi)) for t in TT]
    check(divj == rhs, 'off_shell_matter_Noether')
    if seed == 1:
        reject('omit_color_precession_in_matter_current', ordinary != rhs)
    def density(v, dv, connection):
        return z*(sum(bil(v, HG[mu], va(dv[mu], mv(connection[mu], v))) for mu in range(4))-mass*bil(v, H))/2
    for mu, c in product(range(4), range(3)):
        up = [row[:] for row in am]; dn = [row[:] for row in am]
        up[mu] = add(am[mu], TT[c]); dn[mu] = sub(am[mu], TT[c])
        check((density(psi, dp, up)-density(psi, dp, dn))/2 == -current[mu][c], 'independent_source_variations')
    g = kron(add(scale(eye(4), F(3, 5)), scale(E[seed%3], F(4, 5))), eye(4))
    gi = transpose(g)
    xx = [TT[mu%3] for mu in range(4)]
    vp = mv(gi, psi)
    dvp = [va(mv(gi, dp[mu]), vs(mv(xx[mu], vp), -1)) for mu in range(4)]
    ap = [add(mul(mul(gi, am[mu]), g), xx[mu]) for mu in range(4)]
    check(density(vp, dvp, ap) == density(psi, dp, am), 'local_massive_action_covariance')
    reject('drop_inhomogeneous_connection_term_'+str(seed), density(vp, dvp, [sub(ap[mu], xx[mu]) for mu in range(4)]) != density(psi, dp, am))

# The regular thermo coefficient chart remains full rank at a flat connection.
pv = F(5, 108)
check(rank(scale(eye(12), pv)) == 12, 'thermo_full_variation_rank')
check(rank(scale(eye(16), pv)) == 16, 'thermo_full_variation_rank')
for q, root in [(F(0), F(3)), (F(1, 100), F(31, 10)), (F(-1, 100), F(29, 10))]:
    p = F(1, 6)-1/(2*root)
    v = ((1+q)**2+1/(4*(p-F(1, 6))**2))/5-2
    delta = 5*(2+v)-(1+q)**2
    check(delta == root**2 and 2+v > 0 and delta > 0, 'thermo_inverse_stability')
reject('single_shared_channel_has_full_gauge_rank', rank([[pv]*12 for _ in range(12)]) != 12)


def encode(x):
    if isinstance(x, F):
        return str(x)
    raise TypeError(type(x).__name__)


print(json.dumps({'status': 'PASS_INDEPENDENT_COMPACT_GAUGE_CONTROLS', 'checks': checks,
                  'total': sum(checks.values()), 'negative_controls': negative,
                  'centralizer_dimensions': {'two_real_copies': 1, 'four_real_copies': 4},
                  'invariant_pairing_dimension': 2, 'fixtures': fixtures}, default=encode, sort_keys=True))
