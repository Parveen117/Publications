"""Exact fixtures for this publication; not a replacement native operator engine."""
from fractions import Fraction as Q
from itertools import product


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def zero(n, m=None):
    return [[Q(0) for _ in range(m or n)] for _ in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def add(a, b):
    return [[x+y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(c, a):
    return [[c*x for x in r] for r in a]


def mm(a, b):
    return [[sum(x*y for x, y in zip(r, s)) for s in zip(*b)] for r in a]


def mv(a, v):
    return [sum(x*y for x, y in zip(r, v)) for r in a]


def dot(u, v):
    return sum(x*y for x, y in zip(u, v))


def comm(a, b):
    return add(mm(a, b), scale(-1, mm(b, a)))


def anti(a, b):
    return add(mm(a, b), mm(b, a))


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def kron(a, b):
    return [[x*y for x in ar for y in br] for ar in a for br in b]


def rref(a):
    a = [[Q(x) for x in r] for r in a]
    row = 0
    pivots = []
    for col in range(len(a[0])):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        c = a[row][col]
        a[row] = [x/c for x in a[row]]
        for i in range(len(a)):
            if i != row:
                c = a[i][col]
                a[i] = [x-c*y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    return a, pivots


def rank(a):
    return len(rref(a)[1])


def inverse(a):
    n = len(a)
    rows, _ = rref([list(r)+s for r, s in zip(a, eye(n))])
    if [r[:n] for r in rows] != eye(n):
        raise ValueError('singular matrix')
    return [r[n:] for r in rows]


def det(a):
    a = [[Q(x) for x in r] for r in a]
    answer = Q(1)
    for j in range(len(a)):
        p = next((i for i in range(j, len(a)) if a[i][j]), None)
        if p is None:
            return Q(0)
        if p != j:
            a[j], a[p] = a[p], a[j]
            answer = -answer
        c = a[j][j]
        answer *= c
        for i in range(j+1, len(a)):
            s = a[i][j]/c
            a[i] = [x-s*y for x, y in zip(a[i], a[j])]
    return answer


def fixture():
    i = eye(2)
    r = [[Q(0), Q(-1)], [Q(1), Q(0)]]
    k = [[Q(1), Q(0)], [Q(0), Q(-1)]]
    l = mm(k, r)
    gamma = [kron(r, i), kron(k, k), kron(k, l), kron(l, i)]
    a = [mm(gamma[0], g) for g in gamma[1:]]
    j = mm(mm(a[0], a[1]), a[2])
    return (i, r, k, l), gamma, a, j


def derivative(poly, point, indices=()):
    """Evaluate exact derivatives of a (possibly Laurent) polynomial."""
    result = Q(0)
    for powers, coeff in poly.items():
        powers = list(powers)
        coeff = Q(coeff)
        for i in indices:
            coeff *= powers[i]
            powers[i] -= 1
        if coeff:
            for x, exponent in zip(point, powers):
                coeff *= Q(x)**exponent
            result += coeff
    return result


def ricci(metric, point):
    """Direct real-coordinate Levi-Civita calculation from metric derivatives."""
    n = len(metric)
    ix = range(n)
    g = [[derivative(p, point) for p in row] for row in metric]
    gi = inverse(g)
    dg = [[[derivative(metric[a][b], point, (c,)) for b in ix] for a in ix] for c in ix]
    ddg = [[[[derivative(metric[a][b], point, (c, d)) for b in ix]
             for a in ix] for d in ix] for c in ix]
    dgi = [scale(-1, mm(mm(gi, d), gi)) for d in dg]
    gamma = [[[sum(gi[k][l]*(dg[a][b][l]+dg[b][a][l]-dg[l][a][b])
                   for l in ix)/2 for b in ix] for a in ix] for k in ix]

    def dgamma(c, k, a, b):
        return sum(dgi[c][k][l]*(dg[a][b][l]+dg[b][a][l]-dg[l][a][b])
                   +gi[k][l]*(ddg[c][a][b][l]+ddg[c][b][a][l]-ddg[c][l][a][b])
                   for l in ix)/2

    ric = [[sum(dgamma(k, k, a, b)-dgamma(b, k, a, k) for k in ix)
            +sum(gamma[k][k][l]*gamma[l][a][b]-gamma[k][b][l]*gamma[l][a][k]
                 for k, l in product(ix, repeat=2)) for b in ix] for a in ix]
    scalar = sum(gi[a][b]*ric[a][b] for a, b in product(ix, repeat=2))
    return g, ric, scalar


def centered_periodic(n):
    d = zero(n)
    for i in range(n):
        d[i][(i+1) % n] += Q(1, 2)
        d[i][(i-1) % n] -= Q(1, 2)
    return d


def cayley(generator, state, step):
    i = eye(len(state))
    return mv(inverse(add(i, scale(-step/2, generator))),
              mv(add(i, scale(step/2, generator)), state))
