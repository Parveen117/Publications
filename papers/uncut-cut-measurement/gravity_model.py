"""Exact conditional transport-to-geometry bridge; no physical gravity solver.

All arithmetic is rational. Real K-sector Cayley maps are candidate invertible
maps, not the dagger-unitary flows of native T55. No metric is given to the
invariant-form solver. Smooth charts are separate, explicitly supplied adapters.
"""
from fractions import Fraction as Q


def matrix(rows):
    rows = tuple(tuple(Q(x) for x in row) for row in rows)
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError('Supply a nonempty rectangular rational matrix')
    return rows


def identity(n):
    if type(n) is not int or n < 1:
        raise ValueError('Positive integer dimension required')
    return matrix(tuple(tuple(int(i == j) for j in range(n)) for i in range(n)))


def transpose(a):
    return tuple(zip(*matrix(a)))


def mul(a, b):
    a, b = matrix(a), matrix(b)
    if len(a[0]) != len(b):
        raise ValueError('Matrix composition has incompatible dimensions')
    return tuple(tuple(sum(x * y for x, y in zip(row, col))
                       for col in transpose(b)) for row in a)


def scale(a, factor):
    return tuple(tuple(Q(factor) * x for x in row) for row in matrix(a))


def add(a, b):
    a, b = matrix(a), matrix(b)
    if (len(a), len(a[0])) != (len(b), len(b[0])):
        raise ValueError('Matrix addition has incompatible dimensions')
    return tuple(tuple(x + y for x, y in zip(ar, br)) for ar, br in zip(a, b))


def sub(a, b):
    return add(a, scale(b, -1))


def rref(a):
    a = [list(row) for row in matrix(a)]
    pivots, row = [], 0
    for col in range(len(a[0])):
        pivot = next((k for k in range(row, len(a)) if a[k][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        divisor = a[row][col]
        a[row] = [x / divisor for x in a[row]]
        for k in range(len(a)):
            if k != row:
                factor = a[k][col]
                a[k] = [x - factor * y for x, y in zip(a[k], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    return matrix(a), tuple(pivots)


def nullspace(a):
    reduced, pivots = rref(a)
    free = tuple(i for i in range(len(reduced[0])) if i not in pivots)
    result = []
    for j in free:
        v = [Q(0)] * len(reduced[0])
        v[j] = Q(1)
        for row, pivot in enumerate(pivots):
            v[pivot] = -reduced[row][j]
        result.append(tuple(v))
    return tuple(result)


def inverse(a):
    a = matrix(a)
    n = len(a)
    if len(a[0]) != n:
        raise ValueError('Only square matrices can be inverted')
    aug, pivots = rref(tuple(a[i] + identity(n)[i] for i in range(n)))
    if pivots[:n] != tuple(range(n)):
        raise ValueError('Singular matrix')
    return tuple(row[n:] for row in aug)


def congruence(a, g):
    return mul(mul(transpose(a), g), a)


def commutator(a, b):
    return sub(mul(a, b), mul(b, a))


def generator(n, i, *, rotation=False):
    identity(n)
    if type(i) is not int or not 1 <= i < n:
        raise ValueError('Pair index must lie strictly between zero and n')
    a = [[0] * n for _ in range(n)]
    a[0][i], a[i][0] = (-1 if rotation else 1), 1
    return matrix(a)


def cayley(g, t):
    g = matrix(g)
    unit = identity(len(g))
    return mul(add(unit, scale(g, t)), inverse(sub(unit, scale(g, t))))


def emk_step(n, i, t, *, rotation=False):
    """Closed rational Cayley block in the declared (0,i) plane."""
    generator(n, i, rotation=rotation)
    t = Q(t)
    denominator = 1 + t * t if rotation else 1 - t * t
    if not denominator:
        raise ValueError('Singular K-sector Cayley parameter')
    a = (1 - t * t if rotation else 1 + t * t) / denominator
    b = 2 * t / denominator
    result = [list(row) for row in identity(n)]
    result[0][0] = result[i][i] = a
    result[0][i] = -b if rotation else b
    result[i][0] = b
    return matrix(result)


def invariant_forms(n, transports):
    """Basis of all symmetric G with U^T G U=G; no proposed G is an input."""
    identity(n)
    transports = tuple(matrix(u) for u in transports)
    if any(len(u) != n or len(u[0]) != n for u in transports):
        raise ValueError('Transports must share the declared dimension')
    for u in transports:
        inverse(u)
    pairs = tuple((i, j) for i in range(n) for j in range(i, n))
    basis = []
    for i, j in pairs:
        e = [[0] * n for _ in range(n)]
        e[i][j] = e[j][i] = 1
        basis.append(matrix(e))
    rows = []
    for u in transports:
        defects = tuple(sub(congruence(u, e), e) for e in basis)
        rows.extend(tuple(d[i][j] for d in defects) for i, j in pairs)
    if not rows:
        rows = [tuple(Q(0) for _ in basis)]
    return tuple(tuple(tuple(sum(v[k] * basis[k][i][j] for k in range(len(basis)))
                             for j in range(n)) for i in range(n))
                 for v in nullspace(rows))


def right_inverse(cut):
    cut = matrix(cut)
    n, m = len(cut[0]), len(cut)
    _, columns = rref(cut)
    if len(columns) != m:
        raise ValueError('Cut must be surjective onto its declared target')
    block = tuple(tuple(row[j] for j in columns) for row in cut)
    block_inverse = inverse(block)
    section = [[Q(0)] * m for _ in range(n)]
    for i, col in enumerate(columns):
        section[col] = list(block_inverse[i])
    return matrix(section)


def descend_form(g, cut):
    """Unique h with C^T h C=g for a surjective cut, else None."""
    g, cut = matrix(g), matrix(cut)
    n = len(cut[0])
    if len(g) != n or len(g[0]) != n or transpose(g) != g:
        raise ValueError('A symmetric form on the source carrier is required')
    section = right_inverse(cut)
    target = congruence(section, g)
    return target if congruence(cut, target) == g else None


def descend_transport(u, cut):
    """Exact intertwiner V C=C U, or None; compression alone is insufficient."""
    cut = matrix(cut)
    target = mul(mul(cut, u), right_inverse(cut))
    return target if mul(target, cut) == mul(cut, u) else None


def metric_at_root(root_form, path_transport):
    return congruence(inverse(path_transport), root_form)


def torsion_2d(e, du_e, dv_e, a_u, a_v):
    """Internal components of de + A wedge e on the supplied chart (u,v)."""
    e, du_e, dv_e = matrix(e), matrix(du_e), matrix(dv_e)
    ae_u, ae_v = mul(a_u, e), mul(a_v, e)
    return tuple(du_e[i][1] - dv_e[i][0] + ae_u[i][1] - ae_v[i][0]
                 for i in range(len(e)))
