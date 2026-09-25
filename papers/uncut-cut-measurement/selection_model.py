"""Conditional coupling and coframe constraints; no physical sector selector.

Vertices are declared response components, not spacetime points. A smooth
chart enters only in the separately supplied homogeneous coframe adapter.
"""
from fractions import Fraction as Q

import gravity_model as gm


def coupling_edges(n, edges):
    if type(n) is not int or n < 2:
        raise ValueError('At least two response components are required')
    edges = tuple(tuple(edge) for edge in edges)
    seen = set()
    for edge in edges:
        if len(edge) != 3:
            raise ValueError('Each edge is (i, j, sector)')
        i, j, sector = edge
        if (type(i) is not int or type(j) is not int or not 0 <= i < j < n
                or sector not in (-1, 1) or (i, j) in seen):
            raise ValueError('Use a simple graph with i < j and sector +/-1')
        seen.add((i, j))
    return edges


def pair_generator(n, i, j, sector):
    """sector +1: K; sector -1: R. Upper entry is sector, lower is one."""
    coupling_edges(n, ((i, j, sector),))
    rows = [[0] * n for _ in range(n)]
    rows[i][j], rows[j][i] = sector, 1
    return gm.matrix(rows)


def coupling_form(n, edges):
    """Connected graph: normalized invariant form, or None for cycle conflict.

    This graph propagation algorithm supplies no proposed metric to a solver.
    Disconnected graphs fall outside U21 and are rejected, not classified.
    """
    edges = coupling_edges(n, edges)
    neighbours = [[] for _ in range(n)]
    for i, j, sector in edges:
        neighbours[i].append((j, -sector))
        neighbours[j].append((i, -sector))
    signs, queue, conflict = {0: 1}, [0], False
    for i in queue:
        for j, ratio in neighbours[i]:
            value = ratio * signs[i]
            if j not in signs:
                signs[j] = value
                queue.append(j)
            elif signs[j] != value:
                conflict = True
    if len(signs) != n:
        raise ValueError('U21 requires a connected coupling graph')
    if conflict:
        return None
    return gm.matrix([[signs[i] if i == j else 0 for j in range(n)]
                      for i in range(n)])


def diagonal_signature(form):
    form = gm.matrix(form)
    n = len(form)
    if len(form[0]) != n or any(form[i][j] for i in range(n)
                                for j in range(n) if i != j):
        raise ValueError('A diagonal form is required')
    return tuple(sum((form[i][i] > 0, form[i][i] < 0, form[i][i] == 0)[k]
                     for i in range(n)) for k in range(3))


def homogeneous_data(n, a, b, db, kappa, sector):
    """Pointwise data for e=(a du,b dx^i), A_0=0, A_i=kappa G_(0i).

    a,b are nonzero values; db is db/du. All functions depend only on u.
    No derivatives of a enter de because da wedge du is zero.
    """
    coupling_edges(n, ((0, 1, sector),))
    a, b, db, kappa = map(Q, (a, b, db, kappa))
    if not a or not b or not kappa:
        raise ValueError('Nonzero lapse, coframe scale and coupling required')
    zero = gm.scale(gm.identity(n), 0)
    e = gm.matrix([[a if i == j == 0 else b if i == j else 0
                    for j in range(n)] for i in range(n)])
    de = gm.matrix([[db if i == j and i else 0 for j in range(n)]
                    for i in range(n)])
    connection = (zero,) + tuple(gm.scale(pair_generator(n, 0, i, sector), kappa)
                                  for i in range(1, n))
    h = coupling_form(n, tuple((0, i, sector) for i in range(1, n)))
    torsion = {}
    for mu in range(n):
        for nu in range(mu + 1, n):
            left, right = gm.mul(connection[mu], e), gm.mul(connection[nu], e)
            torsion[f'{mu},{nu}'] = tuple(
                (de[r][nu] if mu == 0 else 0)
                - (de[r][mu] if nu == 0 else 0)
                + left[r][nu] - right[r][mu] for r in range(n))
    return {'coframe': e, 'internal_form': h, 'metric': gm.congruence(e, h),
            'connection': connection, 'torsion': torsion}
