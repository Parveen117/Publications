"""Native EMK grading and exact linear memory for a declared edge assembly.

The edge carrier is a direct sum of native two-component representations.
It is a presented algebraic model, not an identification of the uncut ground.
"""
from fractions import Fraction as Q

import gravity_model as gm
import selection_model as sm


def graph_edges(n, edges):
    edges = tuple(tuple(edge) for edge in edges)
    sm.coupling_edges(n, tuple(tuple(edge)+(1,) for edge in edges))
    return edges


def graph_data(n, edges):
    edges = graph_edges(n, edges)
    adjacent = [[] for _ in range(n)]
    for i, j in edges:
        adjacent[i].append(j)
        adjacent[j].append(i)
    colours, components, bipartite = {}, 0, True
    for root in range(n):
        if root in colours:
            continue
        components += 1
        colours[root], queue = 1, [root]
        for i in queue:
            for j in adjacent[i]:
                if j not in colours:
                    colours[j] = -colours[i]
                    queue.append(j)
                elif colours[j] == colours[i]:
                    bipartite = False
    return {'connected': components == 1, 'bipartite': bipartite,
            'degrees': tuple(len(row) for row in adjacent)}


def rank(a):
    return len(gm.rref(a)[1])


def common_graders(n, edges):
    """Full linear solution space of JK=KJ and JR=-RJ in vertex coordinates.

    No diagonal, symmetry, or involution constraint is assumed by this oracle.
    A zero solution space therefore excludes every possible involution.
    """
    edges = graph_edges(n, edges)
    rows = []
    for i, j in edges:
        for sector in (1, -1):
            g = sm.pair_generator(n, i, j, sector)
            for a in range(n):
                for b in range(n):
                    row = [Q(0)] * (n*n)
                    for r in range(n):
                        row[a*n+r] += g[r][b]
                        row[r*n+b] -= sector*g[a][r]
                    if any(row):
                        rows.append(row)
    if not rows:
        rows = [[0] * (n*n)]
    return tuple(gm.matrix([v[i*n:(i+1)*n] for i in range(n)])
                 for v in gm.nullspace(rows))


def port_datum(n, edges):
    """C sums endpoint copies; J swaps the two copies belonging to each edge."""
    edges = graph_edges(n, edges)
    info = graph_data(n, edges)
    if not info['connected']:
        raise ValueError('The edge-memory theorem requires a connected graph')
    size = 2*len(edges)
    cut = [[0]*size for _ in range(n)]
    grading = [[0]*size for _ in range(size)]
    for e, (i, j) in enumerate(edges):
        cut[i][2*e] = cut[j][2*e+1] = 1
        grading[2*e][2*e+1] = grading[2*e+1][2*e] = 1
    c, j = gm.matrix(cut), gm.matrix(grading)
    return {'cut': c, 'grading': j, 'target': c+gm.mul(c, j),
            'edges': edges, **info}


def minimum_grading_observer(cut, grading):
    """Retain all initial rows and append independent next-grading readouts."""
    cut, grading = gm.matrix(cut), gm.matrix(grading)
    n = len(cut)
    if rank(cut) != n:
        raise ValueError('Initial observer must have independent rows')
    if gm.mul(grading, grading) != gm.identity(len(grading)):
        raise ValueError('Grading must be an involution')
    observer = cut
    for row in gm.mul(cut, grading):
        candidate = observer+(row,)
        if rank(candidate) > len(observer):
            observer = candidate
    descended = gm.descend_transport(grading, observer)
    if descended is None:
        raise RuntimeError('Constructed observer does not close under grading')
    return {'observer': observer, 'descended_grading': descended,
            'extra_scalar_channels': len(observer)-n}


def grading_memory_lift(cut, grading):
    """M=CJQ for the declared Euclidean port presentation; rank is target cost."""
    cut, grading = gm.matrix(cut), gm.matrix(grading)
    section = gm.mul(gm.transpose(cut), gm.inverse(gm.mul(cut, gm.transpose(cut))))
    blind = gm.sub(gm.identity(len(cut[0])), gm.mul(section, cut))
    return gm.mul(gm.mul(cut, grading), blind)


def edge_generator(edge_count, edge_index, sector):
    if (type(edge_count) is not int or edge_count < 1 or
            type(edge_index) is not int or not 0 <= edge_index < edge_count):
        raise ValueError('Valid edge index and positive edge count required')
    return sm.pair_generator(2*edge_count, 2*edge_index, 2*edge_index+1, sector)


def edge_step(edge_count, edge_index, t, sector=-1):
    """Block Cayley step; default R sector. K is an algebraic comparison only."""
    edge_generator(edge_count, edge_index, sector)
    t = Q(t)
    if not t:
        raise ValueError('A nonzero protocol parameter is required')
    block = gm.emk_step(2, 1, t, rotation=(sector == -1))
    rows = [list(row) for row in gm.identity(2*edge_count)]
    for i in range(2):
        for j in range(2):
            rows[2*edge_index+i][2*edge_index+j] = block[i][j]
    return gm.matrix(rows)


def response_rows(cut, steps):
    """Declared transcript: initial readout and each separately addressed step."""
    cut = gm.matrix(cut)
    return cut+tuple(row for step in steps for row in gm.mul(cut, step))
