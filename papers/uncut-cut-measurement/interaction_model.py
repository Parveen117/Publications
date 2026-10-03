"""Declared real inter-edge mixing, exact linear cuts, and response protocols.

Couplings obey the native real anti-self-dagger and EMK cut-grading rules.
The interaction graph and strengths are supplied, not physically selected.
"""
from fractions import Fraction as Q

import gravity_model as gm
import grading_model as gr
import selection_model as sm


def links_for(m, links):
    if type(m) is not int or m < 1:
        raise ValueError('A positive number of two-component channels is required')
    links, seen = tuple(tuple(e) for e in links), set()
    for edge in links:
        if (len(edge) != 2 or any(type(x) is not int for x in edge)
                or not 0 <= edge[0] < edge[1] < m or edge in seen):
            raise ValueError('Use a simple interaction graph with a < b')
        seen.add(edge)
    return links


def components(m, links):
    links = links_for(m, links)
    neighbours = [[] for _ in range(m)]
    for a, b in links:
        neighbours[a].append(b)
        neighbours[b].append(a)
    seen, result = set(), []
    for root in range(m):
        if root in seen:
            continue
        seen.add(root)
        queue = [root]
        for a in queue:
            for b in sorted(neighbours[a]):
                if b not in seen:
                    seen.add(b)
                    queue.append(b)
        result.append(tuple(sorted(queue)))
    return tuple(result)


def kron(a, b):
    a, b = gm.matrix(a), gm.matrix(b)
    return gm.matrix([[a[i][j]*b[r][s] for j in range(len(a[0]))
                       for s in range(len(b[0]))]
                      for i in range(len(a)) for r in range(len(b))])


def native_grading(m):
    links_for(m, ())
    return kron(gm.identity(m), sm.pair_generator(2, 0, 1, 1))


def mixing_generator(m, a, b):
    links_for(m, ((a, b),))
    return kron(sm.pair_generator(m, a, b, -1), gm.identity(2))


def mixing_step(m, a, b, parameter):
    links_for(m, ((a, b),))
    q = Q(parameter)
    if not q:
        raise ValueError('Mixing parameter must be nonzero')
    c, s = (1-q*q)/(1+q*q), 2*q/(1+q*q)
    rows = [list(row) for row in gm.identity(2*m)]
    for r in range(2):
        i, j = 2*a+r, 2*b+r
        rows[i][i] = rows[j][j] = c
        rows[i][j], rows[j][i] = -s, s
    return gm.matrix(rows)


def family(m, links, phase=Q(1, 2), mixing=Q(1, 2)):
    links = links_for(m, links)
    return (tuple(gr.edge_step(m, a, phase) for a in range(m))
            + tuple(mixing_step(m, a, b, mixing) for a, b in links))


def projection(m, channels):
    links_for(m, ())
    channels = tuple(channels)
    if (not channels or len(set(channels)) != len(channels)
            or any(type(a) is not int or not 0 <= a < m for a in channels)):
        raise ValueError('Select distinct channels from the carrier')
    unit = gm.identity(2*m)
    return tuple(unit[2*a+r] for a in channels for r in range(2))


def predicted_observer(m, links, cut):
    """Graph oracle for the real-linear closure theorem; nonzero readout only."""
    cut = gm.matrix(cut)
    if len(cut[0]) != 2*m or gr.rank(cut) == 0:
        raise ValueError('Supply a nonzero readout on the declared real carrier')
    active = []
    for comp in components(m, links):
        if any(row[2*a+r] for row in cut for a in comp for r in range(2)):
            active.extend(comp)
    observed = projection(m, sorted(active))
    return {'observer': observed, 'rank': len(observed),
            'extra_scalar_channels': len(observed)-gr.rank(cut)}


def row_closure(cut, steps):
    """Independent matrix observability closure; knows no interaction graph."""
    cut = gm.matrix(cut)
    steps = tuple(gm.matrix(u) for u in steps)
    n = len(cut[0])
    if any((len(u), len(u[0])) != (n, n) for u in steps):
        raise ValueError('Steps and readout use different carriers')
    reduced, pivots = gm.rref(cut)
    if not pivots:
        raise ValueError('Supply a nonzero readout')
    basis, initial = reduced[:len(pivots)], len(pivots)
    for depth in range(n-initial+1):
        rows = basis+tuple(row for u in steps for row in gm.mul(basis, u))
        reduced, pivots = gm.rref(rows)
        if len(pivots) == len(basis):
            return {'observer': basis, 'depth': depth,
                    'extra_scalar_channels': len(basis)-initial}
        basis = reduced[:len(pivots)]
    raise RuntimeError('Finite row-space closure exceeded its rank bound')


def invariant_form_basis(m, links):
    """Predicted full symmetric-form space: one scalar identity per component."""
    basis = []
    for comp in components(m, links):
        selected = set(comp)
        basis.append(gm.matrix([[int(i == j and i//2 in selected)
                                  for j in range(2*m)] for i in range(2*m)]))
    return tuple(basis)


def optimal_transcript(m, links, root=0, phase=Q(1, 2), mixing=Q(1, 2)):
    """Two scalar response rows per channel via a declared rooted spanning tree.

    row_product lists matrices multiplied after the readout row; chronological
    application to the input is in the reverse order. The matrix rows and the
    explicit chronological words are both returned to avoid convention drift.
    """
    links = links_for(m, links)
    if type(root) is not int or not 0 <= root < m:
        raise ValueError('Invalid observation root')
    if len(components(m, links)) != 1:
        raise ValueError('The optimal full-carrier transcript requires connected mixing')
    phase, mixing = Q(phase), Q(mixing)
    if not phase or not mixing:
        raise ValueError('Both protocol parameters must be nonzero')
    neighbours = [[] for _ in range(m)]
    for a, b in links:
        neighbours[a].append(b)
        neighbours[b].append(a)
    rows_at = {root: (gm.identity(2*m)[2*root],)}
    paths, queue = {root: ()}, [root]
    for a in queue:
        for b in sorted(neighbours[a]):
            if b in paths:
                continue
            edge = (min(a, b), max(a, b))
            paths[b] = paths[a]+(f'B{edge[0]},{edge[1]}',)
            rows_at[b] = gm.mul(rows_at[a], mixing_step(m, *edge, mixing))
            queue.append(b)
    records, rows, words = [], [], []
    for a in queue:
        row = rows_at[a][0]
        turned = gm.mul((row,), gr.edge_step(m, a, phase))[0]
        records.append({'channel': a, 'first_row': row, 'coefficient': row[2*a],
                        'row_product': paths[a]})
        rows.extend((row, turned))
        words.extend((tuple(reversed(paths[a])), tuple(reversed(paths[a]+(f'R{a}',)))))
    return {'rows': tuple(rows), 'records': tuple(records), 'chronological_words': tuple(words),
            'phase': phase, 'mixing': mixing, 'root': root}


def decode_transcript(packet, values):
    """Triangular decoder, independent of generic matrix inversion."""
    records, values = packet['records'], tuple(Q(x) for x in values)
    if len(values) != 2*len(records):
        raise ValueError('Provide exactly two readings per declared channel')
    q = packet['phase']
    c, s = (1-q*q)/(1+q*q), 2*q/(1+q*q)
    result = [Q(0)] * len(values)
    for k, record in enumerate(records):
        a, row, alpha = record['channel'], record['first_row'], record['coefficient']
        first, turned = values[2*k:2*k+2]
        first -= sum(row[2*b]*result[2*b] for b in range(len(records)) if b != a)
        result[2*a] = first/alpha
        result[2*a+1] = ((c-1)*alpha*result[2*a]-(turned-values[2*k]))/(s*alpha)
    return tuple(result)
