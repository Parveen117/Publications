"""Exact finite checks for native grading obstruction and linear seam memory."""
from fractions import Fraction as Q
from itertools import combinations, product

import gravity_model as gm
import grading_model as gr
import selection_model as sm
from gravity_certificate import rational_records


def require(value, message):
    if not value:
        raise RuntimeError(message)


def apply(a, x):
    return tuple(sum(ai*xi for ai, xi in zip(row, x)) for row in a)


def controls():
    edges = ((0, 1), (0, 2))
    star = gr.port_datum(3, edges)
    c, j = star['cut'], star['grading']
    z = (1, 0, -1, 0)
    require(apply(c, z) == (0, 0, 0) and apply(gm.mul(c, j), z) == (0, 1, -1),
            'W22 overlap memory witness failed')
    require(gm.descend_transport(j, c) is None, 'W22 grading falsely descends')
    repaired = gr.minimum_grading_observer(c, j)
    require(repaired['extra_scalar_channels'] == 1, 'W23 tree memory cost')

    triangle = gr.port_datum(3, ((0, 1), (0, 2), (1, 2)))
    ct, jt = triangle['cut'], triangle['grading']
    hidden = (1, -1, -1, 1, 1, -1)
    require(apply(triangle['target'], hidden) == (0,)*6, 'W24 target blindness')
    kstep, rstep = (gr.edge_step(3, 0, Q(1, 2), s) for s in (1, -1))
    kout, rout = (apply(gm.mul(ct, u), hidden) for u in (kstep, rstep))
    require(kout == (Q(-2, 3), Q(2, 3), 0), 'W25 K transcript')
    require(rout == (Q(2, 5), Q(6, 5), 0), 'W25 R transcript')
    require(gr.rank(triangle['target']) == 5 and
            gr.rank(gr.response_rows(ct, [gr.edge_step(3, e, Q(1, 2))
                                         for e in range(3)])) == 6,
            'W24 grading target mistaken for all edge responses')

    # Canonical averaging section for C; it is not an exact descent licence.
    degree_inverse = gm.matrix([[Q(1, d) if i == k else 0
                                 for k in range(3)]
                                for i, d in enumerate(star['degrees'])])
    section = gm.mul(gm.transpose(c), degree_inverse)
    u, v = (gr.edge_step(2, e, Q(1, 2)) for e in range(2))
    cu, cv = (gm.mul(gm.mul(c, a), section) for a in (u, v))
    require(gm.mul(u, v) == gm.mul(v, u), 'W26 block transports should commute')
    require(gm.mul(cu, cv) != gm.mul(cv, cu), 'W26 compression artifact missing')
    require(all(gm.descend_transport(a, c) is None for a in (u, v)),
            'W26 invalid compression passed exact descent')
    h = gm.matrix(((1, 0), (0, -1)))
    k = sm.pair_generator(2, 0, 1, 1)
    r = sm.pair_generator(2, 0, 1, -1)
    require(gm.mul(gm.mul(k, k), k) == k and gm.mul(gm.mul(k, r), k) == gm.scale(r, -1),
            'W21 native grading fixture failed')
    require(gm.mul(gm.mul(h, k), h) == gm.scale(k, -1),
            'W21 metric-sign involution must differ from native grading')
    return {'W21_two_involutions_are_not_identical': {'native_J': k, 'metric_sign_J': h},
            'W22_gluing_loses_grading': {'hidden_state': z, 'initial': apply(c, z),
                                       'after_J': apply(gm.mul(c, j), z)},
            'W23_three_vertex_tree': {'source_rank': 4, 'cut_rank': 3, 'minimum_extra': 1},
            'W24_triangle_targets_differ': {'hidden_state': hidden,
                'grading_observer_rank': 5, 'edge_response_observer_rank': 6},
            'W25_declared_sector_transcripts': {'parameter': Q(1, 2),
                'K_algebraic_candidate': kout, 'R_native_sector': rout,
                'linfinity_separation': Q(16, 15)},
            'W26_compression_invents_commutator': {'native_commutes': True,
                'compressed_commutator': gm.commutator(cu, cv), 'exact_descent': False}}


def build_checks():
    counts = {'graphs_enumerated': 0, 'connected_graphs': 0,
              'bipartite_connected': 0, 'nonbipartite_connected': 0,
              'full_grading_linear_systems': 0, 'minimal_observers_constructed': 0,
              'addressed_R_families': 0, 'local_native_grading_checks': 0,
              'edge_response_reconstructions': 0, 'cut_memory_rank_checks': 0}
    dimensions = []
    for n in range(2, 6):
        pairs = tuple(combinations(range(n), 2))
        connected = 0
        for flags in product((False, True), repeat=len(pairs)):
            edges = tuple(e for e, chosen in zip(pairs, flags) if chosen)
            counts['graphs_enumerated'] += 1
            if not gr.graph_data(n, edges)['connected']:
                continue
            data = gr.port_datum(n, edges)
            c, j, target = data['cut'], data['grading'], data['target']
            m, size = len(edges), 2*len(edges)
            connected += 1
            counts['connected_graphs'] += 1
            counts['bipartite_connected' if data['bipartite'] else 'nonbipartite_connected'] += 1
            expected_rank = 2*n-2 if data['bipartite'] else 2*n-1
            require(gr.rank(c) == n and gr.rank(target) == expected_rank,
                    'U25 rank formula disagrees with rational elimination')
            # This rank is computed from raw C,CJ, not from graph parity.
            if n > 4:
                continue
            basis = gr.common_graders(n, edges)
            require(len(basis) == (1 if n == 2 else 0), 'U24 common grading obstruction')
            counts['full_grading_linear_systems'] += 1
            repaired = gr.minimum_grading_observer(c, j)
            a, dj = repaired['observer'], repaired['descended_grading']
            require(len(a) == expected_rank and gm.mul(dj, dj) == gm.identity(len(a)),
                    'U25 minimal observer or descended involution failed')
            require(gm.mul(dj, a) == gm.mul(a, j), 'U25 grading intertwiner')
            # Every independent added scalar channel is necessary for this target.
            for removed in range(n, len(a)):
                reduced = a[:removed]+a[removed+1:]
                require(gr.rank(reduced+target) > gr.rank(reduced),
                        'U25 claimed repair includes a redundant channel')
            counts['minimal_observers_constructed'] += 1
            memory = gr.grading_memory_lift(c, j)
            require(gr.rank(memory) == expected_rank-n,
                    'U25 cut-return memory rank disagrees with observer minimum')
            counts['cut_memory_rank_checks'] += 1
            steps = [gr.edge_step(m, e, Q(1, 2)) for e in range(m)]
            require(gr.rank(gr.response_rows(c, steps)) == size,
                    'U26 addressed R responses fail to observe full source')
            counts['addressed_R_families'] += 1
            for e, (i, k) in enumerate(edges):
                # Direct column permutation implements J G J independently of mul.
                for sector in (1, -1):
                    g = gr.edge_generator(m, e, sector)
                    conjugated = tuple(tuple(g[r ^ 1][s ^ 1] for s in range(size))
                                       for r in range(size))
                    require(conjugated == gm.scale(g, sector), 'U24 local grading')
                    counts['local_native_grading_checks'] += 1
                # Recover each endpoint pair from the initial/one-step difference.
                delta = gm.sub(gm.mul(c, steps[e]), c)
                block = tuple(tuple(delta[v][2*e+r] for r in range(2)) for v in (i, k))
                selector = gm.mul(gm.inverse(block), (delta[i], delta[k]))
                expected = tuple(tuple(Q(col == 2*e+r) for col in range(size))
                                 for r in range(2))
                require(selector == expected, 'U26 local response inverse failed')
                counts['edge_response_reconstructions'] += 1
        dimensions.append({'vertices': n, 'connected_graphs': connected})
    return rational_records({'counts': counts, 'by_dimension': dimensions,
        'controls': controls(),
        'coverage': 'All simple labelled graphs on 2-5 vertices for C/CJ rank; all connected graphs on 2-4 vertices for full grading equations, minimal observer construction and separately addressed R Cayley responses at t=1/2.',
        'minimum_cost_type': 'Additional independent linear scalar channels; not bits, entropy or physical action',
        'physical_instrument_or_gravity_identification': False,
        'native_multi_pair_interaction_selected': False})
