"""Finite exact evidence for U21-U23; proofs and coverage are separate."""
from fractions import Fraction as Q
from itertools import combinations, product

import gravity_model as gm
import selection_model as sm
from gravity_certificate import in_span, rational_records


def require(value, message):
    if not value:
        raise RuntimeError(message)


def metric_curvature(n, b, db, ddb, sector):
    """Independent generic Christoffel formula from a diagonal metric's jets.

    All derivatives except u derivatives vanish. Does not read a connection,
    a coframe, or the expected closed-form curvature.
    """
    b, db, ddb = map(Q, (b, db, ddb))
    if not b:
        raise ValueError('Degenerate metric')
    g = gm.matrix([[1 if i == j == 0 else -sector * b*b if i == j else 0
                    for j in range(n)] for i in range(n)])
    dg = gm.matrix([[-2*sector*b*db if i == j and i else 0
                     for j in range(n)] for i in range(n)])
    ddg = gm.matrix([[-2*sector*(db*db+b*ddb) if i == j and i else 0
                      for j in range(n)] for i in range(n)])
    inv = gm.inverse(g)
    dinv = gm.scale(gm.mul(gm.mul(inv, dg), inv), -1)

    def jet(order, mu, i, j):
        return (dg if order == 1 else ddg)[i][j] if mu == 0 else Q(0)

    gamma, derivative = [], []
    for mu in range(n):
        rows, drows = [], []
        for a in range(n):
            row, drow = [], []
            for v in range(n):
                value = dvalue = Q(0)
                for r in range(n):
                    first = jet(1, mu, r, v)+jet(1, v, r, mu)-jet(1, r, mu, v)
                    second = jet(2, mu, r, v)+jet(2, v, r, mu)-jet(2, r, mu, v)
                    value += inv[a][r]*first/2
                    dvalue += (dinv[a][r]*first+inv[a][r]*second)/2
                row.append(value)
                drow.append(dvalue)
            rows.append(row)
            drows.append(drow)
        gamma.append(gm.matrix(rows))
        derivative.append(gm.matrix(drows))
    zero = gm.scale(gm.identity(n), 0)
    curvature = {}
    for mu, nu in combinations(range(n), 2):
        differential = gm.sub(derivative[nu] if mu == 0 else zero,
                              derivative[mu] if nu == 0 else zero)
        curvature[f'{mu},{nu}'] = gm.add(differential,
                                        gm.commutator(gamma[mu], gamma[nu]))
    return {'metric': g, 'gamma': tuple(gamma), 'curvature': curvature}


def controls():
    graphs = {
        'K_star': ((0, 1, 1), (0, 2, 1), (0, 3, 1)),
        'K_path': ((0, 1, 1), (1, 2, 1), (2, 3, 1)),
        'R_star': ((0, 1, -1), (0, 2, -1), (0, 3, -1)),
    }
    signatures = {name: sm.diagonal_signature(sm.coupling_form(4, edges))
                  for name, edges in graphs.items()}
    require(signatures == {'K_star': (1, 3, 0), 'K_path': (2, 2, 0),
                           'R_star': (4, 0, 0)}, 'W17 signature tie failed')
    require(sm.coupling_form(3, ((0, 1, 1), (1, 2, 1), (0, 2, 1))) is None,
            'W18 inconsistent triangle accepted')
    require(sm.diagonal_signature(sm.coupling_form(3,
            ((0, 1, 1), (0, 2, 1), (1, 2, -1)))) == (1, 2, 0),
            'W18 balanced triangle rejected')
    # One point can have zero torsion while the proposed function fails on a patch.
    defects = []
    for u in (Q(0), Q(1, 2), Q(1)):
        data = sm.homogeneous_data(3, 1, 1+u*u, 2*u, 1, 1)
        defects.append(data['torsion']['0,1'][1])
    require(defects == [-1, 0, 1], 'W19 fixed-connection rejection failed')
    k = sm.homogeneous_data(3, 1, 3, 1, 1, 1)
    r = sm.homogeneous_data(3, 1, 3, 1, 1, -1)
    require(all(not any(t) for d in (k, r) for t in d['torsion'].values()),
            'W20 both coframe closures must pass')
    require(k['metric'] != r['metric'], 'W20 signatures collapsed')
    return {'W17_minimum_edge_tie': {'edge_count_each': 3, 'signatures': signatures},
            'W18_cycle_gate': {'all_K_triangle_nonzero_form': False,
                              'two_K_one_R_signature': [1, 2, 0]},
            'W19_pointwise_zero_is_not_patch_closure': {
                'u': [0, Q(1, 2), 1], 'torsion_coefficients': defects},
            'W20_coframe_closure_does_not_select_sector': {
                'K_metric': k['metric'], 'R_metric': r['metric']}}


def build_checks():
    counts = {'labelled_graphs_enumerated': 0, 'connected_graphs': 0,
              'disconnected_graphs_outside_theorem': 0, 'balanced_graphs': 0,
              'inconsistent_graphs': 0, 'chart_points': 0,
              'curvature_component_comparisons': 0}
    per_dimension = []
    for n in range(2, 5):
        pairs = tuple(combinations(range(n), 2))
        connected = balanced = 0
        for labels in product((0, -1, 1), repeat=len(pairs)):
            edges = tuple((i, j, s) for (i, j), s in zip(pairs, labels) if s)
            counts['labelled_graphs_enumerated'] += 1
            try:
                form = sm.coupling_form(n, edges)
            except ValueError as exc:
                require(str(exc) == 'U21 requires a connected coupling graph',
                        'Unexpected invalid graph')
                counts['disconnected_graphs_outside_theorem'] += 1
                continue
            connected += 1
            # Oracle solves full symmetric-form equations for finite Cayley maps;
            # it does not read propagated signs or restrict to diagonal forms.
            arrows = tuple(gm.cayley(sm.pair_generator(n, i, j, s), Q(1, 2))
                           for i, j, s in edges)
            basis = gm.invariant_forms(n, arrows)
            if form is None:
                require(not basis, 'U21 conflicting cycle admitted a form')
                counts['inconsistent_graphs'] += 1
            else:
                require(len(basis) == 1 and in_span(form, basis),
                        'U21 propagation disagrees with full invariant solver')
                require(all(gm.congruence(u, form) == form for u in arrows),
                        'U21 form failed a finite transport')
                balanced += 1
                counts['balanced_graphs'] += 1
        counts['connected_graphs'] += connected
        per_dimension.append({'dimension': n, 'connected': connected,
                              'balanced': balanced})

    for n, sector, kappa, u in product(range(2, 6), (-1, 1),
                                      (Q(-1), Q(1, 2)), (Q(-1), Q(0), Q(2))):
        b = 3+kappa*u
        data = sm.homogeneous_data(n, 1, b, kappa, kappa, sector)
        oracle = metric_curvature(n, b, kappa, 0, sector)
        require(data['metric'] == oracle['metric'], 'U23 metric route mismatch')
        require(all(not any(t) for t in data['torsion'].values()), 'U23 torsion')
        zero = gm.scale(gm.identity(n), 0)
        for a in data['connection']:
            require(gm.add(gm.mul(gm.transpose(a), data['internal_form']),
                           gm.mul(data['internal_form'], a)) == zero,
                    'U23 metric compatibility failed')
        for mu, nu in combinations(range(n), 2):
            f = gm.commutator(data['connection'][mu], data['connection'][nu])
            represented = gm.mul(gm.mul(gm.inverse(data['coframe']), f), data['coframe'])
            actual = oracle['curvature'][f'{mu},{nu}']
            require(represented == actual, 'U23 independent LC curvature mismatch')
            if mu == 0:
                require(actual == zero, 'U23 mixed curvature must vanish')
            else:
                expected = [[Q(0)]*n for _ in range(n)]
                expected[mu][nu], expected[nu][mu] = sector*kappa*kappa, -sector*kappa*kappa
                require(actual == gm.matrix(expected), 'U23 spatial curvature sign')
            counts['curvature_component_comparisons'] += 1
        counts['chart_points'] += 1
    return rational_records({'counts': counts, 'graphs_by_dimension': per_dimension,
        'controls': controls(),
        'coverage': 'All absent/K/R simple labelled graphs on 2-4 vertices; connected graphs compared to full Cayley invariant-form nullspace. Charts: n=2..5, both sectors, kappa=-1 or 1/2, b=3+kappa*u, u=-1,0,2.',
        'physical_sector_selection': False,
        'coframe_profile_constrained_within_supplied_ansatz': True,
        'chart_and_homogeneity_derived': False})
