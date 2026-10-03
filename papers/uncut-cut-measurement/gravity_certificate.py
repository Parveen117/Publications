"""Exact finite cross-checks of conditional transport/geometry propositions."""
from fractions import Fraction as Q
from itertools import product

import gravity_model as gm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def flat(a):
    return tuple(x for row in a for x in row)


def rational_records(value):
    if isinstance(value, Q):
        return str(value)
    if isinstance(value, dict):
        return {k: rational_records(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [rational_records(v) for v in value]
    return value


def in_span(target, basis):
    if not basis:
        return not any(flat(target))
    rows = tuple(flat(g) for g in basis)
    return len(gm.rref(rows)[1]) == len(gm.rref(rows + (flat(target),))[1])


def lorentz(n):
    return gm.matrix(tuple(tuple((1 if i == 0 else -1) if i == j else 0
                                 for j in range(n)) for i in range(n)))


def warped_chart(u, coefficients):
    """Independent pointwise routes for a supplied polynomial soldering form."""
    a, b, c = map(Q, coefficients)
    u = Q(u)
    f, df, ddf = a + b * u + c * u * u, b + 2 * c * u, 2 * c
    if not f:
        raise ValueError('Soldering form is singular at this chart point')
    e = gm.matrix(((1, 0), (0, f)))
    zero = gm.matrix(((0, 0), (0, 0)))
    du_e = gm.matrix(((0, 0), (0, df)))
    k = gm.generator(2, 1)
    av = gm.scale(k, df)
    fibre_curvature = gm.scale(k, ddf)
    g = gm.congruence(e, lorentz(2))
    # Independent Christoffel construction from g, inverse(g), and dg.
    derivatives = (gm.matrix(((0, 0), (0, -2 * f * df))), zero)
    gi = gm.inverse(g)
    gamma = tuple(gm.matrix(tuple(tuple(
        sum(gi[r][l] * (derivatives[i][l][j] + derivatives[j][l][i]
                        - derivatives[l][i][j]) / 2 for l in range(2))
        for j in range(2)) for r in range(2))) for i in range(2))
    # Derivative of the exact rational Christoffel entries, explicitly evaluated.
    du_gamma_v = gm.matrix(((0, df * df + f * ddf),
                            ((ddf * f - df * df) / (f * f), 0)))
    riemann = gm.add(du_gamma_v, gm.commutator(gamma[0], gamma[1]))
    represented = gm.mul(gm.mul(gm.inverse(e), fibre_curvature), e)
    require(riemann == represented, 'U20 connection curvature does not intertwine')
    require(gm.torsion_2d(e, du_e, zero, zero, av) == (0, 0), 'U20 nonzero torsion')
    require(gm.add(gm.mul(gm.transpose(av), lorentz(2)),
                   gm.mul(lorentz(2), av)) == zero, 'U20 incompatible connection')
    require(gamma[0] == gm.mul(gm.inverse(e), du_e), 'U20 u connection mismatch')
    require(gamma[1] == gm.mul(gm.mul(gm.inverse(e), av), e), 'U20 v connection mismatch')
    return {'warp': f, 'metric': g, 'fibre_curvature': fibre_curvature,
            'riemann_uv': riemann, 'torsion': (Q(0), Q(0))}


def build_checks():
    counts = {key: 0 for key in ('star_families', 'sector_comparisons',
        'cayley_two_route_checks', 'basis_change_checks', 'integer_single_transports',
        'candidate_form_comparisons', 'graph_edge_checks', 'cut_form_checks',
        'warped_chart_points', 'lawful_memory_cut_examples')}
    grid = (Q(-2), Q(-1, 2), Q(1, 3), Q(1, 2), Q(2))
    for n in range(2, 6):
        for t in grid:
            transports = tuple(gm.emk_step(n, i, t) for i in range(1, n))
            forms = gm.invariant_forms(n, transports)
            require(len(forms) == 1 and in_span(lorentz(n), forms), 'U18 signature selection')
            rotations = tuple(gm.emk_step(n, i, t, rotation=True) for i in range(1, n))
            rforms = gm.invariant_forms(n, rotations)
            require(len(rforms) == 1 and in_span(gm.identity(n), rforms), 'U18 rotation control')
            require(not gm.invariant_forms(n, transports + rotations), 'U18 incompatible sectors')
            for rotation, family in ((False, transports), (True, rotations)):
                for i, u in enumerate(family, 1):
                    require(u == gm.cayley(gm.generator(n, i, rotation=rotation), t),
                            'U18 closed block differs from generic inverse route')
                    require(gm.mul(u, gm.emk_step(n, i, -t, rotation=rotation)) == gm.identity(n),
                            'U18 inverse failed')
                    counts['cayley_two_route_checks'] += 1
            counts['star_families'] += 1
            counts['sector_comparisons'] += 1
        # Changing the frame cannot create the signature conclusion by convention.
        mixed = tuple(gm.emk_step(n, i, Q(i, i + 1)) for i in range(1, n))
        mixed_forms = gm.invariant_forms(n, mixed)
        require(len(mixed_forms) == 1 and in_span(lorentz(n), mixed_forms),
                'U18 unequal generator parameters failed')
        counts['star_families'] += 1
        family = tuple(gm.emk_step(n, i, Q(1, 2)) for i in range(1, n))
        for q in (Q(-1), Q(1), Q(3, 2)):
            s = [list(row) for row in gm.identity(n)]
            s[1][0] = q
            s = gm.matrix(s)
            si = gm.inverse(s)
            changed = tuple(gm.mul(gm.mul(s, u), si) for u in family)
            forms = gm.invariant_forms(n, changed)
            expected = gm.congruence(si, lorentz(n))
            require(len(forms) == 1 and in_span(expected, forms), 'U17 frame covariance')
            counts['basis_change_checks'] += 1

    # Exhaustive single-arrow and candidate-form oracle on an explicit integer grid.
    for entries in product((-1, 0, 1), repeat=4):
        if entries[0] * entries[3] == entries[1] * entries[2]:
            continue
        u = gm.matrix((entries[:2], entries[2:]))
        forms = gm.invariant_forms(2, (u,))
        for a, b, c in product((-1, 0, 1), repeat=3):
            g = gm.matrix(((a, b), (b, c)))
            require((gm.congruence(u, g) == g) == in_span(g, forms),
                    'U17 direct candidate-form oracle mismatch')
            counts['candidate_form_comparisons'] += 1
        require(all(gm.congruence(u, g) == g for g in forms), 'U17 basis contains invalid form')
        counts['integer_single_transports'] += 1

    # Reconstruct node forms from root form and tree transports, then check all edges.
    n, h = 3, lorentz(3)
    b1, b2 = gm.emk_step(3, 1, Q(1, 2)), gm.emk_step(3, 2, Q(1, 3))
    paths = (gm.identity(n), gm.matrix(((1, 1, 0), (0, 1, 0), (0, 0, 2))),
             gm.matrix(((2, 0, 0), (0, 1, 1), (0, 0, 1))))
    metrics = tuple(gm.metric_at_root(h, path) for path in paths)
    for x, y, loop in product(range(3), range(3), (gm.identity(n), b1, b2)):
        u = gm.mul(gm.mul(paths[y], loop), gm.inverse(paths[x]))
        require(gm.congruence(u, metrics[y]) == metrics[x], 'U17 graph compatibility')
        counts['graph_edge_checks'] += 1
    dilation = gm.scale(gm.identity(2), 2)
    require(not gm.invariant_forms(2, (dilation,)), 'U17 nonmetric loop accepted')

    for n in range(2, 6):
        cut = tuple(tuple(int(i == j) for j in range(n)) for i in range(n - 1))
        require(gm.descend_form(lorentz(n), cut) is None, 'U19 lossy nondegenerate descent')
        target = gm.identity(n - 1)
        degenerate = gm.congruence(cut, target)
        require(gm.descend_form(degenerate, cut) == target, 'U19 radical quotient rejected')
        counts['cut_form_checks'] += 2

    # Positive lossy-cut bridge: the discarded memory is exactly the form radical.
    b = gm.emk_step(2, 1, Q(1, 2))
    lift = gm.matrix((b[0] + (0,), b[1] + (0,), (1, 0, 2)))
    memory_step = gm.matrix(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
    source_form = gm.matrix(((1, 0, 0), (0, -1, 0), (0, 0, 0)))
    source_forms = gm.invariant_forms(3, (lift, memory_step))
    memory_cut = gm.matrix(((1, 0, 0), (0, 1, 0)))
    require(len(source_forms) == 1 and in_span(source_form, source_forms),
            'W16 native family did not determine the degenerate source target')
    require(gm.descend_form(source_form, memory_cut) == lorentz(2) and
            gm.descend_transport(lift, memory_cut) == b and
            gm.descend_transport(memory_step, memory_cut) == gm.identity(2),
            'W16 lawful lossy metric/transport bridge failed')
    counts['lawful_memory_cut_examples'] += 1

    cycle = gm.matrix(((0, 0, 1), (1, 0, 0), (0, 1, 0)))
    back = gm.mul(cycle, cycle)
    cut = gm.matrix(((1, 0, 0), (0, 1, 0)))
    section = gm.transpose(cut)
    ca, cb = gm.mul(gm.mul(cut, cycle), section), gm.mul(gm.mul(cut, back), section)
    require(not any(flat(gm.commutator(cycle, back))), 'W13 native pair should commute')
    require(any(flat(gm.commutator(ca, cb))), 'W13 compression should invent a commutator')
    require(gm.descend_transport(cycle, cut) is None and gm.descend_transport(back, cut) is None,
            'W13 illegal cut passed transport descent')

    for coefficients in ((1, 0, 0), (3, 1, 0), (1, 0, 1), (2, -1, 2)):
        for u in (Q(-2), Q(-1), Q(0), Q(1, 2), Q(1), Q(2)):
            warped_chart(u, coefficients)
            counts['warped_chart_points'] += 1
    curved = warped_chart(Q(0), (1, 0, 1))
    zero = gm.matrix(((0, 0), (0, 0)))
    k = gm.generator(2, 1)
    false_torsion = gm.torsion_2d(gm.identity(2), zero, zero, zero, k)
    require(false_torsion == (0, -1), 'W14 torsion gate failed')
    four = tuple(gm.emk_step(4, i, Q(1, 2)) for i in range(1, 4))
    require(len(gm.invariant_forms(4, four[:1])) == 4, 'W12 missing-sector ambiguity')
    loop = gm.mul(gm.mul(gm.mul(b1, b2), gm.inverse(b1)), gm.inverse(b2))
    require(loop != gm.identity(3) and gm.congruence(loop, h) == h, 'W10 metric/nonflat compatibility')
    return rational_records({'checks': counts,
            'rational_encoding': 'Fraction values are exact integer or p/q strings',
            'controls': {
                'W10_EMK_candidate_signature': {'dimension': 4, 'form_up_to_scale': lorentz(4),
                                               'nontrivial_three_component_loop': loop},
                'W11_sector_nonselection': {'K_form': lorentz(4), 'R_form': gm.identity(4),
                                           'combined_invariant_dimension': 0},
                'W12_incomplete_family': {'dimension': 4, 'one_K_pair_invariant_dimension': 4},
                'W13_compression_artifact': {'native_commutator': gm.commutator(cycle, back),
                                             'compressed_commutator': gm.commutator(ca, cb),
                                             'exact_descent': False},
                'W14_torsion_blocks_identification': {'coframe': gm.identity(2), 'at_u': 1,
                                                      'fibre_curvature': k, 'metric_curvature': zero,
                                                      'torsion': false_torsion},
                'W15_valid_curved_adapter_at_zero': curved,
                'W16_lawful_memory_cut': {'source_form': source_form, 'cut': memory_cut,
                                         'native_lift': lift, 'memory_step': memory_step,
                                         'quotient_form': lorentz(2), 'quotient_transport': b}},
            'scope': 'Conditional linear transport and supplied smooth-chart checks. No physical gravity law, continuum emergence, selected dimension, quantum readout rule or empirical validation.'})
