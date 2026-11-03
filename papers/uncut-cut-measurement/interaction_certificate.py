"""Finite exact evidence for a declared native-compatible interacting family."""
from fractions import Fraction as Q
from itertools import combinations, product

import gravity_model as gm
import grading_model as gr
import interaction_model as im
from grading_certificate import apply
from gravity_certificate import in_span, rational_records


def require(value, message):
    if not value:
        raise RuntimeError(message)


def controls():
    u = im.mixing_step(2, 0, 1, Q(1, 2))
    v = gr.edge_step(2, 0, Q(1, 2))
    x = (0, 0, 0, 1)
    uv, vu = apply(gm.mul(u, v), x)[0], apply(gm.mul(v, u), x)[0]
    require((uv, vu) == (0, Q(16, 25)), 'W27 ordered response witness failed')
    loop = gm.mul(gm.mul(gm.mul(u, v), gm.transpose(u)), gm.transpose(v))
    require(loop != gm.identity(4), 'W27 actual native loop is trivial')
    b, g = im.mixing_generator(2, 0, 1), gr.edge_generator(2, 0, -1)
    symmetric_link = gm.matrix(((0, 1), (1, 0)))
    r = gm.matrix(((0, -1), (1, 0)))
    require(gm.commutator(b, g) == im.kron(symmetric_link, r),
            'U27 generator commutator sign')

    family = im.family(3, ((0, 1),))
    cut = im.projection(3, (0, 1))
    targets = tuple(gm.descend_transport(a, cut) for a in family)
    require(all(a is not None for a in targets), 'W28 component cut rejected')
    require(gm.descend_transport(im.native_grading(3), cut) == im.native_grading(2),
            'W28 grading does not descend')
    target_form = im.invariant_form_basis(3, ((0, 1),))[0]
    require(gm.descend_form(target_form, cut) == gm.identity(4), 'W28 form target')
    native_comm = gm.commutator(family[-1], family[0])
    target_comm = gm.commutator(targets[-1], targets[0])
    require(target_comm != gm.scale(gm.identity(4), 0) and
            gm.mul(target_comm, cut) == gm.mul(cut, native_comm),
            'W28 noncommutator was lost or manufactured')

    partial = im.projection(2, (0,))
    require(gm.descend_transport(u, partial) is None, 'W29 partial component cut accepted')
    closed = im.row_closure(partial, im.family(2, ((0, 1),)))
    require(len(closed['observer']) == 4 and closed['extra_scalar_channels'] == 2,
            'W29 insufficient repair')
    cancellation = gm.matrix(((1, 0, -1, 0),))
    cancelled_state = (1, 0, 1, 0)
    require(apply(cancellation, cancelled_state) == (0,), 'W30 initial cancellation')
    response = apply(gm.mul(cancellation, v), cancelled_state)
    require(response == (Q(-2, 5),) and
            len(im.row_closure(cancellation, im.family(2, ()))['observer']) == 4,
            'W30 mixed-component readout falsely closes')
    forms = gm.invariant_forms(4, im.family(2, ((0, 1),)))
    require(len(forms) == 1 and in_span(gm.identity(4), forms), 'W31 invariant form line')
    scalar = (gm.identity(4)[0],)
    require(len(im.row_closure(scalar, (u,))['observer']) == 2,
            'W32 mixing alone should leave internal second coordinates unseen')
    require(len(im.row_closure(scalar, im.family(2, ()))['observer']) == 2,
            'W32 local phase alone should leave the other channel unseen')
    return {
        'W27_actual_order_response': {'initial_state': x, 'R_then_mixing': uv,
            'mixing_then_R': vu, 'difference': uv-vu, 'group_commutator': loop},
        'W28_lossy_interacting_quotient': {'source_dimension': 6, 'target_dimension': 4,
            'lost_component': [2], 'grading_and_form_descend': True,
            'nonzero_target_commutator': target_comm},
        'W29_cut_inside_connected_component': {'exact_descent': False, 'minimum_extra': 2},
        'W30_cancellation_is_not_blind_component': {'initial': 0, 'after_local_R': response,
            'closed_observer_rank': 4},
        'W31_connected_invariant_form': {'space_dimension': 1,
            'nonzero_forms_definite': True, 'Lorentz_signature_derived': False},
        'W32_missing_control': {'mixing_only_rank': 2, 'local_phases_only_rank': 2,
            'full_connected_family_rank': 4}}


def build_checks():
    counts = {'interaction_graphs': 0, 'connected_graphs': 0,
              'graded_orthogonal_steps': 0, 'quarter_turn_commutations': 0,
              'readout_closure_comparisons': 0,
              'component_projection_intertwiners': 0, 'optimal_transcripts': 0,
              'decoder_basis_reconstructions': 0, 'chronological_word_checks': 0,
              'full_invariant_form_solves': 0}
    dimensions = []
    for m in range(1, 5):
        pairs = tuple(combinations(range(m), 2))
        graphs = connected = 0
        for flags in product((False, True), repeat=len(pairs)):
            links = tuple(e for e, included in zip(pairs, flags) if included)
            steps = im.family(m, links)
            comp = im.components(m, links)
            unit, j = gm.identity(2*m), im.native_grading(m)
            quarter_turn = im.kron(gm.identity(m), ((0, -1), (1, 0)))
            graphs += 1
            counts['interaction_graphs'] += 1
            if len(comp) == 1:
                connected += 1
                counts['connected_graphs'] += 1
            for index, u in enumerate(steps):
                require(gm.mul(gm.transpose(u), u) == unit, 'U27 energy preservation')
                require(gm.mul(gm.mul(j, u), j) == (gm.transpose(u) if index < m else u),
                        'U27 common native grading law')
                counts['graded_orthogonal_steps'] += 1
                require(gm.mul(u, quarter_turn) == gm.mul(quarter_turn, u),
                        'U27 represented quarter-turn not preserved')
                counts['quarter_turn_commutations'] += 1
            # First-coordinate readout of each channel and one aggregate readout.
            cuts = tuple((unit[2*a],) for a in range(m))+(gm.matrix(((1,)*(2*m),)),)
            for cut in cuts:
                actual = im.row_closure(cut, steps)
                predicted = im.predicted_observer(m, links, cut)
                a, p = actual['observer'], predicted['observer']
                require(len(a) == len(p) and gr.rank(a+p) == len(a),
                        'U28 component formula disagrees with matrix closure')
                require(actual['extra_scalar_channels'] == predicted['extra_scalar_channels'],
                        'U28 memory count')
                counts['readout_closure_comparisons'] += 1
            for component in comp:
                cut = im.projection(m, component)
                for u in steps:
                    target = gm.descend_transport(u, cut)
                    require(target is not None and gm.mul(target, cut) == gm.mul(cut, u),
                            'U28 exact component projection rejected')
                    counts['component_projection_intertwiners'] += 1
            if m <= 3:
                full_basis = gm.invariant_forms(2*m, steps)
                expected = im.invariant_form_basis(m, links)
                require(len(full_basis) == len(expected) and
                        all(in_span(h, full_basis) for h in expected),
                        'U28 full symmetric form classification failed')
                counts['full_invariant_form_solves'] += 1
            if len(comp) == 1:
                named = {f'R{a}': steps[a] for a in range(m)}
                named.update({f'B{a},{b}': steps[m+k] for k, (a, b) in enumerate(links)})
                for root in range(m):
                    packet = im.optimal_transcript(m, links, root)
                    rows = packet['rows']
                    require(len(rows) == 2*m and gr.rank(rows) == 2*m,
                            'U29 transcript is not minimal full-rank')
                    require(max(map(len, packet['chronological_words'])) <= m,
                            'U29 word-length bound')
                    for column in range(2*m):
                        readings = tuple(row[column] for row in rows)
                        require(im.decode_transcript(packet, readings) == unit[column],
                                'U29 triangular decoder failed a basis vector')
                        counts['decoder_basis_reconstructions'] += 1
                    probe = tuple(Q((-1)**a*(a+1)) for a in range(2*m))
                    for row, word in zip(rows, packet['chronological_words']):
                        state = probe
                        for label in word:
                            state = apply(named[label], state)
                        require(state[2*root] == apply((row,), probe)[0],
                                'U29 chronological word convention mismatch')
                        counts['chronological_word_checks'] += 1
                    counts['optimal_transcripts'] += 1
        dimensions.append({'channels': m, 'graphs': graphs, 'connected': connected})
    return rational_records({'counts': counts, 'by_dimension': dimensions,
        'controls': controls(),
        'coverage': 'All simple interaction graphs on 1-4 two-component channels, phase=mixing=1/2. Every channel first-coordinate readout plus one aggregate row. Every root on connected graphs for minimal transcripts. Full symmetric-form nullspaces on 1-3 channels.',
        'native_rules_used': 'Real anti-self-dagger generators, common J=I tensor K, exact Cayley steps; parameters are protocol labels, not a physical clock.',
        'interaction_graph_and_strengths_supplied': True,
        'physical_gravity_identification': False,
        'quantum_readout_or_instrument_realized': False})
