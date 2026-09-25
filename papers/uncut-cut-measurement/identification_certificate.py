"""Finite exact inverse, calibration and error-box evidence for U33-U35."""
from fractions import Fraction as Q
from itertools import combinations, product

import gravity_model as gm
import grading_model as gr
import interaction_model as im
import identification_model as ident
import response_model as rm
from gravity_certificate import rational_records


def require(value, message):
    if not value:
        raise RuntimeError(message)


def frames():
    return (gm.identity(4), gm.matrix(((2, 0, 0, 0), (0, 3, 0, 0),
                                      (0, 0, 5, 0), (0, 0, 0, 7))),
            gm.matrix(((1, 1, 0, 0), (0, 1, 0, 0), (0, 0, 1, -2), (0, 0, 0, 1))))


def controls():
    t, s = Q(1, 2), Q(1, 2)
    u, v = ident.pair(t, s)
    ur, vr = ident.pair(t, 1/s)
    ell, ellr = rm.order_response(u, v), rm.order_response(ur, vr)
    first, second = ident.protocol(u, v)['reported'], ident.protocol(ur, vr)['reported']
    require(u != ur and ell == ellr and first[:2] == second[:2], 'W39 reciprocal blindness')
    require(first == (Q(-8, 25), Q(16, 25), Q(3, 5)) and second[2] == Q(-3, 5), 'W39 reported values')
    require(rm.tensor(ell, rm.quarter_turn(2)) == rm.tensor(ellr, rm.quarter_turn(2)), 'W39 tensor blindness')
    require(ident.compare_exact(first, second) == 'SHARED_PAIR_REJECTED', 'W40 direct readout ignored')
    up, vp = ident.pair(Q(1, 3), Q(2, 3))
    ellp = rm.order_response(up, vp)
    fitted = gm.mul(ell, gm.inverse(ellp))
    require(ellp != ell and gm.mul(fitted, ellp) == ell, 'W41 fitted-decoder control')
    u6 = im.mixing_step(3, 0, 1, s)
    v6a = gm.mul(gr.edge_step(3, 0, t), gr.edge_step(3, 2, Q(1, 3)))
    v6b = gm.mul(gr.edge_step(3, 0, t), gr.edge_step(3, 2, Q(2, 3)))
    corner = lambda a: tuple(row[:4] for row in a[:4])
    require(ident.protocol(corner(u6), corner(v6a)) == ident.protocol(corner(u6), corner(v6b)), 'W42 hidden extension visible')
    require(v6a[5][4] != v6b[5][4], 'W42 held-out channel failed to differ')
    narrow, broad = (0, 0, Q(1, 10)), (0, 0, Q(3, 5))
    require(ident.compare_boxes(first, narrow, second, narrow) == 'SHARED_PAIR_REJECTED', 'W43 robust separation')
    require(ident.compare_boxes(first, broad, second, broad) == 'UNRESOLVED', 'W43 overlap promoted to equality')
    weak = ident.forward_formula(Q(1, 100), s)
    require(ident.decode_box(weak, (Q(1, 100),)*3)['status'] == 'INSUFFICIENT_RESOLUTION', 'W44 weak-response margin')
    large = ident.forward_formula(t, Q(1000))
    require(ident.decode_box(large, (Q(1, 100000),)*3)['status'] == 'INSUFFICIENT_RESOLUTION', 'W44 mixing denominator')
    # Internal agreement alone cannot select even one numerical pair.
    require(ident.compare_exact(first, first) == ident.compare_exact(ident.forward_formula(Q(1, 3), Q(2, 3)), ident.forward_formula(Q(1, 3), Q(2, 3))) == 'EXACT_SHARED_PAIR', 'NP0 nonselection control')
    return {'W39_reciprocal_blindness': {'phase': t, 'mixings': (s, 1/s),
                'same_complete_order_matrix': True, 'same_tensor': True,
                'reported_first': first, 'reported_second': second},
            'W40_augmented_probe_rejection': 'SHARED_PAIR_REJECTED',
            'W41_posthoc_output_fit': {'different_native_response': True, 'fitted_residual_zero': True,
                'independent_two_sided_calibration': False},
            'W42_hidden_native_extension': {'same_five_raw_readings': True,
                'extra_channel_R_readings': (v6a[5][4], v6b[5][4])},
            'W43_bounded_error': {'direct_gap': Q(6, 5), 'radius_one_tenth': 'SHARED_PAIR_REJECTED',
                'radius_three_fifths': 'UNRESOLVED'},
            'W44_no_uniform_inverse_margin': {'weak_phase': 'INSUFFICIENT_RESOLUTION',
                'large_mixing': 'INSUFFICIENT_RESOLUTION'},
            'NP0_not_physical_selection': {'two_different_pairs_pass_internal_agreement': True,
                'gravity_identified': False}}


def build_checks():
    parameters = tuple(sign*x for sign in (-1, 1) for x in (Q(1, 3), Q(1, 2), Q(1), Q(2), Q(3)))
    counts = {'parameter_pairs': 0, 'held_out_basis_predictions': 0, 'calibrated_frames': 0,
              'pairwise_order_comparisons': 0, 'distinct_reciprocal_collisions': 0,
              'pairwise_augmented_verdicts': 0, 'noisy_enclosures': 0, 'phase_error_bounds': 0}
    packets = []
    eps = Q(1, 10000)
    for t, s in product(parameters, repeat=2):
        u, v = ident.pair(t, s)
        record = ident.protocol(u, v)['reported']
        require(record == ident.forward_formula(t, s), 'U33 matrix protocol/closed formula disagreement')
        recovered = ident.decode_exact(record)
        require((recovered['phase'], recovered['mixing']) == (t, s), 'U33 inverse mismatch')
        prediction = ident.pair(recovered['phase'], recovered['mixing'])
        for actual, predicted in zip((u, v), prediction):
            # e2/e4 preparations are not used by the fitting protocol (e1/e3).
            for column in (1, 3):
                require(tuple(row[column] for row in actual) == tuple(row[column] for row in predicted), 'U34 held-out matrix response')
                counts['held_out_basis_predictions'] += 1
        for f in frames():
            fi = gm.inverse(f)
            up, vp = (gm.mul(gm.mul(f, a), fi) for a in (u, v))
            require(ident.protocol(up, vp, f)['reported'] == record, 'U34 calibration covariance')
            counts['calibrated_frames'] += 1
        ell = rm.order_response(u, v)
        packets.append((t, s, record, ell))
        for offsets in product((-1, 0, 1), repeat=3):
            measured = tuple(x+offset*eps for x, offset in zip(record, offsets))
            decoded = ident.decode_box(measured, (eps,)*3)
            require(decoded['status'] == 'CONDITIONAL_ENCLOSURE', 'U35 compatible box incorrectly rejected')
            require(ident.contains(decoded['phase'], t) and ident.contains(decoded['mixing'], s), 'U35 true parameters outside enclosure')
            require(ident.candidate_compatible(t, s, measured, (eps,)*3), 'U35 known true candidate rejected')
            bound = ident.ratio_error_bound(measured[0], measured[1], eps, eps)
            require(abs(-measured[0]/measured[1]-t) <= bound, 'U35 phase ratio bound')
            counts['noisy_enclosures'] += 1
            counts['phase_error_bounds'] += 1
        counts['parameter_pairs'] += 1
    for a, b in combinations(packets, 2):
        t, s, record, ell = a
        tp, sp, other, ellp = b
        same_order = t == tp and (s == sp or s*sp == 1)
        require((ell == ellp) == same_order, 'U33 complete-fibre classification')
        counts['pairwise_order_comparisons'] += 1
        counts['distinct_reciprocal_collisions'] += same_order
        require(ident.compare_exact(record, other) == 'SHARED_PAIR_REJECTED', 'U34 distinct parameters became same pair')
        counts['pairwise_augmented_verdicts'] += 1
    return rational_records({'status': 'PASS_FINITE_CHECKS', 'counts': counts,
        'domains': {'parameters': parameters, 'frames': frames(), 'reported_error_radius': eps,
            'noise_offsets_per_statistic': [-1, 0, 1], 'reported_statistics': 3,
            'raw_scalar_readings': 5, 'generated_data': 'synthetic exact native matrices'},
        'controls': controls(), 'scope': 'U33-U35 under the declared oriented two-channel family; neither instrument realization nor physical gravitational universality is established'})
