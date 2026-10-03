"""Independent finite comparisons for U30-U32; exact rational arithmetic."""
from fractions import Fraction as Q
from itertools import combinations, product

import finite_model as fm
import gravity_model as gm
import grading_model as gr
import interaction_model as im
import response_model as rm
from grading_certificate import apply
from gravity_certificate import rational_records


def require(value, message):
    if not value:
        raise RuntimeError(message)


def paired_gram_oracle(y):
    """Direct Gaussian-rational inner products, independent of Z multiplication."""
    cols = tuple(zip(*y))
    real, imag = [], []
    for u in cols:
        real.append(tuple(sum(u[2*a]*v[2*a]+u[2*a+1]*v[2*a+1]
                              for a in range(len(u)//2)) for v in cols))
        imag.append(tuple(sum(u[2*a]*v[2*a+1]-u[2*a+1]*v[2*a]
                              for a in range(len(u)//2)) for v in cols))
    return tuple(real), tuple(imag)


def pairwise_discard(points, weights, labels):
    """Weighted pair-distance formula; does not compute conditional means."""
    d = len(points[0])
    out = rm.zero(d)
    for i, j in combinations(range(len(points)), 2):
        if labels[i] != labels[j]:
            continue
        mass = sum(weights[k] for k in range(len(points)) if labels[k] == labels[i])
        delta = tuple(points[i][a]-points[j][a] for a in range(d))
        outer = gm.matrix([[a*b for b in delta] for a in delta])
        out = gm.add(out, gm.scale(outer, weights[i]*weights[j]/mass))
    return out


def loop_datum(m, a, b, phase, mixing):
    u = im.mixing_step(m, a, b, mixing)
    v = gr.edge_step(m, a, phase)
    return u, v, rm.order_response(u, v)


def controls():
    z = rm.quarter_turn(1)
    split = rm.split_tensor(gm.identity(2), ((1, 0),), z)
    gap = rm.tensor_sub(split['full'], split['visible'])
    determinant = gap[0][0][0]*gap[0][1][1]-gap[0][0][1]**2-gap[1][0][1]**2
    require(determinant == -1 and split['cross'][1] != rm.zero(2), 'W33 false complex monotonicity')
    require(rm.paired_repair(((1, 0),), z)['extra_scalar_channels'] == 1, 'W33 repair')
    u, v, response = loop_datum(2, 0, 1, Q(1, 2), Q(1, 2))
    c = im.projection(2, (0,))
    require(rm.split_tensor(gm.identity(4), c, rm.quarter_turn(2))['global_ledger_closes'], 'W34 paired cut')
    require(gm.descend_transport(u, c) is None, 'W34 tensor closure confused with dynamics')
    signed = apply(response, (0, 0, 0, 1))
    require(signed == (Q(-16, 25), Q(-8, 25), 0, 0), 'W35 signed order')
    require(gm.mul(gm.transpose(response), response) == gm.scale(gm.identity(4), Q(64, 125)), 'W35 response scale')
    require(rm.tensor(response, rm.quarter_turn(2)) == rm.tensor(gm.scale(response, -1), rm.quarter_turn(2)), 'W37 sign blindness')
    ensemble = (signed, tuple(-a for a in signed))
    ledger = rm.conditional_ledger(ensemble, (Q(1, 2), Q(1, 2)), (0, 0))
    require(ledger['mean'] == (0, 0, 0, 0) and ledger['discarded'] != rm.zero(4), 'W37 mean blindness')
    u3, v3, r3 = loop_datum(3, 0, 1, Q(1, 2), Q(1, 2))
    c3 = im.projection(3, (0, 1))
    ubar, vbar = gm.descend_transport(u3, c3), gm.descend_transport(v3, c3)
    require(gm.mul(c3, r3) == gm.mul(rm.order_response(ubar, vbar), c3), 'W36 quotient response')
    require(gr.rank(r3) == 4, 'W36 proper lost component')
    return {'W33_real_cut_complex_seam': {'gap_determinant': determinant,
                'cross_imaginary': split['cross'][1], 'minimum_extra': 1},
            'W34_paired_cut_not_dynamically_closed': {'tensor_ledger': True, 'mixing_descends': False},
            'W35_signed_response': {'response': signed, 'squared_size': Q(64, 125),
                'full_target_rank': 4, 'norm_target_also_requires_rank': 4},
            'W36_proper_lossy_response_quotient': {'source_dimension': 6, 'target_dimension': 4},
            'W37_zero_average_and_sign_blindness': {'mean': ledger['mean'],
                'discarded_covariance': ledger['discarded'], 'opposite_responses_same_tensor': True}}


def build_checks():
    counts = {'tensor_transport_checks': 0, 'coordinate_cuts': 0,
              'paired_coordinate_cuts': 0, 'general_readout_repairs': 0,
              'loop_scale_and_quotient_checks': 0, 'statistical_ledgers': 0,
              'nested_statistical_ledgers': 0}
    # All graphs from U27, but an independent entrywise complex Gram oracle.
    for m in range(1, 5):
        n, z = 2*m, rm.quarter_turn(m)
        y = gm.matrix([[Q((i+1)*(j+2) % 7 - 3, j+1) for j in range(3)] for i in range(n)])
        t = rm.tensor(y, z)
        require(t == paired_gram_oracle(y), 'U30 independent complex Gram mismatch')
        jy = gm.mul(im.native_grading(m), y)
        require(rm.tensor(jy, z) == (t[0], gm.scale(t[1], -1)), 'U30 grading conjugation')
        pairs = tuple(combinations(range(m), 2))
        for flags in product((0, 1), repeat=len(pairs)):
            links = tuple(e for e, yes in zip(pairs, flags) if yes)
            for u in im.family(m, links):
                require(paired_gram_oracle(gm.mul(u, y)) == t, 'U30 transport failed')
                counts['tensor_transport_checks'] += 1
        # Include the zero and identity cuts, independently count complete pairs.
        for flags in product((0, 1), repeat=n):
            cut = tuple(gm.identity(n)[i] for i, yes in enumerate(flags) if yes) or (tuple(Q(0) for _ in range(n)),)
            split = rm.split_tensor(gm.identity(n), cut, z)
            require(split['full'] == rm.tensor_add(rm.tensor_add(split['visible'], split['hidden']), split['cross']), 'U31 ledger failed')
            paired = all(flags[2*a] == flags[2*a+1] for a in range(m))
            require(split['global_ledger_closes'] == paired, 'U31 paired-cut oracle')
            require((split['cross'][1] == rm.zero(n)) == paired, 'U31 iff detected by full catalogue')
            repair = rm.paired_repair(cut, z)
            missing = sum(flags[2*a] != flags[2*a+1] for a in range(m))
            require(repair['extra_scalar_channels'] == missing, 'U31 minimum coordinate repair')
            require(rm.split_tensor(gm.identity(n), repair['observer'], z)['global_ledger_closes'], 'U31 repaired projector')
            counts['coordinate_cuts'] += 1
            counts['paired_coordinate_cuts'] += paired
        for k in range(1, min(n, 4)+1):
            cut = gm.matrix([[Q((i+2)*(j+1) % 5-2, i+1) for j in range(n)] for i in range(k)])
            repair = rm.paired_repair(cut, z)
            if gr.rank(cut):
                iterative = im.row_closure(cut, (z,))
                require(gr.rank(repair['observer']+iterative['observer']) == repair['rank'] == len(iterative['observer']), 'U31 independent row closure')
            require(rm.split_tensor(gm.identity(n), repair['observer'], z)['global_ledger_closes'], 'U31 noncoordinate repair')
            counts['general_readout_repairs'] += 1
    parameters = (Q(-2, 3), Q(1, 2), Q(1))
    for m in range(2, 5):
        for a, b in combinations(range(m), 2):
            c = im.projection(m, (a, b))
            for phase, mixing in product(parameters, repeat=2):
                u, v, ell = loop_datum(m, a, b, phase, mixing)
                kappa = 16*phase**2*mixing**2/((1+phase**2)*(1+mixing**2)**2)
                expected = gm.scale(gm.mul(gm.transpose(c), c), kappa)
                require(gm.mul(gm.transpose(ell), ell) == expected and gr.rank(ell) == 4, 'U32 response scale/rank')
                ub, vb = gm.descend_transport(u, c), gm.descend_transport(v, c)
                require(gm.mul(c, ell) == gm.mul(rm.order_response(ub, vb), c), 'U32 exact response square')
                loop = gm.sub(gm.mul(gm.mul(gm.mul(u, v), gm.transpose(u)), gm.transpose(v)), gm.identity(2*m))
                require(loop == gm.mul(gm.mul(ell, gm.transpose(u)), gm.transpose(v)), 'U32 loop convention')
                counts['loop_scale_and_quotient_checks'] += 1
    weights_list = ((Q(1, 4),)*4, (Q(1, 10), Q(2, 10), Q(3, 10), Q(4, 10)),
                    (Q(1, 2), Q(1, 6), Q(1, 6), Q(1, 6)))
    cuts = tuple(fm.partitions(4))
    inputs = ((1, 0, 0, 0), (-1, 0, 0, 0), (0, 0, 0, 1), (0, 0, 0, 1))
    for phase, mixing in product(parameters, repeat=2):
        _, _, ell = loop_datum(2, 0, 1, phase, mixing)
        points = tuple(apply(ell, x) for x in inputs)
        for weights in weights_list:
            for coarse in cuts:
                ledger = rm.conditional_ledger(points, weights, coarse)
                require(ledger['total'] == gm.add(ledger['recognized'], ledger['discarded']), 'U32 covariance ledger')
                require(ledger['discarded'] == pairwise_discard(points, weights, coarse), 'U32 independent pair-distance oracle')
                faithful = all(coarse[i] != coarse[j] or points[i] == points[j] for i, j in combinations(range(4), 2))
                require((ledger['discarded'] == rm.zero(4)) == faithful, 'U32 support-faithfulness iff')
                counts['statistical_ledgers'] += 1
                for fine in cuts:
                    if any(fine[i] == fine[j] and coarse[i] != coarse[j] for i, j in combinations(range(4), 2)):
                        continue
                    detail = rm.conditional_ledger(points, weights, fine)
                    block_labels = tuple(coarse[fine.index(k)] for k in detail['block_labels'])
                    extra = rm.conditional_ledger(detail['block_means'], detail['block_weights'], block_labels)
                    require(ledger['discarded'] == gm.add(detail['discarded'], extra['discarded']), 'U32 tower ledger')
                    counts['nested_statistical_ledgers'] += 1
    return rational_records({'status': 'PASS_FINITE_CHECKS', 'counts': counts,
        'domains': {'channels': [1, 2, 3, 4], 'graphs': 'all simple graphs through four channels',
            'cuts': 'all real coordinate subsets through eight dimensions, including zero and identity',
            'noncoordinate_readouts': 'specified rational modular-entry fixtures',
            'loop_parameters': parameters, 'ensemble_size': 4, 'ensemble_weights': weights_list,
            'ensemble_partitions': len(cuts), 'nested_partitions': 'all refinement-related pairs'},
        'controls': controls(), 'scope': 'U30-U32 written proofs are separate from these finite comparisons; no physical gravity or quantum readout law'} )
