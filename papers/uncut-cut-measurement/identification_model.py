"""Identification of a declared native control pair; no gravity identification.

The protocol has three reported statistics from five raw scalar readings.
Uncertainty boxes give conditional enclosures, not a proof that a model exists.
"""
from fractions import Fraction as Q

import gravity_model as gm
import grading_model as gr
import interaction_model as im


def pair(phase, mixing):
    return im.mixing_step(2, 0, 1, mixing), gr.edge_step(2, 0, phase)


def protocol(u, v, frame=None):
    """Predeclared F maps canonical preparations to the probe carrier.

    Readouts use F^-1. F must be independently calibrated, never fitted here.
    Raw order: first/second coordinates after V,U and U,V on e3; Ue1 first.
    """
    u, v = gm.matrix(u), gm.matrix(v)
    if (len(u), len(u[0]), len(v), len(v[0])) != (4, 4, 4, 4):
        raise ValueError('This protocol requires the declared four-real-coordinate carrier')
    if frame is not None:
        f = gm.matrix(frame)
        if (len(f), len(f[0])) != (4, 4):
            raise ValueError('Calibration frame must be invertible and four dimensional')
        fi = gm.inverse(f)
        u, v = gm.mul(gm.mul(fi, u), f), gm.mul(gm.mul(fi, v), f)
    uv, vu = gm.mul(u, v), gm.mul(v, u)
    raw = uv[0][2], vu[0][2], uv[1][2], vu[1][2], u[0][0]
    return {'raw': raw, 'reported': (raw[0]-raw[1], raw[2]-raw[3], raw[4])}


def forward_formula(phase, mixing):
    t, s = Q(phase), Q(mixing)
    if not t or not s:
        raise ValueError('Both finite Cayley parameters must be nonzero')
    k = 2*s/(1+s*s)
    return -2*k*t*t/(1+t*t), 2*k*t/(1+t*t), (1-s*s)/(1+s*s)


def order_invariants(d, e):
    d, e = Q(d), Q(e)
    if not d or not e:
        raise ValueError('Nonzero parameters require nonzero signed components d and e')
    t, k = -d/e, -(d*d+e*e)/(2*d)
    if abs(k) > 1:
        raise ValueError('Response is outside the declared native pair family')
    return {'phase': t, 'mixing_sine': k}


def decode_exact(record):
    values = tuple(Q(x) for x in record)
    if len(values) != 3:
        raise ValueError('Supply the three calibrated reported statistics')
    d, e, c = values
    found = order_invariants(d, e)
    k = found['mixing_sine']
    if not -1 < c < 1 or c*c+k*k != 1:
        raise ValueError('Direct response fails the native circle consistency gate')
    s = k/(1+c)
    if forward_formula(found['phase'], s) != values:
        raise RuntimeError('Exact inverse failed its forward identity')
    return {'phase': found['phase'], 'mixing': s, 'mixing_sine': k}


def compare_exact(first, second):
    """Agreement only within the declared two-control family and calibration."""
    try:
        a, b = decode_exact(first), decode_exact(second)
    except ValueError:
        return 'MODEL_CLASS_REJECTED'
    return ('EXACT_SHARED_PAIR' if (a['phase'], a['mixing']) == (b['phase'], b['mixing'])
            else 'SHARED_PAIR_REJECTED')


def interval(center, radius=0):
    center, radius = Q(center), Q(radius)
    if radius < 0:
        raise ValueError('Error bounds must be nonnegative and independently justified')
    return center-radius, center+radius


def iadd(a, b):
    return a[0]+b[0], a[1]+b[1]


def ineg(a):
    return -a[1], -a[0]


def imul(a, b):
    p = [x*y for x in a for y in b]
    return min(p), max(p)


def isquare(a):
    return (Q(0) if a[0] <= 0 <= a[1] else min(x*x for x in a), max(x*x for x in a))


def idiv(a, b):
    if b[0] <= 0 <= b[1]:
        raise ValueError('Division interval contains zero')
    return imul(a, (Q(1)/b[1], Q(1)/b[0]))


def contains(a, x):
    return a[0] <= x <= a[1]


def disjoint(a, b):
    return a[1] < b[0] or b[1] < a[0]


def boxes(record, errors):
    record, errors = tuple(record), tuple(errors)
    if len(record) != 3 or len(errors) != 3:
        raise ValueError('Three reported values and three error radii are required')
    return tuple(interval(x, r) for x, r in zip(record, errors))


def decode_box(record, errors):
    d, e, c = boxes(record, errors)
    result = {'measurement_boxes': (d, e, c)}
    if d == (0, 0) or e == (0, 0) or c[0] >= 1 or c[1] <= -1:
        return dict(result, status='MODEL_CLASS_REJECTED')
    if contains(d, 0) or contains(e, 0) or contains(iadd((Q(1), Q(1)), c), 0):
        return dict(result, status='INSUFFICIENT_RESOLUTION')
    t = ineg(idiv(d, e))
    k = ineg(idiv(iadd(isquare(d), isquare(e)), imul((Q(2), Q(2)), d)))
    circle = iadd(iadd(isquare(k), isquare(c)), (Q(-1), Q(-1)))
    if not contains(circle, 0):
        return dict(result, status='MODEL_CLASS_REJECTED', circle_residual=circle)
    s = idiv(k, iadd((Q(1), Q(1)), c))
    exact = all(a == b for a, b in (d, e, c))
    return dict(result, status='EXACT_MODEL_MATCH' if exact else 'CONDITIONAL_ENCLOSURE',
                phase=t, mixing=s, mixing_sine=k, circle_residual=circle)


def compare_boxes(first, first_errors, second, second_errors):
    a, b = decode_box(first, first_errors), decode_box(second, second_errors)
    if 'MODEL_CLASS_REJECTED' in (a['status'], b['status']):
        return 'MODEL_CLASS_REJECTED'
    if any(disjoint(x, y) for x, y in zip(a['measurement_boxes'], b['measurement_boxes'])):
        return 'SHARED_PAIR_REJECTED'
    if 'phase' in a and 'phase' in b:
        if disjoint(a['phase'], b['phase']) or disjoint(a['mixing'], b['mixing']):
            return 'SHARED_PAIR_REJECTED'
        if a['status'] == b['status'] == 'EXACT_MODEL_MATCH':
            return 'EXACT_SHARED_PAIR'
    return 'UNRESOLVED'


def candidate_compatible(phase, mixing, record, errors):
    return all(contains(a, x) for a, x in zip(boxes(record, errors), forward_formula(phase, mixing)))


def ratio_error_bound(numerator, denominator, numerator_error, denominator_error):
    """Bound on a/b relative to measured numerator/denominator, if separated."""
    a, b, ea, eb = map(Q, (numerator, denominator, numerator_error, denominator_error))
    if ea < 0 or eb < 0 or abs(b) <= eb:
        raise ValueError('A positive denominator margin is required')
    return (ea*abs(b)+abs(a)*eb)/(abs(b)*(abs(b)-eb))
