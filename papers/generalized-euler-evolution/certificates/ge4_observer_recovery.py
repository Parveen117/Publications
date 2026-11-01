"""GE4 native observer-memory recovery. Python 3.12 only.

Exact finite controls support GE4_MINIMAL_OBSERVER_MEMORY.md; they do not
certify an unknown physical source or the interacting row. --check is
read-only; --write changes only GE4 evidence.
"""
import argparse
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import sys

import ge3_returning_memory as prior

a, q, y, p = prior.a, prior.q, prior.y, prior.p
old, native, Iv, refine = prior.old, prior.native, prior.Iv, prior.refine
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PINS = HERE/'GE4_SOURCE_PINS.json'
RESULT = HERE/'GE4_RESULT.json'
EXPECTED = HERE/'EXPECTED_GE4.sha256'


def moments(values):
    values = tuple(q.exact(x) for x in values)
    if not values or values[0] != 1:
        raise ValueError('normalized exact moments with m0=1 required')
    return values


def hankel(values, n, shift=0):
    values = moments(values); old.integer(n, 1); old.integer(shift)
    if len(values) < 2*n-1+shift:
        raise ValueError('insufficient moments for this word Gram matrix')
    return [[values[i+j+shift] for j in range(n)] for i in range(n)]


def psd_rank(H):
    if not a.psd(H):
        raise ValueError('moment Gram is not positive semidefinite')
    return a.y37.inertia(H)[0]


def model_moments(L, G, count):
    old.integer(count)
    L, G = old.pairing(L, G)
    if a.multiply(G, L) != a.transpose(a.multiply(G, L)):
        raise ValueError('self-adjoint model required')
    v = [F(i == 0) for i in range(len(L))]
    out = []
    for _ in range(count+1):
        out.append(a.dot(G[0], v))
        v = a.apply(L, v)
    return tuple(out)


def realize(values, dimension, ceiling):
    values = moments(values); old.integer(dimension, 1)
    ceiling = q.exact(ceiling)
    if ceiling <= 0:
        raise ValueError('positive declared generator ceiling required')
    G = hankel(values, dimension)
    extended = hankel(values, dimension+1)
    if psd_rank(G) != dimension or psd_rank(extended) != dimension:
        raise ValueError('exact flat-rank certificate absent at this dimension')
    N = hankel(values, dimension, 1)
    if not a.psd(N) or not a.psd(a.add(a.scale(G, ceiling), N, -1)):
        raise ValueError('shifted moments violate the positive generator interval')
    L = a.multiply(a.inverse(G), N)
    if model_moments(L, G, len(values)-1) != values:
        raise ValueError('later supplied moments contradict flat reconstruction')
    return {'G': G, 'L': L, 'dimension': dimension, 'memory_dimension': dimension-1}


def three_jet(values, ceiling):
    values = moments(values); ceiling = q.exact(ceiling)
    if len(values) < 4 or ceiling <= 0:
        raise ValueError('three derivatives and positive ceiling required')
    A, m2, m3 = values[1:4]
    variance = m2-A*A
    if not 0 <= A <= ceiling or variance < 0 or m2 > ceiling*A:
        raise ValueError('moments violate positive bounded-source constraints')
    if variance == 0:
        if m3 != A**3:
            raise ValueError('zero variance contradicts the supplied third moment')
        return {'a': A, 'q': variance, 'ceiling_exact': True, 'G': [[F(1)]],
                'L': [[A]], 'memory_dimension': 0, 'Xi': F(0)}
    if A == ceiling:
        raise ValueError('positive variance at the ceiling is impossible')
    r = (ceiling*A-m2)/(ceiling-A)
    Xi = -m3+(ceiling+2*r)*m2-(2*ceiling*r+r*r)*A+ceiling*r*r
    if not 0 <= r < A or Xi < 0:
        raise ValueError('ceiling localizer is inadmissible')
    delta = (m3-2*A*m2+A**3)/variance
    G = a.diag((1, variance))
    L = [[A, variance], [F(1), delta]]
    GL = a.multiply(G, L)
    if not a.psd(GL) or not a.psd(a.add(a.scale(G, ceiling), GL, -1)):
        raise ValueError('two-word compression violates the source interval')
    return {'a': A, 'q': variance, 'r': r, 'Xi': Xi, 'delta': delta,
            'ceiling_exact': Xi == 0, 'G': G, 'L': L,
            'squared_error_coefficient': ceiling*Xi/variance}


def native_data(C, poly, count=8):
    norm = y.phi(p.mul(poly, poly))
    if norm <= 0:
        raise ValueError('nonzero native probe required')
    out, word = [], poly
    for _ in range(count+1):
        out.append(y.phi(p.mul(poly, word))/norm)
        word = q.lap(word, C)
    return tuple(out)


def diagonal_words():
    return tuple(p.add(p.mul(p.COORD[0], p.COORD[0]), p.mul(p.COORD[j], p.COORD[j]),
                       p.scale(y.RADIUS, F(-1, 2))) for j in (1, 2, 3))


def scalar_model_enclosure(model, time, ceiling=F(1), order=36):
    time, ceiling = q.exact(time), q.exact(ceiling)
    old.integer(order, 1)
    if time < 0 or ceiling <= 0:
        raise ValueError('nonnegative clock and positive ceiling required')
    mm = model_moments(model['L'], model['G'], order)
    center = sum(((-time)**k*mm[k]/factorial(k) for k in range(order+1)), F(0))
    radius = refine.factorial_tail(time*ceiling, order)
    return Iv(center-radius, center+radius)


def native_controls():
    C0 = a.diag((0, F(1, 2), F(1, 2)))
    _, _, probe, _ = prior.sector()
    m = native_data(C0, probe)
    assert m[:5] == (F(1), F(3, 4), F(5, 8), F(9, 16), F(17, 32))
    recovered = three_jet(m[:4], 1)
    assert recovered['ceiling_exact'] and recovered['q'] == F(1, 16)
    assert recovered['r'] == F(1, 2) and recovered['delta'] == F(3, 4)
    full = realize(m, 2, 1)
    assert full['memory_dimension'] == 1
    assert model_moments(recovered['L'], recovered['G'], 8) == m
    words = diagonal_words(); C1 = a.diag((F(1, 6), F(1, 3), F(1, 2)))
    for j, word in enumerate(words):
        assert q.lap(word, C1) == p.scale(word, 1-C1[j][j])
        for k, other in enumerate(words):
            assert y.phi(p.mul(word, other)) == (F(1, 12) if j == k else 0)
    combined = p.add(*words)
    m3 = native_data(C1, combined)
    rates = (F(5, 6), F(2, 3), F(1, 2))
    assert m3 == tuple(sum((x**k for x in rates), F(0))/3 for k in range(9))
    exact = realize(m3, 3, 1)
    approx = three_jet(m3[:4], 1)
    assert exact['memory_dimension'] == 2
    assert approx['q'] == F(1, 54) and approx['Xi'] == F(5, 972)
    assert approx['r'] == F(11, 18) and not approx['ceiling_exact']
    # Independent leakage check in the three actual eigenword coordinates.
    J = [[F(1), rate-approx['a']] for rate in rates]
    residual = a.add(a.multiply(a.diag(rates), J), a.multiply(J, approx['L']), -1)
    source_G = a.scale(a.identity(3), F(1, 3))
    gram_residual = a.multiply(a.multiply(a.transpose(residual), source_G), residual)
    assert a.psd(a.add(a.scale(approx['G'], approx['squared_error_coefficient']), gram_residual, -1))
    rows = []
    for t in (F(1, 4), F(1), F(2)):
        target = sum((refine.ex(-rate*t) for rate in rates), Iv(F(0)))/3
        enclosure = scalar_model_enclosure(approx, t)
        error = old.abs_ceiling(enclosure-target)
        bound2 = t*t*approx['squared_error_coefficient']
        assert error*error <= bound2
        rows.append({'tau': str(t), 'error_upper': str(error), 'bound_squared': str(bound2)})
    return {'GE3_three_jet_recovery': {k: str(recovered[k]) for k in ('a', 'q', 'r', 'Xi', 'delta')},
            'native_minimum_memory_dimensions': [1, 2],
            'three_mode_moments': [str(x) for x in m3[:7]],
            'native_flat_Gram_ranks': [[psd_rank(hankel(m, n)) for n in (2, 3)],
                                       [psd_rank(hankel(m3, n)) for n in (3, 4)]],
            'two_state_approximation_enclosures': rows}


def boundary_controls():
    L = [[F(3, 4), F(1, 4), F(0)], [F(1, 4), F(3, 4), F(1, 8)],
         [F(0), F(1, 8), F(3, 4)]]
    I = a.identity(3)
    assert a.psd(L) and not a.psd(a.add(I, L, -1))
    m = model_moments(L, I, 8)
    assert m[:4] == (F(1), F(3, 4), F(5, 8), F(9, 16))
    assert m[4] != F(17, 32) and realize(m, 3, 2)['memory_dimension'] == 2
    # The ceiling is an independently warranted input, not inferable from 3 jets.
    assert three_jet(m[:4], 1)['ceiling_exact']
    assert not three_jet(m[:4], 2)['ceiling_exact']
    # Xi>0 does not imply that more than one memory coordinate is needed.
    below = tuple((F(1, 4)**k+F(3, 4)**k)/2 for k in range(9))
    assert three_jet(below[:4], 1)['Xi'] > 0
    assert realize(below, 2, 1)['memory_dimension'] == 1
    base = prior.L
    base_m = model_moments(base, a.identity(2), 8)
    invisible = []
    for eps in (F(1, 8), F(1, 100), F(1, 10000)):
        extended = [base[0]+[F(0)], base[1]+[F(0)], [F(0), F(0), eps]]
        assert a.psd(extended) and a.psd(a.add(I, extended, -1))
        assert model_moments(extended, I, 8) == base_m
        assert realize(base_m, 2, 1)['memory_dimension'] == 1
        invisible.append(str(eps))
    return {'missing_ceiling_counterexample_fourth_moment': str(m[4]),
            'same_first_three_moments': [str(x) for x in m[:4]],
            'positive_localizer_can_still_have_one_memory_coordinate': True,
            'invisible_appended_rates': invisible,
            'all_moment_blindness': 'proved by direct-sum identity; finite checks through order eight'}


def refusals():
    good = tuple((F(1, 2)**k+F(1)**k)/2 for k in range(9))
    tri = tuple(sum((x**k for x in (F(1, 2), F(2, 3), F(5, 6))), F(0))/3 for k in range(9))
    cases = {'float_moment': lambda: three_jet((1, 0.75, F(5, 8), F(9, 16)), 1),
             'unnormalized_probe': lambda: three_jet((2, 1, 1, 1), 1),
             'insufficient_derivatives': lambda: three_jet((1, 1, 1), 1),
             'zero_ceiling': lambda: three_jet(good[:4], 0),
             'negative_variance': lambda: three_jet((1, F(3, 4), F(1, 4), F(1, 8)), 1),
             'inconsistent_zero_variance': lambda: three_jet((1, F(1, 2), F(1, 4), F(1, 4)), 1),
             'false_flat_rank': lambda: realize(tri, 2, 1),
             'later_moment_tamper': lambda: realize(good[:-1]+(good[-1]+1,), 2, 1),
             'unsupported_ceiling': lambda: realize(good, 2, F(3, 4))}
    out = {}
    for key, call in cases.items():
        try:
            call()
        except ValueError:
            out[key] = 'rejected'
        else:
            raise AssertionError('invalid data accepted: '+key)
    return out


def source_checks(root=ROOT, pins=None):
    pins = json.loads(PINS.read_text()) if pins is None else pins
    checked = prior.source_checks(root, pins)
    if root == ROOT:
        checked[str(PINS.relative_to(ROOT))] = old.digest(PINS)
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or not __debug__:
        raise RuntimeError('Python 3.12 without optimization is required')
    return {'certificate_type': 'GE4_MINIMAL_NATIVE_OBSERVER_MEMORY', 'verdict': 'PASS',
            'input_sha256': source_checks(), 'native_recovery': native_controls(),
            'identification_boundaries': boundary_controls(), 'refusals': refusals(),
            'scope': {'written_results': ['GE4-T1', 'GE4-T2', 'GE4-T3', 'GE4-T4'],
                      'source': 'one declared finite positive self-adjoint core and calibrated native clock',
                      'ceiling': 'must be independently warranted; not inferred from three derivatives',
                      'minimum': 'observable linear cyclic realization, not all ambient hidden states',
                      'certification': 'scoped written proofs and exact finite checks, not formal or expert certification',
                      'open': ['noisy experimental derivative inference', 'actual interacting row',
                               'moving and unbounded cuts', 'physical selection', '4D Yang-Mills', 'Clay']}}


def check(cert, result=RESULT, expected=EXPECTED):
    sha = y.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or expected.read_text().strip() != sha:
        raise ValueError('GE4 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true'); group.add_argument('--check', action='store_true')
    args = parser.parse_args(); cert = run(); sha = y.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n'); EXPECTED.write_text(sha+'\n')
    else:
        check(cert)
    print('GE4 PASS', sha)
    print('4 written results; native minimum memory counts 1 and 2; ambient gap and Clay OPEN')


if __name__ == '__main__':
    main()
