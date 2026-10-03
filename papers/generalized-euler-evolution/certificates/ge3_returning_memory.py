"""GE3 exact native returning-memory controls. Python 3.12 only.

General fixed-core proofs are in GE3_RETURNING_MEMORY.md. These finite
controls do not evaluate the actual interacting infinite-chain row cut.
Default/--check is read-only; --write regenerates GE3 evidence only.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import ge2_clock_free as prior

a, q, y, p = prior.a, prior.q, prior.y, prior.p
old, native, Iv, refine = prior.old, prior.native, prior.Iv, prior.refine
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PINS = HERE/'GE3_SOURCE_PINS.json'
RESULT = HERE/'GE3_RESULT.json'
EXPECTED = HERE/'EXPECTED_GE3.sha256'
L = [[F(3, 4), F(-1, 4)], [F(-1, 4), F(3, 4)]]


def sector():
    e = native.EVEN
    f = p.add(p.mul(p.COORD[0], p.COORD[0]), p.mul(p.COORD[1], p.COORD[1]),
              p.scale(y.RADIUS, F(-1, 2)))
    u, v = p.add(e, f), p.add(e, p.scale(f, -1))
    return e, f, u, v


def word_product(word):
    out = a.identity(2)
    for W in word:
        W = old.square_matrix(W)
        if len(W) != 2:
            raise ValueError('this finite witness adapter requires two modes')
        out = a.multiply(W, out)
    return out


def eliminated_step(word, n, visible, hidden_initial):
    old.integer(n)
    if n >= len(word) or len(visible) < n+1:
        raise ValueError('step or visible history missing')
    W = [old.square_matrix(M) for M in word[:n+1]]
    if any(len(M) != 2 for M in W):
        raise ValueError('two-mode scalar block adapter required')
    hidden_initial = q.exact(hidden_initial)
    visible = [q.exact(x) for x in visible]
    tail = F(1)
    for j in range(n):
        tail *= W[j][1][1]
    out = W[n][0][0]*visible[n]+W[n][0][1]*tail*hidden_initial
    for j in range(n):
        middle = F(1)
        for k in range(j+1, n):
            middle *= W[k][1][1]
        out += W[n][0][1]*middle*W[j][1][0]*visible[j]
    return out


def rational_step(r):
    r = q.exact(r)
    c = prior.turn((r, 0, 0))[0]
    alpha, beta = c*c, 2*c*c-1
    mean, crossing = (alpha+beta)/2, (alpha-beta)/2
    return [[mean, crossing], [crossing, mean]], alpha, beta, r*r/2


def memory_rhs_bound(coupling, hidden_floor, age, source_bound):
    coupling, hidden_floor, age, source_bound = map(q.exact,
        (coupling, hidden_floor, age, source_bound))
    if min(coupling, age, source_bound) < 0 or hidden_floor <= 0:
        raise ValueError('nonnegative data and strictly positive hidden floor required')
    return refine.ex(-hidden_floor*age)*(coupling*coupling*source_bound/hidden_floor)


def native_sector_controls():
    e, f, u, v = sector()
    C = a.diag((0, F(1, 2), F(1, 2)))
    assert q.lap(e, C) == p.scale(e, F(1, 2)) and q.lap(f, C) == f
    assert y.phi(e) == y.phi(f) == 0
    assert y.phi(p.mul(e, e)) == y.phi(p.mul(f, f)) == F(1, 12)
    assert y.phi(p.mul(e, f)) == 0
    embedding = a.transpose([y.vector(u, 2), y.vector(v, 2)])
    assert a.multiply(q.operator(2, C), embedding) == a.multiply(embedding, L)
    phi_gram, _ = native.gram(2)
    pulled = a.multiply(a.multiply(a.transpose(embedding), phi_gram), embedding)
    assert pulled == a.scale(a.identity(2), F(1, 6))
    rows = []
    for r in (F(1, 2), F(1, 3), F(1, 8), F(3)):
        protocol = prior.protocol(2, ((F(1, 2), (0, r, 0)), (F(1, 2), (0, 0, r))))
        W, alpha, beta, h = rational_step(r)
        assert protocol['h'] == h and protocol['C'] == C
        assert a.multiply(protocol['S'], embedding) == a.multiply(embedding, W)
        assert a.psd(a.add(a.identity(2), a.multiply(a.transpose(W), W), -1))
        rows.append({'turn_size': str(r), 'h': str(h), 'e_factor': str(alpha),
                     'f_factor': str(beta), 'return_corner': str(W[0][1])})
    return {'native_eigenvalues': ['1/2', '1'], 'Phi_pairing_in_visible_hidden_basis': 'I/6',
            'native_substitution_controls': rows}


def word_controls():
    word = [rational_step(r)[0] for r in (F(1, 2), F(1, 3), F(1, 4), F(1, 5))]
    checks = 0
    for initial in ((F(1), F(0)), (F(0), F(1)), (F(2), F(-3))):
        states = [list(initial)]
        for W in word:
            states.append(a.apply(W, states[-1]))
        visible = [z[0] for z in states]
        for n in range(len(word)):
            assert eliminated_step(word, n, visible[:n+1], initial[1]) == visible[n+1]
            checks += 1
        assert a.apply(word_product(word), list(initial)) == states[-1]
    # Positive self-adjoint steps can have negative *nonstationary* return terms.
    plus = [[F(3, 5), F(1, 10)], [F(1, 10), F(1, 2)]]
    minus = [[F(3, 5), F(-1, 10)], [F(-1, 10), F(1, 2)]]
    for W in (plus, minus):
        assert a.psd(W) and a.psd(a.add(a.identity(2), W, -1))
    return_term = minus[0][1]*plus[1][0]
    assert return_term == F(-1, 100)
    assert word_product((plus, minus))[0][0]-minus[0][0]*plus[0][0] == return_term
    return {'nonuniform_recursion_checks': checks, 'hidden_initial_sources': 3,
            'nonstationary_positive_step_negative_return': str(return_term)}


def refinement_controls():
    T = F(13, 72)
    retained_target = (refine.ex(-T/2)+refine.ex(-T))/2
    reset_target = refine.ex(-3*T/4)
    assert (retained_target-reset_target).lo > 0
    rows = []
    for n in (1, 2, 4, 8):
        first, second = rational_step(F(1, 2*n)), rational_step(F(1, 3*n))
        count = n*n
        alpha = (first[1]*second[1])**count
        beta = (first[2]*second[2])**count
        retained = (alpha+beta)/2
        reset = (first[0][0][0]*second[0][0][0])**count
        tau = count*(first[3]+second[3])
        assert tau == T
        word = [first[0], second[0]]*count
        assert word_product(word)[0][0] == retained
        bound = prior.constants(2)[1]*F(1, 4*n*n)*T
        er = old.abs_ceiling(Iv(retained)-retained_target)
        ez = old.abs_ceiling(Iv(reset)-reset_target)
        assert er <= bound and ez <= bound and retained > reset
        rows.append({'n': n, 'steps': 2*count, 'tau': str(T),
                     'retained_error_upper': str(er), 'reset_error_upper': str(ez),
                     'bound': str(bound)})
    return {'unequal_step_enclosures': rows, 'different_limit_separation_lower': str((retained_target-reset_target).lo)}


def continuum_controls():
    A, B, D = F(3, 4), F(-1, 4), F(3, 4)
    P = a.diag((1, 0)); Q = a.diag((0, 1))
    marker = a.add(a.multiply(a.multiply(P, a.multiply(L, L)), P),
                   a.multiply(a.multiply(a.multiply(P, L), P), a.multiply(L, P)), -1)
    assert marker == a.diag((F(1, 16), 0))
    commutator = a.add(a.multiply(P, L), a.multiply(L, P), -1)
    assert a.multiply(a.transpose(commutator), commutator) == a.scale(a.identity(2), F(1, 16))
    # Integrate the two exponential modes algebraically. The kernel's D-rate cancels.
    kernel_coefficients = {F(1, 2): B*B/(2*(D-F(1, 2))), F(1): B*B/(2*(D-1))}
    hidden_rate_coefficient = -sum(kernel_coefficients.values())
    assert hidden_rate_coefficient == 0
    assert kernel_coefficients == {F(1, 2): F(1, 8), F(1): F(-1, 8)}
    for rate, memory_coefficient in kernel_coefficients.items():
        assert -rate/2 == -A/2+memory_coefficient
    resolvents = []
    for z in (F(1, 10), F(1), F(3)):
        full = a.inverse(a.add(a.scale(a.identity(2), z), L))
        sigma = B*B/(z+D)
        assert full[0][0] == 1/(z+A-sigma)
        assert full[0][0] == (1/(z+F(1, 2))+1/(z+1))/2
        assert sigma > 0
        resolvents.append({'z': str(z), 'memory_self_energy': str(sigma)})
    rows = []
    for t in (F(1, 8), F(1, 2), F(1), F(2)):
        ea, eb = refine.ex(-t/2), refine.ex(-t)
        defect = (ea-eb)*(ea-eb)/4
        assert defect.lo > 0
        # The concrete tail is (1/12) exp(-3R/4) for a unit visible bound.
        bound = memory_rhs_bound(F(1, 4), F(3, 4), t, 1)
        expected = refine.ex(-3*t/4)/12
        assert bound.lo == expected.lo and bound.hi == expected.hi
        rows.append({'tau': str(t), 'closure_defect_lower': str(defect.lo),
                     'unit_source_memory_tail_upper': str(bound.hi)})
    return {'closure_marker': '1/16', 'initial_visible_rate': '3/4', 'retained_asymptotic_rate': '1/2',
            'kernel_coefficient': '1/16', 'kernel_decay': '3/4', 'resolvent_checks': resolvents,
            'strict_defect_and_tail_enclosures': rows}


def refusals():
    cases = {'float_turn': lambda: rational_step(0.5),
             'wrong_dimension': lambda: word_product((a.identity(3),)),
             'missing_history': lambda: eliminated_step((a.identity(2),), 0, [], 0),
             'boolean_step': lambda: eliminated_step((a.identity(2),), True, [1, 1], 0),
             'zero_hidden_floor': lambda: memory_rhs_bound(1, 0, 1, 1),
             'negative_age': lambda: memory_rhs_bound(1, 1, -1, 1)}
    out = {}
    for key, call in cases.items():
        try:
            call()
        except ValueError:
            out[key] = 'rejected'
        else:
            raise AssertionError('invalid input accepted: '+key)
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
    return {'certificate_type': 'GE3_NATIVE_RETURNING_MEMORY', 'verdict': 'PASS',
            'input_sha256': source_checks(), 'native_sector': native_sector_controls(),
            'clock_free_words': word_controls(), 'retained_and_reset_limits': refinement_controls(),
            'continuum_memory': continuum_controls(), 'refusals': refusals(),
            'scope': {'written_results': ['GE3-T1', 'GE3-T2', 'GE3-T3', 'GE3-T4'],
                'carrier': 'fixed finite native coefficient core; explicit centered quadratic sector',
                'certificate': 'written scoped proofs and exact finite controls, not formal or expert certification',
                'new': 'intrinsic unequal-step retained/reset limits and actual native polynomial witness',
                'inherited': 'block elimination, native heat, closure criterion and Schur calculus',
                'open': ['actual interacting row cut', 'general unbounded or moving cuts',
                         'physical protocol and seconds', '4D Yang-Mills', 'Clay']}}


def check(cert, result=RESULT, expected=EXPECTED):
    sha = y.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or expected.read_text().strip() != sha:
        raise ValueError('GE3 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run(); sha = y.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        EXPECTED.write_text(sha+'\n')
    else:
        check(cert)
    print('GE3 PASS', sha)
    print('4 written results; actual native quadratic memory; interacting row and Clay OPEN')


if __name__ == '__main__':
    main()
