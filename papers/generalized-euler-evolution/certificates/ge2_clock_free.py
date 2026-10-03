"""GE2: exact controls for clock-free recognition and intrinsic refinement.

Python 3.12 only. Reuses the pinned native quaternion, polynomial, pairing,
matrix and interval engines. Written general proofs are separate evidence.
Default/--check is read-only; --write changes only GE2 evidence.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys

import ge1_evolution as old

q, a, y, p = old.q, old.a, old.y, old.p
native, Iv, refine = old.native, old.Iv, old.refine
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PINS = HERE/'GE2_SOURCE_PINS.json'
RESULT = HERE/'GE2_RESULT.json'
EXPECTED = HERE/'EXPECTED_GE2.sha256'


def turn(v):
    if len(v) != 3:
        raise ValueError('three native imaginary coefficients required')
    v = tuple(q.exact(x) for x in v)
    r2 = a.dot(v, v)
    denominator = 1+r2/16
    return ((1-r2/16)/denominator,)+tuple(-x/(2*denominator) for x in v)


def action(d, v):
    old.integer(d)
    return y.action_matrix(d, y.qmatrix(turn(v)))


def dagger(A, source_G, target_G=None):
    target_G = source_G if target_G is None else target_G
    if len(A) != len(target_G) or any(len(row) != len(source_G) for row in A):
        raise ValueError('arrow and pairing dimensions disagree')
    return a.multiply(a.multiply(a.inverse(source_G), a.transpose(A)), target_G)


def protocol(d, records):
    """Each (w,v) labels the symmetric pair +v,-v of weights w/2 each."""
    old.integer(d)
    records = tuple((w, v) for w, v in q.protocol(records) if w)
    Q, m2, m4 = q.moments(records)
    h = m2/2
    C = a.scale(Q, 1/m2) if m2 else None
    labels = tuple((w/2, tuple(sign*x for x in v))
                   for w, v in records for sign in (1, -1))
    actions = tuple((w, action(d, v)) for w, v in labels)
    S = old.zero(len(y.basis(d)))
    for w, U in actions:
        S = a.add(S, a.scale(U, w))
    eps2 = max(a.dot(v, v) for _, v in records)
    return {'records': records, 'labels': labels, 'actions': actions,
            'Q': Q, 'h': h, 'm4': m4, 'C': C, 'eps2': eps2, 'S': S}


def normalized_shape(item):
    if not item['h']:
        raise ValueError('idle arrow has no normalized shape')
    return item['C']


def constants(d):
    old.integer(d)
    return F(d*d*(d*d+2), 384), F(d*d*(2*d*d+1), 96)


def memory(S, G):
    return a.add(a.identity(len(S)), a.multiply(dagger(S, G), S), -1)


def norm_bounded(R, G, bound):
    bound = q.exact(bound)
    if bound < 0:
        raise ValueError('nonnegative norm bound required')
    return a.psd(a.add(a.scale(G, bound*bound),
                      a.multiply(a.multiply(a.transpose(R), G), R), -1))


def lift(weighted_actions, G):
    """Finite adapter only: constructs T24's typed source/target and cut."""
    weighted_actions = tuple(weighted_actions)
    if not weighted_actions:
        raise ValueError('nonempty record list required')
    actions = []
    for w, U in weighted_actions:
        w = q.exact(w)
        U, checked_G = old.pairing(U, G)
        if w <= 0 or a.multiply(a.multiply(a.transpose(U), checked_G), U) != checked_G:
            raise ValueError('strict positive record weights and isometries required')
        actions.append((w, U))
    if sum(w for w, _ in actions) != 1:
        raise ValueError('record weights must sum to one')
    n, count = len(G), len(actions)
    I = a.identity(n)
    A = [row[:] for _, U in actions for row in U]
    J = [row[:] for _ in actions for row in I]
    Jdag = [[actions[s][0]*I[i][j] for s in range(count) for j in range(n)]
            for i in range(n)]
    W = [[actions[s][0]*G[i][j] if s == t else F(0)
          for t in range(count) for j in range(n)]
         for s in range(count) for i in range(n)]
    P = a.multiply(J, Jdag)
    S = a.multiply(Jdag, A)
    return {'A': A, 'J': J, 'Jdag': Jdag, 'W': W, 'P': P, 'S': S}


def step_controls():
    rows, jets = [], []
    fixtures = (
        ((F(1), (F(1, 4), F(0), F(0))),),
        ((F(1, 2), (0, F(1, 3), 0)), (F(1, 2), (0, 0, F(1, 3)))),
        ((F(1, 3), (F(1, 2), F(1, 3), 0)),
         (F(2, 3), (0, F(1, 4), F(1, 2)))),
    )
    for d in (1, 2, 3):
        G = y.degree_data(d)['G']
        I = a.identity(len(G))
        B, K = constants(d)
        for records in fixtures:
            z = protocol(d, records)
            S, h = z['S'], z['h']
            L = q.operator(d, normalized_shape(z))
            M = memory(S, G)
            for _, v in z['labels']:
                u = turn(v)
                assert y.qmul(u, y.qdagger(u)) == y.ONE
                assert turn(tuple(-x for x in v)) == y.qdagger(u)
            for _, U in z['actions']:
                assert a.multiply(a.multiply(a.transpose(U), G), U) == G
            assert dagger(S, G) == S and a.psd(a.multiply(G, M))
            assert sum(z['C'][i][i] for i in range(3)) == 1
            residual = a.add(a.add(S, I, -1), a.scale(L, h))
            assert norm_bounded(residual, G, B*z['m4'])
            density = a.add(a.scale(M, 1/(2*h)), L, -1)
            bound = B*z['m4']/h+F(d**4, 32)*h
            assert norm_bounded(density, G, bound)
            assert bound <= K*z['eps2']
            assert norm_bounded(a.add(I, S, -1), G, h*d*d/4)
            rows.append({'degree': d, 'h': str(h), 'm4': str(z['m4']),
                         'step_bound': str(B*z['m4']), 'density_bound': str(bound),
                         'uniform_density_bound': str(K*z['eps2'])})
        # Independent interpolation of the exact rational substitution numerator.
        # Its degree is 2d, so 2d+1 exact evaluations recover the first two jets.
        v = (F(2, 3), F(-1, 2), F(1, 5))
        r2 = a.dot(v, v)
        nodes = list(range(-d, d+1))
        inverse = a.inverse([[F(s)**k for k in range(2*d+1)] for s in nodes])
        values = [a.scale(action(d, tuple(s*x for x in v)), (1+s*s*r2/16)**d)
                  for s in nodes]
        coefficients = []
        for k in (1, 2):
            out = old.zero(len(G))
            for coefficient, value in zip(inverse[k], values):
                out = a.add(out, a.scale(value, coefficient))
            coefficients.append(out)
        D = old.zero(len(G))
        for x, Dj in zip(v, q.generators(d)):
            D = a.add(D, a.scale(Dj, x))
        assert coefficients[0] == D
        second = a.add(coefficients[1], a.scale(I, d*r2/16), -1)
        assert second == a.scale(a.multiply(D, D), F(1, 2))
        jets.append({'degree': d, 'exact_interpolation_nodes': len(nodes)})
    return {'exact_operator_bound_controls': rows, 'independent_phase_and_second_jets': jets}


def record_controls():
    d = 1
    G = y.degree_data(d)['G']
    n = len(G)
    z1 = protocol(d, ((F(1), (F(1, 2), 0, 0)),))
    z2 = protocol(d, ((F(1), (0, F(1, 3), 0)),))
    first = lift(z1['actions'], G)
    # Direct enumeration of all four independent words.
    words = tuple((w*v, a.multiply(V, U))
                  for w, U in z1['actions'] for v, V in z2['actions'])
    whole = lift(words, G)
    append = [[V[i][j] if s == t else F(0)
               for t in range(2) for j in range(n)]
              for s in range(2) for _, V in z2['actions'] for i in range(n)]
    A1, W1, P1 = first['A'], first['W'], first['P']
    A2, W2, P2 = whole['A'], whole['W'], whole['P']
    assert a.multiply(append, A1) == A2
    assert a.multiply(a.multiply(a.transpose(append), W2), append) == W1
    no_return = a.multiply(a.multiply(P2, append), a.add(a.identity(2*n), P1, -1))
    assert not any(x for row in no_return for x in row)
    for model in (first, whole):
        A, W, P = model['A'], model['W'], model['P']
        assert a.multiply(a.multiply(a.transpose(A), W), A) == G
        assert a.multiply(P, P) == P
        assert a.multiply(a.transpose(P), W) == a.multiply(W, P)
        delta = a.add(a.multiply(a.multiply(dagger(A, G, W), P), A), a.identity(n), -1)
        assert delta == a.scale(memory(model['S'], G), -1)
    S1, S2 = z1['S'], z2['S']
    assert whole['S'] == a.multiply(S2, S1)
    pulled_step = a.multiply(a.multiply(dagger(append, W1, W2), P2), append)
    delta_step = a.add(pulled_step, P1, -1)
    delta_first = a.scale(memory(S1, G), -1)
    delta_total = a.add(delta_first, a.multiply(a.multiply(dagger(A1, G, W1), delta_step), A1))
    assert delta_total == a.scale(memory(whole['S'], G), -1)
    telescope = a.add(memory(S1, G), a.multiply(a.multiply(dagger(S1, G), memory(S2, G)), S1))
    assert telescope == memory(whole['S'], G)
    return {'retained_words': len(words), 'typed_chain_law': True,
            'no_return_on_old_zero_mean_records': True, 'memory_telescope': True,
            'first_source_dimension': n, 'final_event_dimension': len(W2)}


def refinement_controls():
    rows = []
    # Unequal squared turn sizes; their sum is exactly fixed across refinements.
    for d in (1, 2):
        G = y.degree_data(d)['G']
        poly = p.COORD[0] if d == 1 else native.EVEN
        v = y.vector(poly, d)
        C = a.diag((0, F(1, 2), F(1, 2)))
        rate = F(1, 4) if d == 1 else F(1, 2)
        assert a.apply(q.operator(d, C), v) == [rate*x for x in v]
        for n in (1, 2, 4, 8):
            # n^2 copies of turns 1/(2n) and 1/(3n) => tau=13/72.
            factors, steps = [], []
            for denominator in (2*n, 3*n):
                r = F(1, denominator)
                z = protocol(d, ((F(1, 2), (0, r, 0)), (F(1, 2), (0, 0, r))))
                # Independent scalar substitution formulas, not matrix eigenvalues.
                c, s = turn((r, 0, 0))[:2]
                scalar = c if d == 1 else c*c
                assert c*c+s*s == 1
                assert a.apply(z['S'], v) == [scalar*x for x in v]
                factors.append(scalar**(n*n))
                steps.append(z)
            tau = n*n*sum(z['h'] for z in steps)
            assert tau == F(13, 72)
            eps2 = max(z['eps2'] for z in steps)
            product = factors[0]*factors[1]
            error = old.abs_ceiling(Iv(product)-refine.ex(-tau*rate))
            bound = constants(d)[1]*eps2*tau
            assert error <= bound
            rows.append({'degree': d, 'n': n, 'steps': 2*n*n, 'tau': str(tau),
                         'epsilon_squared': str(eps2), 'error_upper': str(error),
                         'bound_per_unit_norm': str(bound)})
    phase_rows = []
    for n in (1, 2, 4, 8):
        U = action(1, (F(1, n), 0, 0))
        W = a.power(U, n)
        # Degree-one rotation entries: independently enclose cos(1/2), sin(1/2).
        c, s = old.finite_phase_reference(F(1, 2), 1)
        D = q.generators(1)[0]
        I = a.identity(4)
        error_square = F(0)
        for i in range(4):
            expected = c*I[i][0]+s*(2*D[i][0])
            error_square += old.abs_ceiling(Iv(W[i][0])-expected)**2
        bound = F(1, 96*n*n)
        assert error_square <= bound*bound
        phase_rows.append({'n': n, 'error_square_upper': str(error_square), 'bound': str(bound)})
    return {'irregular_heat_enclosures': rows, 'deterministic_phase_enclosures': phase_rows}


def distinction_controls():
    U = [[F(3, 5), F(-4, 5)], [F(4, 5), F(3, 5)]]
    P = a.diag((1, 0)); Q = a.diag((0, 1))
    full = a.multiply(a.multiply(a.multiply(P, a.transpose(U)), U), P)
    compressed = a.multiply(a.multiply(a.multiply(a.multiply(P, a.transpose(U)), P), U), P)
    returned = a.multiply(a.multiply(a.multiply(a.multiply(P, a.transpose(U)), Q), U), P)
    assert full == P and compressed == a.scale(P, F(9, 25))
    assert returned == a.scale(P, F(16, 25)) and a.add(compressed, returned) == full
    assert memory(U, a.identity(2)) == old.zero(2) and U != a.identity(2)
    assert memory(a.scale(U, -1), a.identity(2)) == memory(U, a.identity(2))
    v = (F(1, 2), 0, 0)
    assert y.qmul(turn(v), turn(tuple(-x for x in v))) == y.ONE
    loop_budget = a.dot(v, v)  # two successive deterministic record arrows
    assert loop_budget == F(1, 4) > 0
    C0 = a.diag((0, F(1, 2), F(1, 2)))
    rotation = q.rotation((F(3, 5), 0, 0, F(4, 5)))
    C1 = q.rotate_tensor(C0, rotation)
    L0, L1 = q.operator(2, C0), q.operator(2, C1)
    commutator = a.add(a.multiply(L1, L0), a.multiply(L0, L1), -1)
    witness = next((i, j, x) for i, row in enumerate(commutator) for j, x in enumerate(row) if x)
    records = ((F(1, 2), (0, F(1, 3), 0)), (F(1, 2), (0, 0, F(1, 3))))
    rotated = tuple((w, tuple(a.apply(rotation, v))) for w, v in records)
    left, right = protocol(2, records), protocol(2, rotated)
    assert left['h'] == right['h'] and left['C'] == C0 and right['C'] == C1
    assert a.multiply(left['S'], right['S']) != a.multiply(right['S'], left['S'])
    f = p.COORD[0]
    defect = p.add(p.scale(q.lap(p.mul(f, f), C0), -1), p.scale(p.mul(f, q.lap(f, C0)), 2))
    assert defect == p.scale(native.gamma(f, f, C0), 2) and defect
    expected_gamma = p.scale(p.add(p.mul(p.COORD[2], p.COORD[2]),
                                  p.mul(p.COORD[3], p.COORD[3])), F(1, 8))
    assert native.gamma(f, f, C0) == expected_gamma
    return {'returning_memory': '16/25', 'compressed_inverse_product': '9/25',
            'closed_endpoint_positive_ledger': str(loop_budget),
            'ordered_generator_commutator_entry': [witness[0], witness[1], str(witness[2])],
            'rational_word_order_detected': True, 'heat_is_not_a_derivation': True,
            'zero_memory_nontrivial_phase': True}


def curvature_controls():
    rows = []
    # Actual rational shape vectors: G_ij=2 dot(z_i,z_j), f=-cross(z_1,z_2)/2.
    for z1, z2 in (((F(1), 0), (0, F(1))), ((F(1), 0), (F(3, 5), F(4, 5))),
                   ((F(1), 0), (F(4, 5), F(3, 5))), ((F(1), 0), (F(1), 0))):
        G = [[2*a.dot(u, v) for v in (z1, z2)] for u in (z1, z2)]
        T = G[0][0]+G[1][1]
        f = -(z1[0]*z2[1]-z1[1]*z2[0])/2
        assert G[0][0]*G[1][1]-G[0][1]*G[1][0] == 16*f*f
        kappa = 64*f*f/(T*T)
        # Equal unit length fixture gives rational eigenvalues 1 +/- dot(z1,z2).
        b = (1-a.dot(z1, z2))/2
        c = 1-b
        assert 4*b*c == kappa and 0 <= kappa <= 1 and b >= kappa/4
        rate = native.rates((F(0), b, c))[0]
        assert rate == min(F(1, 4), b)
        assert native.rates((F(0), 2*b, 2*c))[0] == 2*rate
        rows.append({'T': str(T), 'f': str(f), 'kappa': str(kappa),
                     'normalized_eigenvalues': ['0', str(b), str(c)], 'record_clock_rate': str(rate)})
    return {'native_response_fixtures': rows, 'spectral_input': 'unchanged YM52/YM54'}


def refusal_controls():
    I = a.identity(2)
    cases = {
        'float_turn': lambda: turn((0.5, 0, 0)),
        'boolean_turn': lambda: turn((True, 0, 0)),
        'wrong_turn_dimension': lambda: turn((1, 0)),
        'boolean_degree': lambda: action(True, (1, 0, 0)),
        'negative_weight': lambda: protocol(1, ((-1, (1, 0, 0)), (2, (0, 1, 0)))),
        'unnormalized_weights': lambda: protocol(1, ((2, (1, 0, 0)),)),
        'empty_records': lambda: protocol(1, ()),
        'idle_shape_division': lambda: normalized_shape(protocol(1, ((1, (0, 0, 0)),))),
        'nonisometric_record': lambda: lift(((1, a.scale(I, 2)),), I),
        'zero_weight_target_pairing': lambda: lift(((0, I), (1, I)), I),
        'indefinite_source_pairing': lambda: lift(((1, I),), a.diag((1, -1))),
        'negative_norm_bound': lambda: norm_bounded(I, I, -1),
    }
    result = {}
    for label, call in cases.items():
        try:
            call()
        except ValueError:
            result[label] = 'rejected'
        else:
            raise AssertionError('invalid input accepted: '+label)
    return result


def source_checks(root=ROOT, pins=None):
    pins = json.loads(PINS.read_text()) if pins is None else pins
    checked = old.source_checks(root, pins)
    # old.source_checks also hashes GE1's pins at the real repository root.
    if root == ROOT:
        checked[str(PINS.relative_to(ROOT))] = old.digest(PINS)
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or not __debug__:
        raise RuntimeError('Python 3.12 without optimization is required')
    return {'certificate_type': 'GE2_CLOCK_FREE_RECOGNITION_GENERATOR', 'verdict': 'PASS',
            'runtime_policy': 'Python 3.12 only', 'input_sha256': source_checks(),
            'rational_steps': step_controls(), 'typed_records': record_controls(),
            'intrinsic_refinement': refinement_controls(), 'necessary_distinctions': distinction_controls(),
            'curvature_clock_corollary': curvature_controls(), 'refusals': refusal_controls(),
            'evidence_scope': {'written_theorems': ['GE2-T'+str(i) for i in range(1, 7)],
                'finite_controls': 'exact arithmetic and outward scalar enclosures; not infinite formal certification',
                'new_bridge': 'typed independent record arrows and intrinsic quadratic-ledger refinement',
                'inherited': 'native algebra, cut law, reference, variance/energy, heat and compact rates',
                'open': ['physical clock and protocol selection', 'general feedback and moving cuts',
                         'interacting row closure', '4D Yang-Mills', 'Clay', 'RH', 'quantum gravity']}}


def check(cert, result=RESULT, expected=EXPECTED):
    sha = y.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or expected.read_text().strip() != sha:
        raise ValueError('GE2 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = y.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        EXPECTED.write_text(sha+'\n')
    else:
        check(cert)
    print('GE2 PASS', sha)
    print(json.dumps({'written_results': 6, 'degree_bound_controls': len(cert['rational_steps']['exact_operator_bound_controls']),
                      'heat_enclosures': len(cert['intrinsic_refinement']['irregular_heat_enclosures']),
                      'phase_enclosures': len(cert['intrinsic_refinement']['deterministic_phase_enclosures']),
                      'refusals': len(cert['refusals']), '4D_and_Clay': 'OPEN'}, sort_keys=True))


if __name__ == '__main__':
    main()
