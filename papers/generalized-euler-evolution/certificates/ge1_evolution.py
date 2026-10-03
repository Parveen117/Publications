"""GE1 exact finite controls for the written domain/evolution theorems.

Python 3.12 only. Uses the unchanged native YM coefficient, rational matrix
and interval engines. Finite controls do not prove infinite-domain claims.
Default/--check is read-only; only --write updates this packet's evidence.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import factorial, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
ROOT = PACKET.parents[1]
YM = ROOT / 'papers/yang-mills-certified-benchmark/certificates'
sys.path.insert(0, str(YM))
import ym52_energy_observer as native

q, a, y, p = native.q, native.a, native.y, native.p
Iv = y.Iv
refine = y.y44
RESULT = HERE / 'GE1_RESULT.json'
EXPECTED = HERE / 'EXPECTED_GE1.sha256'
PINS = HERE / 'GE1_SOURCE_PINS.json'


def integer(n, lower=0):
    if isinstance(n, bool) or not isinstance(n, int) or n < lower:
        raise ValueError('integer outside the admitted range')
    return n


def square_matrix(A):
    if not A or any(len(row) != len(A) for row in A):
        raise ValueError('nonempty square matrix required')
    return [[q.exact(x) for x in row] for row in A]


def pairing(A, G):
    A, G = square_matrix(A), square_matrix(G)
    if len(A) != len(G) or G != a.transpose(G):
        raise ValueError('compatible symmetric pairing required')
    if a.y37.inertia(G) != (len(G), 0, 0):
        raise ValueError('strictly positive finite pairing required')
    return A, G


def zero(n):
    return a.scale(a.identity(n), 0)


def cayley(D, h, G):
    D, G = pairing(D, G)
    h = q.exact(h)
    GD = a.multiply(G, D)
    if a.add(a.transpose(GD), GD) != zero(len(D)):
        raise ValueError('skew generator required for a reversible phase step')
    I = a.identity(len(D))
    return a.multiply(a.inverse(a.add(I, a.scale(D, -h/2))),
                      a.add(I, a.scale(D, h/2)))


def resolvent(L, h, G):
    L, G = pairing(L, G)
    h = q.exact(h)
    if h < 0 or not a.psd(a.multiply(G, L)):
        raise ValueError('nonnegative heat step and positive symmetric generator required')
    return a.inverse(a.add(a.identity(len(L)), a.scale(L, h)))


def phase_budget(t, n, third_graph_norm):
    integer(n, 1)
    t, third_graph_norm = map(q.exact, (t, third_graph_norm))
    if third_graph_norm < 0:
        raise ValueError('nonnegative graph bound required')
    return abs(t)**3 * third_graph_norm / (12*n*n)


def heat_budget(t, n, second_graph_norm):
    integer(n, 1)
    t, second_graph_norm = map(q.exact, (t, second_graph_norm))
    if min(t, second_graph_norm) < 0:
        raise ValueError('nonnegative time and graph bound required')
    return 3*t*t*second_graph_norm / (2*n)


def change_residual(k, real, imag):
    if not k or len(k) != len(real) or len(k) != len(imag):
        raise ValueError('matched nonempty frequency and integer mode lists required')
    if any(isinstance(v, bool) or not isinstance(v, int) for v in k):
        raise ValueError('integer character labels required')
    real, imag = [q.exact(v) for v in real], [q.exact(v) for v in imag]
    return (1-sum((j*b for j, b in zip(k, imag)), F(0)),
            sum((j*r for j, r in zip(k, real)), F(0)))


def coefficient(k):
    integer(k, 1)
    r = isqrt(k)
    return F(1, 2**(r+(r*r != k)))


def forward_growth(m):
    integer(m, 1)
    n = m*m
    return (1+F(n*n, n*n))**n * coefficient(n)**2


def pell_pair(m):
    integer(m, 1)
    u, v = 1, 1
    for _ in range(1, m):
        u, v = u+2*v, u+v
    return u, v


def pell_gap_bracket(m):
    u, v = pell_pair(m)
    # 1 < sqrt(2) < 2 and |u^2-2v^2|=1, without subtractive cancellation.
    return Iv(F(1, (u+2*v)**2), F(1, (u+v)**2))


def finite_phase_reference(alpha, t, terms=40):
    alpha, t = map(q.exact, (alpha, t))
    integer(terms, 1)
    x = alpha*t
    c = sum(((-1)**j*x**(2*j)/factorial(2*j)
             for j in range(terms//2+1)), F(0))
    s = sum(((-1)**j*x**(2*j+1)/factorial(2*j+1)
             for j in range((terms-1)//2+1)), F(0))
    tail = refine.factorial_tail(abs(x), terms)
    return Iv(c-tail, c+tail), Iv(s-tail, s+tail)


def abs_ceiling(value):
    return max(abs(value.lo), abs(value.hi))


def phase_controls():
    exact_checks, rows = 0, []
    G = a.identity(2)
    for alpha in (F(-3), F(-1, 2), F(0), F(2)):
        D = [[F(0), -alpha], [alpha, F(0)]]
        for h in (F(-1), F(0), F(1, 4), F(3)):
            V = cayley(D, h, G)
            assert a.multiply(a.transpose(V), V) == G
            assert a.multiply(V, cayley(D, -h, G)) == G
            exact_checks += 2
        for t in (F(-1, 2), F(1, 2)):
            c, s = finite_phase_reference(alpha, t)
            for n in (2, 8, 32):
                W = a.power(cayley(D, t/n, G), n)
                err2 = abs_ceiling(Iv(W[0][0])-c)**2 + abs_ceiling(Iv(W[1][0])-s)**2
                bound = phase_budget(t, n, abs(alpha)**3)
                assert err2 <= bound*bound
                rows.append({'frequency': str(alpha), 'time': str(t), 'n': n,
                             'error_square_upper': str(err2), 'bound': str(bound)})
    return {'exact_isometry_and_reversal_checks': exact_checks,
            'independent_scalar_enclosures': rows}


def domain_controls():
    assert change_residual((1,), (0,), (1,)) == (0, 0)
    assert change_residual((-1,), (0,), (1,)) == (2, 0)
    real_checks = 0
    for k in range(-8, 9):
        r = change_residual((k,), (F(2, 3),), (0,))
        assert r[0] == 1
        real_checks += 1
    growth, taylor = [], []
    for m in (4, 8, 16, 32):
        g = forward_growth(m)
        assert g == F(2)**(m*m-2*m)
        growth.append({'m': m, 'cutoff_and_step_count': m*m,
                       'fixed_vector_coefficient_square_after_forward_steps': str(g)})
        lower = F(m, 2)**m
        exact_term = F(m*m)**m * coefficient(m*m)/factorial(m)
        assert exact_term >= lower
        taylor.append({'order': m, 'single_coefficient_term': str(exact_term),
                       'lower_bound': str(lower)})
    # beta_k=k for imaginary omega=iota: at t=1, k=-m grows as exp(m).
    complex_rows = []
    for m in (1, 4, 16):
        minus, plus = refine.ex(m), refine.ex(-m)
        assert minus.lo > 1 and 0 < plus.lo <= plus.hi < 1
        complex_rows.append({'negative_mode': -m, 'growth_lower': str(minus.lo),
                             'positive_mode': m, 'decay_upper': str(plus.hi)})
    return {'correct_complex_kernel_mode': 1, 'wrong_mode_residual': ['2', '0'],
            'real_kernel_exclusion_checks': real_checks,
            'fixed_smooth_vector_forward_failure': growth,
            'smooth_nonanalytic_taylor_failure': taylor,
            'bilateral_complex_domain_controls': complex_rows}


def native_controls():
    counts = {'native_noncommuting_brackets': 0, 'native_cayley_isometries': 0,
              'native_positive_resolvents': 0, 'inverse_equations': 0,
              'reference_pairing_contractions': 0}
    C0 = a.diag((0, 1, 1))
    # Rotate about the third axis so rank-two eigenspaces acquire off-diagonal entries.
    R = q.rotation((F(3, 5), F(0), F(0), F(4, 5)))
    C1 = q.rotate_tensor(C0, R)
    assert C1[0][1] != 0 and a.y37.inertia(C1) == (2, 0, 1)
    for d in (1, 2, 3):
        G = y.degree_data(d)['G']
        D = q.generators(d)
        assert a.add(a.multiply(D[0], D[1]), a.multiply(D[1], D[0]), -1) == D[2]
        assert D[2] != zero(len(G))
        counts['native_noncommuting_brackets'] += 1
        for h in (F(1, 8), F(2)):
            V = cayley(D[0], h, G)
            assert a.multiply(a.multiply(a.transpose(V), G), V) == G
            counts['native_cayley_isometries'] += 1
            for C in (C0, C1):
                L = q.operator(d, C)
                B = resolvent(L, h, G)
                assert a.multiply(a.add(a.identity(len(G)), a.scale(L, h)), B) == a.identity(len(G))
                assert a.multiply(B, a.add(a.identity(len(G)), a.scale(L, h))) == a.identity(len(G))
                counts['inverse_equations'] += 2
                assert a.psd(a.multiply(G, B))
                assert a.psd(a.add(G, a.multiply(a.multiply(a.transpose(B), G), B), -1))
                counts['native_positive_resolvents'] += 1
                H, _ = native.gram(d)
                assert a.psd(a.add(H, a.multiply(a.multiply(a.transpose(B), H), B), -1))
                counts['reference_pairing_contractions'] += 1
    return {'checks': counts, 'rotated_rank_two_tensor': [[str(v) for v in r] for r in C1],
            'degree_range': [1, 3], 'reused_compact_gap_for_diagonal_rank_two': '1/2'}


def heat_controls():
    rows, exact_modes = [], 0
    for C, rate1, rate2 in ((a.diag((0, 1, 1)), F(1, 2), F(1)),
                            (a.identity(3), F(3, 4), F(2))):
        for d, poly, rate in ((1, p.COORD[0], rate1), (2, native.EVEN, rate2)):
            v = y.vector(poly, d)
            L, G = q.operator(d, C), y.degree_data(d)['G']
            assert a.apply(L, v) == [rate*x for x in v]
            for t in (F(1, 4), F(1)):
                target = refine.ex(-t*rate)
                for n in (1, 4, 16, 64):
                    factor = (1+t*rate/n)**(-n)
                    B = resolvent(L, t/n, G)
                    assert a.apply(B, v) == [x/(1+t*rate/n) for x in v]
                    exact_modes += 1
                    err = Iv(factor)-target
                    bound = heat_budget(t, n, rate*rate)
                    assert 0 <= err.lo <= err.hi <= bound
                    rows.append({'degree': d, 'rate': str(rate), 'time': str(t), 'n': n,
                                 'error_upper': str(err.hi), 'bound_per_unit_norm': str(bound)})
    # Norm-stable midpoint heat steps need not preserve pointwise positivity.
    L3 = [[F(2) if i == j else F(-1) for j in range(3)] for i in range(3)]
    G3 = a.identity(3)
    mid = a.multiply(a.inverse(a.add(G3, a.scale(L3, 2))), a.add(G3, a.scale(L3, -2)))
    backward = resolvent(L3, F(4), G3)
    assert mid[0][0] == F(-1, 7) and all(v >= 0 for r in backward for v in r)
    assert a.psd(a.add(G3, a.multiply(a.transpose(mid), mid), -1))
    assert all(sum(r) == 1 for r in backward)
    return {'native_eigenmode_inverse_checks': exact_modes,
            'independent_heat_enclosures': rows,
            'midpoint_heat_negative_diagonal': '-1/7',
            'positive_resolvent_row_sums': ['1', '1', '1']}


def resonance_controls():
    rows, previous = [], None
    for m in range(1, 25):
        u, v = pell_pair(m)
        assert u*u-2*v*v == (-1)**m
        bound = pell_gap_bracket(m)
        assert 0 < bound.lo < bound.hi
        if previous is not None:
            assert bound.hi < previous.lo
        previous = bound
        rows.append({'m': m, 'p': u, 'q': v, 'pell_residue': (-1)**m,
                     'gap_lower': str(bound.lo), 'gap_upper': str(bound.hi)})
    assert rows[-1]['gap_upper'] != '0'
    # An exact rational-frequency average witness: alpha=0 retains a nonconstant mode.
    rational_stationary_mode = (2, -1)
    assert sum(k*w for k, w in zip(rational_stationary_mode, (1, 2))) == 0
    return {'near_resonant_modes': rows, 'nonconstant_rational_stationary_mode': [2, -1],
            'infinite_gap_statement': 'zero by the written Pell sequence proof; not finite sampling'}


def refusal_controls():
    G = a.identity(2)
    D = [[F(0), F(-1)], [F(1), F(0)]]
    cases = {
        'float_frequency': lambda: change_residual((1,), (0.5,), (0,)),
        'noninteger_mode': lambda: change_residual((F(1, 2),), (1,), (0,)),
        'boolean_count': lambda: phase_budget(1, True, 1),
        'zero_count': lambda: heat_budget(1, 0, 1),
        'negative_heat_time': lambda: heat_budget(-1, 2, 1),
        'negative_graph_bound': lambda: phase_budget(1, 2, -1),
        'non_square_matrix': lambda: cayley([[1, 0]], 1, G),
        'indefinite_pairing': lambda: cayley(D, 1, a.diag((1, -1))),
        'null_pairing': lambda: cayley(D, 1, a.diag((1, 0))),
        'heat_as_phase': lambda: cayley(G, 1, G),
        'phase_as_heat': lambda: resolvent(D, 1, G),
        'negative_generator': lambda: resolvent(a.scale(G, -1), 1, G),
        'negative_resolvent_step': lambda: resolvent(G, -1, G),
        'zero_pell_index': lambda: pell_pair(0),
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


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_checks(root=ROOT, pins=None):
    pins = json.loads(PINS.read_text()) if pins is None else pins
    checked = {}
    for path, expected in pins['upstream_sha256'].items():
        actual = digest(root/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in pins['local_inputs']:
        checked[path] = digest(root/path)
    if root == ROOT:
        checked[str(PINS.relative_to(ROOT))] = digest(PINS)
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or not __debug__:
        raise RuntimeError('Python 3.12 without optimization is required')
    return {'certificate_type': 'GE1_GENERALIZED_EULER_DOMAINS_AND_NATIVE_HEAT_DOCK',
            'verdict': 'PASS', 'runtime_policy': 'Python 3.12 only',
            'input_sha256': source_checks(),
            'phase_stability': phase_controls(), 'domain_counterexamples': domain_controls(),
            'actual_native_YM_actions': native_controls(), 'heat_approximation': heat_controls(),
            'near_resonance': resonance_controls(), 'refusals': refusal_controls(),
            'evidence_scope': {'written_theorems': ['GE1-T1', 'GE1-T2', 'GE1-T3', 'GE1-T4', 'GE1-T5', 'GE1-T6'],
                'finite_controls': 'exact rational algebra and outward enclosures, not an infinite proof assistant',
                'drafts': 'partially corrected and scoped; not certified wholesale',
                'YM_application': 'existing free compact heat generator, graph domain and resolvent approximation',
                'open': ['general phase-ratio/seam dynamics', 'physical state and clock selection',
                         'Born law', 'RH', 'interacting row closure', '4D Yang-Mills', 'Clay', 'quantum gravity']}}


def check(cert, result=RESULT, expected=EXPECTED):
    sha = y.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or expected.read_text().strip() != sha:
        raise ValueError('GE1 fresh certificate/pin mismatch')
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
    print('GE1 PASS', sha)
    print(json.dumps({'phase_enclosures': len(cert['phase_stability']['independent_scalar_enclosures']),
                      'heat_enclosures': len(cert['heat_approximation']['independent_heat_enclosures']),
                      'near_resonant_modes': len(cert['near_resonance']['near_resonant_modes']),
                      'refusals': len(cert['refusals']), '4D_and_Clay': 'OPEN'}, sort_keys=True))


if __name__ == '__main__':
    main()
