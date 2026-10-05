"""YM83: exact controls for actual conditional gaps and global response.

Python 3.12 only. This checks finite arithmetic and source bytes, not the
universal analytic proofs. Actual Yang--Mills correlation tails remain open.
Default/--check are read-only; --write updates this stage's result and pin.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / 'YM83_RESULT.json'
PIN = HERE / 'EXPECTED_YM83.sha256'
SOURCES = HERE / 'YM83_SOURCE_MANIFEST.json'
BASE_COMMIT = 'a26fc6019102643f522d421d8a0ebce7b1714ec1'
BASE_TREE = '637f1596f03f8d4b9c827fac63bf63ab8ce34920'
PARENT_HASH = '0d5afd0ea881d9ceae39d7076ec3b82b4c73910ed55aa45f2b8f5f9421711304'
NSB = ROOT / 'papers/nonlinear-spectral-blindness/certificates'
MODEL = dict(
    actual_ground_state_measure_used=True,
    actual_conditional_gap_at_every_finite_amplitude=True,
    conditional_floor_uniform_in_exterior_and_total_volume=True,
    conditional_floor_uniform_in_amplitude=False,
    global_gap_requires_conditional_correlation_row_margin=True,
    finite_response_compression_alone_bounds_correlation_from_above=False,
    omitted_tail_and_exterior_uniformity_are_required=True,
    comparison_binary_models_are_Yang_Mills_vacua=False,
    quantum_Hamiltonian_is_conditional_resampling_generator=False,
    original_homogeneous_weak_coupling_correlation_bound_proved=False,
    continuum_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def exact(x):
    if isinstance(x, bool) or not isinstance(x, (int, str, F)):
        raise TypeError('Use an integer, rational string or Fraction')
    return F(x)


def nonnegative(x):
    x = exact(x)
    if x < 0:
        raise ValueError('Expected a nonnegative rational')
    return x


def positive(x):
    x = nonnegative(x)
    if not x:
        raise ValueError('Expected a positive rational')
    return x


def size(n, minimum=1):
    if type(n) is not int or n < minimum:
        raise ValueError('Invalid integer size')
    return n


def sqrt_interval(x, digits=12):
    x = nonnegative(x)
    size(digits)
    scale = 10 ** digits
    k = isqrt(x.numerator * scale * scale // x.denominator)
    lower = F(k, scale)
    return dict(lower=lower, upper=lower if lower * lower == x else F(k + 1, scale))


def exp_negative_interval(x, terms=32):
    """For moderate x, invert a positive Taylor sum with a geometric tail."""
    x = nonnegative(x)
    size(terms)
    if x >= terms + 2:
        raise ValueError('Exponential tail ratio is not below one')
    term = total = F(1)
    for k in range(1, terms + 1):
        term *= x / k
        total += term
    tail = term * x / (terms + 1) / (1 - x / (terms + 2))
    return dict(lower=1 / (total + tail), upper=1 / total)


def symbolic_exp_floor(prefactor, exponent_upper):
    """Strictly positive rational floor, encoded without enormous integers.

    e < 3 implies prefactor * exp(-x) >= prefactor * 3**(-ceil(x)).
    The expression is exact; it is not a numerical underflow or log proxy.
    """
    prefactor, exponent_upper = positive(prefactor), nonnegative(exponent_upper)
    return dict(prefactor=prefactor, denominator_base=3,
                denominator_power=exponent_upper.__ceil__())


def conditional_budget(links, omega):
    size(links)
    omega = nonnegative(omega)
    if omega == 0:
        return dict(links=links, omega=omega, exact_free_gap=F(3, 4))
    root = sqrt_interval(3 * links * omega)
    exponent_upper = 4 * F(22, 7) * root['upper']  # pi < 22/7
    pref9 = F(3, 4) * F(19, 21) ** (2 * links)
    prefopt = F(3, 4 * 64 ** links)
    return dict(links=links, omega=omega,
                fixed_time=dict(time=9, prefactor=pref9, exponent=18 * omega,
                                rational_floor=symbolic_exp_floor(pref9, 18 * omega)),
                optimized=dict(prefactor=prefopt, sqrt_argument=3 * links * omega,
                               sqrt_bracket=root, exponent_upper=exponent_upper,
                               rational_floor=symbolic_exp_floor(prefopt, exponent_upper)),
                actual_vacuum=True, exterior_uniform=True, total_volume_uniform=True,
                amplitude_uniform=False)


def original_model_budget(links, incident_faces, theta):
    size(incident_faces)
    theta = nonnegative(theta)
    return dict(theta=theta, incident_faces=incident_faces,
                budget=conditional_budget(links, 2 * incident_faces * theta))


def identity(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    bt = transpose(b)
    return [[sum((x * y for x, y in zip(row, col)), F(0)) for col in bt] for row in a]


def add(a, b, factor=1):
    return [[x + factor * y for x, y in zip(ar, br)] for ar, br in zip(a, b)]


def psd(a):
    """Exact symmetric Schur elimination, including zero pivots."""
    if not a or any(len(row) != len(a) for row in a):
        raise ValueError('Expected a nonempty square matrix')
    a = [[exact(x) for x in row] for row in a]
    if a != transpose(a):
        raise ValueError('Expected a symmetric matrix')
    while a:
        pivot = a[0][0]
        if pivot < 0:
            return False
        if not pivot:
            if any(a[0][1:]):
                return False
            a = [row[1:] for row in a[1:]]
        else:
            a = [[a[i][j] - a[i][0] * a[0][j] / pivot
                  for j in range(1, len(a))] for i in range(1, len(a))]
    return True


def response_budget(matrix, op_upper, chi_square_upper, measurement_error=0):
    """Exact marked matrix plus a total HS squared upper bound.

    measurement_error is an additional operator uncertainty. With measured
    coefficients, users must supply a separate omitted-tail bound instead
    of subtracting the measured Frobenius norm from chi_square_upper.
    """
    if not matrix or not matrix[0] or any(len(row) != len(matrix[0]) for row in matrix):
        raise ValueError('Expected a nonempty rectangular matrix')
    matrix = [[exact(x) for x in row] for row in matrix]
    op_upper = nonnegative(op_upper)
    error = nonnegative(measurement_error)
    total = nonnegative(chi_square_upper)
    gram = multiply(transpose(matrix), matrix)
    if not psd(add([[op_upper ** 2 * x for x in row] for row in identity(len(gram))], gram, -1)):
        raise ValueError('Proposed finite operator-norm upper bound is invalid')
    captured = sum((x * x for row in matrix for x in row), F(0))
    if captured > total:
        raise ValueError('Captured exact response exceeds total squared kernel bound')
    tail = sqrt_interval(total - captured)
    return dict(matrix=matrix, finite_operator_norm_upper=op_upper,
                captured_squared_norm=captured, total_squared_norm_upper=total,
                omitted_squared_norm_upper=total - captured,
                omitted_norm_bracket=tail, additional_operator_error=error,
                correlation_upper=min(F(1), op_upper + error + tail['upper']),
                actual_vacuum_estimate_supplied=False)


def global_budget(local_floors, correlations, kappa=1):
    floors = [positive(x) for x in local_floors]
    if not floors or len(correlations) != len(floors):
        raise ValueError('Local floors and correlation matrix disagree')
    c = [[nonnegative(x) for x in row] for row in correlations]
    if any(len(row) != len(floors) for row in c) or c != transpose(c):
        raise ValueError('Correlation matrix must be square and symmetric')
    if any(c[i][i] for i in range(len(c))) or any(x > 1 for row in c for x in row):
        raise ValueError('Zero diagonal and correlations in [0,1] are required')
    kappa = positive(kappa)
    r = max(sum(row) for row in c)
    accepted = r < 1
    bare = min(floors) * (1 - r) if accepted else None
    return dict(row_sum=r, row_margin=1 - r, accepted=accepted,
                conditional_bare_gap_lower=bare,
                conditional_physical_gap_lower=kappa * bare if accepted else None,
                actual_vacuum_correlation_hypothesis_discharged=False,
                failed_test_proves_gaplessness=False)


def probability_input(states, weights):
    if not states or len(states) != len(weights) or len(set(states)) != len(states):
        raise ValueError('Invalid finite states')
    dim = len(states[0])
    if not dim or any(len(s) != dim for s in states):
        raise ValueError('Invalid state dimension')
    axes = [sorted(set(s[i] for s in states)) for i in range(dim)]
    if set(product(*axes)) != set(states):
        raise ValueError('Strictly positive full Cartesian support is required')
    weights = [positive(w) for w in weights]
    if sum(weights) != 1:
        raise ValueError('Weights must sum exactly to one')
    return weights


def conditional_projection(states, weights, block):
    weights = probability_input(states, weights)
    dim = len(states[0])
    if not block or len(set(block)) != len(block) or any(type(i) is not int or i < 0 or i >= dim for i in block):
        raise ValueError('Invalid nonempty block')
    outside = [i for i in range(dim) if i not in block]
    keys = [tuple(s[i] for i in outside) for s in states]
    totals = {key: sum(w for w, k in zip(weights, keys) if k == key) for key in set(keys)}
    return [[weights[j] / totals[keys[i]] if keys[i] == keys[j] else F(0)
             for j in range(len(states))] for i in range(len(states))]


def resampling_generator(states, weights, blocks):
    flat = [i for block in blocks for i in block]
    if sorted(flat) != list(range(len(states[0]))):
        raise ValueError('Blocks must partition all coordinates')
    out = [[F(0) for _ in states] for _ in states]
    for block in blocks:
        out = add(out, add(identity(len(states)), conditional_projection(states, weights, block), -1))
    return out


def parity_law(n, t):
    size(n, 3)
    t = exact(t)
    if abs(t) >= 1:
        raise ValueError('Strict positivity requires |t|<1')
    states = list(product((-1, 1), repeat=n))
    weights = [(1 + t * prod(state)) / 2 ** n for state in states]
    return states, weights


def prod(values):
    out = 1
    for value in values:
        out *= value
    return out


def parity_control(n, t):
    states, weights = parity_law(n, t)
    k = resampling_generator(states, weights, [[i] for i in range(n)])
    for mask in range(1, 2 ** n - 1):
        subset = [i for i in range(n) if mask & (1 << i)]
        other = [i for i in range(n) if i not in subset]
        v = [prod(s[i] for i in subset) for s in states]
        vc = [prod(s[i] for i in other) for s in states]
        assert [sum(x * y for x, y in zip(row, v)) for row in k] == [len(subset) * (x - t * y) for x, y in zip(v, vc)]
        assert sum(w * value for w, value in zip(weights, v)) == 0
    eta = 1 - (n - 1) * abs(t)
    if eta > 0:
        defect = add(multiply(k, k), k, -eta)
        assert psd([[weights[i] * x for x in row] for i, row in enumerate(defect)])
    root = sqrt_interval((n - 2) ** 2 + 4 * t * t * (n - 1))
    return dict(kind='finite_probability_control_not_Yang_Mills', n=n, t=t,
                every_proper_marginal_uniform=True, conditional_pair_correlation=abs(t),
                exact_gap_interval=dict(lower=(n - root['upper']) / 2, upper=(n - root['lower']) / 2),
                row_test=global_budget([1] * n, [[0 if i == j else abs(t) for j in range(n)] for i in range(n)]))


def marker_controls(t=F(3, 4)):
    states, weights = parity_law(4, t)
    matrix = [[sum(w * s[i] * s[j] for s, w in zip(states, weights)) for j in (2, 3)] for i in (0, 1)]
    marked = [[sum(w * prod(s) for s, w in zip(states, weights))]]
    xi = sum((16 * w - 1) ** 2 for w in weights) / 16
    assert xi == t * t and all(not x for row in matrix for x in row)
    return dict(kind='finite_probability_control_not_Yang_Mills',
                invisible=response_budget(matrix, 0, xi),
                repaired=response_budget(marked, abs(t), xi))


def native_incidence(side):
    size(side, 4)
    def shift(v, axis):
        return tuple((x + (axis == i)) % side for i, x in enumerate(v))
    faces = []
    for v in product(range(side), repeat=3):
        for i, j in combinations(range(3), 2):
            faces.append(frozenset(((v, i), (shift(v, i), j), (shift(v, j), i), (v, j))))
    link = ((0, 0, 0), 0)
    face = faces[0]
    return dict(side=side, total_faces=len(faces),
                link_incident_faces=sum(link in p for p in faces),
                plaquette_block_links=len(face),
                plaquette_block_incident_faces=sum(bool(face & p) for p in faces))


def controls():
    heat_first = exp_negative_interval(F(27, 4))
    heat_ratio = F(9, 4) * exp_negative_interval(F(45, 4))['upper']
    assert heat_ratio < F(1, 2) and 8 * heat_first['upper'] < F(1, 20)
    incidence = [native_incidence(n) for n in (4, 6, 8)]
    assert all(row['link_incident_faces'] == 4 and row['plaquette_block_incident_faces'] == 13 for row in incidence)
    theta = F(1, 16000)
    small_floor = F(3, 4) * F(19, 21) ** 2 * (1 - 144 * theta)
    assert small_floor > F(3, 5)
    return dict(
        heat_time_nine=dict(deviation_upper=8 * heat_first['upper'], tail_ratio_upper=heat_ratio),
        actual_model_local_bounds=[original_model_budget(n, faces, theta) for n, faces in ((1, 4), (4, 13)) for theta in (F(0), F(1, 16000), F(1), F(25))],
        actual_small_coupling_conditional_link_floor=small_floor,
        native_lattice_incidence=incidence,
        conditional_global_example=global_budget([small_floor] * 3, [[0, F(1, 8), F(1, 8)], [F(1, 8), 0, F(1, 8)], [F(1, 8), F(1, 8), 0]], kappa=2),
        parity=[parity_control(n, t) for n in (3, 4) for t in (F(0), F(1, 4), F(-1, 4), F(3, 4), F(99, 100))],
        nonlinear_marker_repair=marker_controls())


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {str(k): encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def canonical_sha(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'YM83' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM83 source base mismatch')
    if len(manifest['upstream_sha256']) != 354 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM83 source manifest shape mismatch')
    if (NSB / 'EXPECTED_NSB1.sha256').read_text().strip() != PARENT_HASH or canonical_sha(json.loads((NSB / 'NSB1_RESULT.json').read_text())) != PARENT_HASH:
        raise ValueError('NSB1 result/pin mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = sha256((ROOT / path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed upstream input: ' + path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = sha256((ROOT / path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return encode(dict(certificate_type='YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE', verdict='PASS',
                       runtime_policy='Python 3.12 only', model=MODEL,
                       input_sha256=source_checks(), controls=controls(),
                       claim_status='ACTUAL_LOCAL_CONDITIONAL_FLOOR__GLOBAL_RESPONSE_CONDITIONAL__CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM83_ACTUAL_VACUUM_CONDITIONAL_RESPONSE':
        raise ValueError('YM83 scope mismatch')
    digest = canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('YM83 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--write', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2) + '\n')
        PIN.write_text(digest + '\n')
    else:
        check(cert)
    print('YM83 PASS', digest)
    print('Actual conditional floors: proved; global weak-coupling response and continuum gap: open')


if __name__ == '__main__':
    main()
