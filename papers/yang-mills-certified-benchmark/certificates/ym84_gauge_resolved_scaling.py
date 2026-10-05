"""YM84: exact controls for gauge sectors, response tails and physical scaling.

Python 3.12 only. Finite arithmetic controls are not a continuum proof.
The positive polynomial densities and Z2 cycle are comparison fixtures,
not numerical representations of the full SU(2) Yang--Mills vacuum.
Default/--check are read-only; --write updates only this stage's result/pin.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import comb
from pathlib import Path
import sys

import ym83_actual_vacuum_response as prev

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / 'YM84_RESULT.json'
PIN = HERE / 'EXPECTED_YM84.sha256'
SOURCES = HERE / 'YM84_SOURCE_MANIFEST.json'
BASE_COMMIT = '455a073936efa9b956b756e8feeccacadb3caa87'
BASE_TREE = '88babd8c6f523596f46903adc839829b1a5fb23c'
PARENT_HASH = '156fdb6d34878c00c6b5d3cb8dc64b78f9dd2eeb2ec08f5011840c56e50989ae'
MODEL = dict(
    actual_homogeneous_vacuum_open_path_response_identity=True,
    actual_open_path_rank_obstruction_uniform_in_volume=True,
    open_path_markers_annihilated_by_Gauss_projection=True,
    unconditioned_response_substituted_for_exterior_conditioned_response=False,
    one_Wilson_expectation_is_a_maximal_correlation_upper_bound=False,
    physical_sector_projection_comparison_proved=True,
    isolated_original_SU2_plaquette_physical_quantum_gap_bounds_proved=True,
    isolated_plaquette_vacuum_is_full_lattice_conditional_vacuum=False,
    polynomial_and_Z2_controls_are_actual_Yang_Mills_vacua=False,
    full_lattice_physical_local_sqrt_theta_floor_proved=False,
    full_lattice_physical_response_margin_of_order_a_proved=False,
    nonzero_continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def polynomial(values):
    p = [prev.exact(x) for x in values]
    if not p:
        raise ValueError('A polynomial needs at least one coefficient')
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def poly_add(a, b, factor=1):
    a, b, factor = polynomial(a), polynomial(b), prev.exact(factor)
    return polynomial([(a[i] if i < len(a) else 0)
                       + factor * (b[i] if i < len(b) else 0)
                       for i in range(max(len(a), len(b)))])


def poly_multiply(a, b):
    a, b = polynomial(a), polynomial(b)
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return polynomial(out)


def poly_value(p, x):
    x = prev.exact(x)
    result = F(0)
    for coefficient in reversed(polynomial(p)):
        result = result * x + coefficient
    return result


def haar_moment(degree):
    """W = Tr(U)/2 under normalized SU(2) Haar measure."""
    prev.size(degree, 0)
    if degree % 2:
        return F(0)
    k = degree // 2
    return F(comb(2 * k, k), (k + 1) * 4 ** k)


def haar_integral(p):
    return sum((x * haar_moment(i) for i, x in enumerate(polynomial(p))), F(0))


def character(n):
    """Spin n/2 character U_n(W), with Casimir n(n+2)/4."""
    prev.size(n, 0)
    before, current = [F(0)], [F(1)]
    for _ in range(n):
        before, current = current, poly_add(poly_multiply([0, 2], current), before, -1)
    return current


def central_response(unnormalized_density, cutoff=0):
    """Exact finite polynomial convolution; positivity is supplied separately.

    This routine does not certify an arbitrary polynomial as nonnegative.
    The calling fixtures provide explicit strictly positive constructions.
    """
    prev.size(cutoff, 0)
    p = polynomial(unnormalized_density)
    z = prev.positive(haar_integral(p))
    q = [x / z for x in p]
    coefficients = {n: haar_integral(poly_multiply(q, character(n))) / (n + 1)
                    for n in range(1, len(q))}
    captured = sum(((n + 1) ** 2 * value ** 2
                    for n, value in coefficients.items() if n <= cutoff), F(0))
    full = haar_integral(poly_multiply(q, q)) - 1
    spectral = sum(((n + 1) ** 2 * value ** 2 for n, value in coefficients.items()), F(0))
    assert full == spectral
    return dict(normalization=z, normalized_density=q,
                lambda_by_twice_spin=coefficients,
                correlation=max((abs(x) for x in coefficients.values()), default=F(0)),
                delta=1 - haar_integral(poly_multiply(q, [0, 1])),
                total_hs_squared=full, captured_hs_squared=captured,
                omitted_hs_squared=full - captured, cutoff=cutoff,
                positivity_certified_by_this_routine=False)


def concentrated_control(k, eta=F(3, 4), cutoff=2):
    prev.size(k)
    eta = prev.exact(eta)
    if not 0 < eta < 1:
        raise ValueError('Use 0 < eta < 1 for strict positivity')
    p = [F(1)]
    for _ in range(2 * k):
        p = poly_multiply(p, [1, eta])
    out = central_response(p, cutoff)
    out.update(kind='positive_SU2_polynomial_comparison_not_Yang_Mills',
               k=k, eta=eta, density_lower=(1 - eta) ** (2 * k) / out['normalization'])
    assert out['density_lower'] > 0
    for n, value in out['lambda_by_twice_spin'].items():
        assert 0 <= 1 - value <= F(n * (n + 2), 3) * out['delta']
    return out


def hidden_channel_control(n=6):
    prev.size(n, 2)
    amplitude = F(1, 2 * (n + 1))
    out = central_response(poly_add([1], character(n), amplitude), n - 1)
    # |chi_n| <= n+1 on SU(2), hence q >= 1/2 everywhere.
    assert out['normalization'] == 1 and out['lambda_by_twice_spin'][1] == 0
    assert out['correlation'] == amplitude / (n + 1)
    assert out['captured_hs_squared'] == 0 and out['omitted_hs_squared'] == amplitude ** 2
    out.update(kind='positive_SU2_polynomial_comparison_not_Yang_Mills',
               density_lower=F(1, 2), invisible_twice_spin=n)
    return out


def channel_dimension(n):
    prev.size(n, 0)
    return n * (2 * n * n + 9 * n + 13) // 6


def tail_budget(n, delta_upper, rank):
    """Universal consequence if the supplied delta is an actual upper bound."""
    prev.size(n)
    prev.size(rank, 0)
    delta = prev.nonnegative(delta_upper)
    if n * (n + 2) * delta > F(3, 2):
        raise ValueError('The retained channels are not certified above 1/2')
    dimension = channel_dimension(n)
    return dict(cutoff=n, delta_upper=delta, dimension=dimension, rank=rank,
                omitted_hs_squared_lower=F(max(0, dimension - rank), 4))


def actual_fourth_power_budget(s, rank=0):
    """Exact theta=s^4 sample of the written theorem for all theta>=256."""
    s = prev.exact(s)
    if s < 4:
        raise ValueError('Requires theta >= 256')
    n = (s / 2).__floor__()
    out = tail_budget(n, F(9, 4) / s ** 2, rank)
    assert n >= s / 4 and out['dimension'] >= s ** 3 / 192
    out.update(theta=s ** 4, xi_lower=s ** 3 / 768,
               all_theta_proof_is_analytic=True, actual_homogeneous_vacuum=True)
    return out


def oscillator_constants():
    """Rational comparisons used in the analytic SU(2) radial proof."""
    pi_lower, pi_upper = F(3), F(22, 7)
    root_lower, root_upper = F(140, 99), F(10, 7)
    c0_squared_upper = 2 * pi_upper ** 2 / 3 - 1
    assert root_lower ** 2 < 2 < root_upper ** 2
    assert c0_squared_upper < F(19, 8) ** 2
    ground_trial_admissibility = pi_lower ** 2 / 6 - F(1, 4) - 1
    two_mode_trial_admissibility = pi_lower ** 2 - 8
    lower = 7 * root_lower / pi_upper - F(19, 8)
    upper = 2 * root_upper * pi_upper
    assert ground_trial_admissibility > 0 and two_mode_trial_admissibility > 0
    assert lower > F(3, 4) and upper < 9
    return dict(pi_bracket=[pi_lower, pi_upper], sqrt_two_bracket=[root_lower, root_upper],
                ground_constant_squared_upper=c0_squared_upper,
                gap_coefficient_lower=lower, gap_lower_slack=lower - F(3, 4),
                gap_coefficient_upper=upper, gap_upper_slack=9 - upper,
                trial_admissibility_slacks=[ground_trial_admissibility, two_mode_trial_admissibility],
                radial_minmax_proof_is_analytic=True,
                model='isolated_original_four_link_SU2_cycle_not_full_lattice')


def apply(a, v):
    return [sum((x * y for x, y in zip(row, v)), F(0)) for row in a]


def cycle_matrices(t):
    """Exact Z2 four-cycle gauge action, separate from SU(2) claims."""
    t = prev.exact(t)
    states, weights = prev.parity_law(4, t)
    index = {s: i for i, s in enumerate(states)}
    gauge = [[F(0) for _ in states] for _ in states]
    for i, s in enumerate(states):
        for g in product((-1, 1), repeat=4):
            changed = tuple(g[k] * s[k] * g[(k + 1) % 4] for k in range(4))
            gauge[i][index[changed]] += F(1, 16)
    p1 = prev.conditional_projection(states, weights, [0, 1])
    p2 = prev.conditional_projection(states, weights, [2, 3])
    k = prev.add(prev.add([[2 * x for x in row] for row in prev.identity(16)], p1, -1), p2, -1)
    return dict(states=states, weights=weights, gauge=gauge, p1=p1, p2=p2,
                p0=[list(weights) for _ in states], k=k)


def cycle_control(t):
    t = prev.exact(t)
    data = cycle_matrices(t)
    p1, p2, p0, gauge, k = (data[x] for x in ('p1', 'p2', 'p0', 'gauge', 'k'))
    assert prev.multiply(gauge, gauge) == gauge
    assert prev.multiply(p1, gauge) == prev.multiply(gauge, p1) == p0
    assert prev.multiply(p2, gauge) == prev.multiply(gauge, p2) == p0
    assert prev.multiply(k, gauge) == [[2 * x for x in row] for row in prev.add(gauge, p0, -1)]
    states = data['states']
    left, right = ([s[0] * s[1] for s in states], [s[2] * s[3] for s in states])
    assert not any(apply(gauge, left)) and not any(apply(gauge, right))
    for sign in (-1, 1):
        v = [x + sign * y for x, y in zip(left, right)]
        assert apply(k, v) == [(1 - sign * t) * x for x in v]
    loop_centered = [prev.prod(s) - t for s in states]
    assert apply(k, loop_centered) == [2 * x for x in loop_centered]
    return dict(kind='finite_Z2_gauge_cycle_comparison_not_SU2_vacuum', t=t,
                full_resampling_gap=1 - abs(t), physical_resampling_gap=2,
                open_path_projection_zero=True, physical_pair_response_zero=True,
                gauge_projection_rank=sum(gauge[i][i] for i in range(16)))


def scale_budget(a, c, d, eta, energy_upper, variance_lower):
    """Conditional compatibility; this supplies no actual lattice inputs."""
    a, c, d, eta, energy, variance = [prev.positive(x) for x in
                                    (a, c, d, eta, energy_upper, variance_lower)]
    if eta > 1:
        raise ValueError('A positive row margin is at most one')
    lower, upper = c * d * eta / a, energy / variance
    return dict(a=a, c=c, d=d, eta=eta,
                conditional_gap_lower=lower, survival_gap_upper=upper,
                margin_upper_for_survival=energy * a / (c * d * variance),
                hypotheses_compatible=lower <= upper,
                actual_full_lattice_hypotheses_discharged=False)


def controls():
    return dict(
        positive_polynomial_comparisons=[concentrated_control(k) for k in (1, 3, 6, 12)],
        invisible_higher_channel=hidden_channel_control(),
        actual_rank_budgets=[actual_fourth_power_budget(s, 8) for s in (F(4), F(9, 2), F(8), F(16), F(32))],
        isolated_SU2_oscillator_constants=oscillator_constants(),
        finite_gauge_cycles=[cycle_control(t) for t in (F(0), F(1, 4), F(-3, 4), F(99, 100))],
        conditional_physical_scaling=[scale_budget(a, 2, F(3, 8), eta, 4, 1)
                                      for a in (F(1, 2), F(1, 16), F(1, 256))
                                      for eta in (a, F(1, 2))])


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'YM84' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM84 source base mismatch')
    if len(manifest['upstream_sha256']) != 361 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM84 source manifest shape mismatch')
    if (HERE / 'EXPECTED_YM83.sha256').read_text().strip() != PARENT_HASH or prev.canonical_sha(json.loads((HERE / 'YM83_RESULT.json').read_text())) != PARENT_HASH:
        raise ValueError('YM83 result/pin mismatch')
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
    return prev.encode(dict(certificate_type='YM84_GAUGE_RESOLVED_CONTINUUM_SCALING',
                            verdict='PASS', runtime_policy='Python 3.12 only', model=MODEL,
                            input_sha256=source_checks(), controls=controls(),
                            claim_status='ACTUAL_OPEN_PATH_TOWER__ISOLATED_PHYSICAL_STIFFNESS__FULL_CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM84_GAUGE_RESOLVED_CONTINUUM_SCALING':
        raise ValueError('YM84 scope mismatch')
    digest = prev.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('YM84 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--write', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = prev.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2) + '\n')
        PIN.write_text(digest + '\n')
    else:
        check(cert)
    print('YM84 PASS', digest)
    print('Actual open-path response and isolated physical stiffness: proved; full 4D continuum gap: open')


if __name__ == '__main__':
    main()
