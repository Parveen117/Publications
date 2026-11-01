"""YM93: local Bochner control and a fourfold Wilson continuation window.

Analytic proofs concern the actual full Wilson vacuum. Exact polynomial
controls concern comparison potentials, not computed Wilson ground states.
Python 3.12, stdlib rational arithmetic; default/--check are read-only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from math import factorial
from pathlib import Path
import sys

import ym92_wilson_flux_corrector as previous

p, poly, geometry, local = previous.p, previous.poly, previous.geometry, previous.local
D = geometry.D
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM93_RESULT.json'
PIN = HERE/'EXPECTED_YM93.sha256'
SOURCES = HERE/'YM93_SOURCE_MANIFEST.json'
BASE_COMMIT = '15f9062ec73d563b5fcb33d50135ec5b2f14994c'
BASE_TREE = 'a7a37917b2e0a201a79c97e07d604392137d81d4'
PARENT_HASH = '73d8a5c2a3382f055c03ecc85eced3e6fabf9ccbb889d1bb702d98573b0f067b'
THETA_MAX = F(1, 100)
BARRIER = F(13, 500)
MODEL = dict(
    actual_Wilson_local_log_gradient_upper='4 theta, every finite theta >= 0',
    actual_Wilson_gap_window='0 <= theta <= 1/100',
    previous_gap_window='0 <= theta <= 1/400',
    coupling_window_enlargement_factor=4,
    actual_diagonal_covariant_Hessian_upper='15 theta',
    actual_mixed_first_response_remainder_upper='260 theta^2',
    actual_gap_lower='1/2 - 32 theta - 520 theta^2',
    endpoint_bare_gap_lower='16/125',
    endpoint_relative_memory_upper='418609/2073600',
    continuation_from_the_free_vacuum_selects_small_branch=True,
    all_nonlinear_residual_modes_included=True,
    uniform_in_finite_volume_and_block_factor=True,
    older_numeric_helpers_used_outside_their_windows=False,
    original_Wilson_Hamiltonian_modified=False,
    comparison_polynomials_are_computed_Wilson_vacua=False,
    all_coupling_gradient_bound_implies_all_coupling_gap=False,
    weak_coupling_continuum_trajectory_controlled=False,
    nonzero_continuum_states_constructed=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def theta_value(theta):
    theta = p.nonnegative(theta)
    if theta > THETA_MAX:
        raise ValueError('Outside the new proved gap window 0 <= theta <= 1/100')
    return theta


def gradient_upper(theta):
    """All finite nonnegative theta; this alone is not a gap certificate."""
    return 4*p.nonnegative(theta)


def exp_lower(x, degree):
    x, n = p.nonnegative(x), p.size(degree, 0)
    return sum((x**k/F(factorial(k)) for k in range(n+1)), F(0))


def exp_upper(x):
    """Positive Taylor series plus its full geometric tail, for 0 <= x < 3."""
    x = p.nonnegative(x)
    if x >= 3:
        raise ValueError('The exponential enclosure requires x < 3')
    return 1+x+x*x/(2*(1-x/3))


def squared_geometric_tail(q, start):
    q, n = p.nonnegative(q), p.size(start, 0)
    if q >= 1:
        raise ValueError('A convergent geometric tail requires q < 1')
    return q**n*((n+1)**2/(1-q)+2*(n+1)*q/(1-q)**2+q*(1+q)/(1-q)**3)


def heat_control():
    # At native time four: Casimir n(n+2)/4, dimension n+1.
    # n=1 uses exp(3)>20; n>=2 uses n(n+2)>=4n and exp(4)>50.
    lower3, lower4 = exp_lower(3, 8), exp_lower(4, 7)
    assert lower3 > 20 and lower4 > 50
    tail = squared_geometric_tail(F(1, 50), 2)
    total = F(1, 5)+tail
    assert total == F(599222, 2941225) < F(21, 100)
    return dict(native_heat_time=4, fundamental_Casimir=F(3, 4),
                exp3_lower=lower3, exp4_lower=lower4,
                channels_n_at_least_2_full_tail=tail,
                complete_kernel_deviation_upper=total,
                chosen_kernel_deviation=F(21, 100),
                infinite_representation_tail_included=True)


def diagonal_budget(theta):
    t = theta_value(theta)
    M = gradient_upper(t)
    X = exp_upper(16*t)
    # Diameter=2*pi; h/g in [exp(-2*pi*M), exp(2*pi*M)].
    # Therefore r^2 <= exp(4*pi*M) <= exp((352/7)*theta).
    R2 = exp_upper(F(352, 7)*t)
    denominator = 2-F(121, 100)*X
    assert denominator >= F(102957, 177500) > 0
    second = R2*(X-1)*(F(1, 4)+M)/denominator
    delta = second+M*M
    assert delta <= 15*t
    return dict(theta=t, local_potential_size=M, log_gradient_upper=M,
                heat_exp_upper=X, square_Harnack_ratio_upper=R2,
                Duhamel_contraction_upper=1-denominator,
                Duhamel_denominator=denominator,
                relative_second_wave_derivative_upper=second,
                adaptive_diagonal_Hessian_upper=delta,
                diagonal_Hessian_upper=15*t,
                positivity_of_an_interacting_gap_assumed=False)


def barrier_polynomial(theta, error):
    t, e = theta_value(theta), p.nonnegative(error)
    return 78*t*t-(1-64*t)*e+2*e*e


def remainder_budget(theta):
    t = theta_value(theta)
    barrier_value = barrier_polynomial(t, BARRIER)
    assert barrier_value <= -F(13, 62500) < 0
    denominator = 1-64*t-2*BARRIER
    refined = 78*t*t/denominator
    assert denominator >= F(77, 250)
    assert refined <= F(19500, 77)*t*t <= 260*t*t
    return dict(theta=t, strict_continuation_barrier=BARRIER,
                polynomial_at_barrier=barrier_value,
                refinement_denominator=denominator,
                refined_remainder_upper=refined,
                mixed_remainder_upper=260*t*t,
                requires_finite_volume_continuity_from_zero=True,
                quadratic_inequality_alone_selects_small_branch=False)


def actual_budget(theta, factor=2, damping=0):
    t, b, s = theta_value(theta), p.size(factor, 2), p.nonnegative(damping)
    G, delta, eps = gradient_upper(t), 15*t, 260*t*t
    beta, hrow = t+eps, delta+t+eps
    d = F(1, 2)-2*hrow
    j = 8*t/9
    residual = 2*eps+2*hrow*j+F(8, 3)*t*G
    amplitude = j+residual/d
    assert amplitude == (F(4, 9)*t+F(1592, 3)*t*t)/d
    eta = amplitude*amplitude
    sigma2 = 4*beta*beta/d
    alpha = min(eta, sigma2/(d+s))
    assert d >= F(16, 125) and eta <= F(418609, 2073600) < F(1, 4)
    assert eta <= 4*beta*beta/(d*d)
    return dict(theta=t, factor=b, damping=s, log_gradient_upper=G,
                diagonal_Hessian_upper=delta, mixed_remainder_upper=eps,
                mixed_Hessian_row_upper=beta, full_Hessian_row_upper=hrow,
                bare_full_conditional_and_coarse_gap_lower=d,
                trial_flux_norm_over_sqrt_b=j,
                residual_gradient_over_sqrt_b=residual,
                completed_flux_norm_over_sqrt_b=amplitude,
                conditional_tensor_norm_upper=b*eta, relative_memory_upper=eta,
                source_energy_factor=sigma2, frequency_memory_upper=alpha,
                relative_resolvent_error=alpha/(1-alpha))


def endpoint_control():
    end = actual_budget(THETA_MAX)
    assert end['completed_flux_norm_over_sqrt_b'] == F(647, 1440)
    assert end['relative_resolvent_error'] == F(418609, 1654991)
    old_t = F(1, 400)
    # This comparison calls the parent's helper only inside its own domain.
    old = previous.actual_budget(old_t)
    now = actual_budget(old_t)
    improvement = old['relative_memory_upper']/now['relative_memory_upper']
    assert now['bare_full_conditional_and_coarse_gap_lower'] == F(1667, 4000)
    assert now['relative_memory_upper'] == F(635209, 5627250225) < F(1, 8800)
    assert improvement > F(29, 5)
    return dict(new_endpoint=end, old_endpoint_with_new_bounds=now,
                old_endpoint_YM92=old, same_theta_memory_upper_improvement=improvement,
                endpoint_window_factor=THETA_MAX/old_t,
                comparison_is_of_certificates_not_computed_memory=True)


def comparison_chi(links):
    if type(links) is not int or links not in (1, 2, 3):
        raise ValueError('Comparison controls support one, two or three links')
    x = [geometry.variable(i, links) for i in range(4*links)]
    chi = poly.add(poly.scale(x[0], F(2, 3)), poly.scale(poly.mul(x[1], x[1]), F(1, 7)))
    for k in range(1, links):
        chi = poly.add(chi, poly.scale(poly.mul(x[0], x[4*k+1]), F(3, k+8)))
    return chi


def covariant_one_form(vector, direction):
    n = len(vector)
    return [poly.add(D(vector[i], direction),
                     *(poly.scale(vector[v], -geometry.connection(direction, i, v))
                       for v in range(n))) for i in range(n)]


@lru_cache(maxsize=None)
def vector_control(links):
    chi = comparison_chi(links)
    n = 3*links
    U = [D(chi, i) for i in range(n)]
    V = poly.scale(poly.add(*(D(U[k], k) for k in range(n)),
                            *(poly.mul(u, u) for u in U)), -1)
    DU = [covariant_one_form(U, k) for k in range(n)]
    DDU = [covariant_one_form(DU[k], k) for k in range(n)]
    rough = [poly.add(*(DDU[k][i] for k in range(n))) for i in range(n)]
    drift = [poly.add(*(poly.mul(U[k], DU[k][i]) for k in range(n))) for i in range(n)]
    wrong_scalar, wrong_Ricci, wrong_drift = 0, 0, 0
    for i in range(n):
        lhs = poly.add(poly.scale(rough[i], -1), poly.scale(drift[i], -2), poly.scale(U[i], F(1, 2)))
        assert lhs == D(V, i)
        wrong_scalar += rough[i] != poly.add(*(D(D(U[i], k), k) for k in range(n)))
        wrong_Ricci += bool(U[i])
        wrong_drift += bool(drift[i])
    assert wrong_scalar == wrong_Ricci == wrong_drift == n
    # Independently evaluate the scalar partial Bochner identity on each factor.
    for e in range(links):
        inds = range(3*e, 3*e+3)
        Q = poly.add(*(poly.mul(U[i], U[i]) for i in inds))
        LQ = poly.add(*(D(D(Q, k), k) for k in range(n)),
                      *(poly.scale(poly.mul(U[k], D(Q, k)), 2) for k in range(n)))
        rhs = poly.add(*(poly.mul(DU[k][i], DU[k][i]) for k in range(n) for i in inds),
                       poly.scale(Q, F(1, 2)),
                       *(poly.scale(poly.mul(U[i], D(V, i)), -1) for i in inds))
        assert poly.scale(LQ, F(1, 2)) == rhs
    return dict(kind='engineered_positive_ground_state_exp_chi_not_Wilson_vacuum',
                links=links, chi_sha256=p.canonical_sha(p.encode(chi)),
                potential_sha256=p.canonical_sha(p.encode(V)),
                covariant_one_form_components=n, scalar_partial_Bochner_identities=links,
                components_changed_by_scalar_Laplacian=wrong_scalar,
                components_changed_by_dropping_Ricci=wrong_Ricci,
                components_changed_by_dropping_drift=wrong_drift,
                native_Ricci_damping=F(1, 2))


@lru_cache(maxsize=None)
def derivative_control():
    W = local.model.loop_polynomial()
    nvars = 16
    one = {(0,)*nvars:F(1)}
    wave = poly.add(one, poly.scale(W, F(2, 5)))
    directions = ((F(3, 5), F(4, 5), 0), (0, F(5, 13), F(12, 13)))
    checked = 0
    for e in range(4):
        for u in directions:
            deriv = lambda f:local.derivative(f, e, u)
            second = lambda f:deriv(deriv(f))
            assert second(W) == poly.scale(W, -F(1, 4))
            assert local.model.electric(second(wave), 4) == second(local.model.electric(wave, 4))
            lhs = poly.add(second(poly.mul(W, wave)), poly.scale(poly.mul(W, second(wave)), -1))
            rhs = poly.add(poly.mul(second(W), wave), poly.scale(poly.mul(deriv(W), deriv(wave)), 2))
            assert lhs == rhs
            checked += 1
    return dict(unit_direction_cases=checked, second_plaquette_derivative_factor=F(-1, 4),
                electric_commutators_and_full_product_rules=True,
                test_wave_is_actual_Wilson_ground_state=False)


def controls():
    return dict(native_vector=[vector_control(n) for n in (1, 2, 3)],
                second_derivative_equation=derivative_control(), heat_kernel=heat_control(),
                endpoint=endpoint_control(),
                diagonal_budgets=[diagonal_budget(t) for t in (0, F(1, 400), F(1, 200), THETA_MAX)],
                continuation_budgets=[remainder_budget(t) for t in (0, F(1, 400), F(1, 200), THETA_MAX)],
                actual_budgets=[actual_budget(t, b, s) for t in (0, F(1, 400), THETA_MAX)
                                for b, s in ((2, 0), (3, F(1, 2)), (8, 10))],
                all_coupling_gradient_budgets={t:gradient_upper(t) for t in (0, F(1, 100), 1, 100)},
                inherited_covariant_plaquette_response=previous.linear_resolvent_control(),
                inherited_native_fluxes=[previous.flux_control(r) for r in ((0,), (0, 1))],
                inherited_full_residual_Gram_control=previous.finite_flux_control(),
                inherited_frequency_counterexample=previous.previous.previous.insufficient_zero_frequency_control())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'YM93' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM93 source base mismatch')
    if len(manifest['upstream_sha256']) != 434 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM93 source manifest shape mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM92 parent evidence mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = sha256((ROOT/path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed upstream input: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return p.encode(dict(certificate_type='YM93_BOCHNER_CONTINUATION', verdict='PASS',
                         runtime_policy='Python 3.12 only', model=MODEL,
                         input_sha256=source_checks(), controls=controls(),
                         claim_status='ACTUAL_WILSON_FOURFOLD_WINDOW__WEAK_COUPLING_CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM93_BOCHNER_CONTINUATION':
        raise ValueError('YM93 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('YM93 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(digest+'\n')
    else:
        check(cert)
    print('YM93 PASS', digest)
    print('Actual Wilson gap >= 16/125 and memory < 1/4 through theta=1/100; 4D continuum open')


if __name__ == '__main__':
    main()
