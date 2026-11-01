"""YM92: actual Wilson flux with a complete nonlinear residual bound.

Exact finite controls accompany analytic full-space bounds, not numerical
Wilson vacua or formal verification. Python 3.12 only, active assertions.
Default/--check are read-only; --write changes only this result and pin.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys

import ym91_conditional_memory_flux as previous
import ym88_spatial_score_contraction as geometry
import ym63_spatial_blocking_scale as blocking

p, a, poly, local = previous.p, previous.a, previous.poly, previous.local
D = geometry.D
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM92_RESULT.json'
PIN = HERE/'EXPECTED_YM92.sha256'
SOURCES = HERE/'YM92_SOURCE_MANIFEST.json'
BASE_COMMIT = 'bec73061e92b56f4772d5a31063dde4da0bcaa65'
BASE_TREE = '81db1458bed2ca1b15dd5b0309b7f2ba5ea8a2e3'
NSB2_HASH = 'f60024412988385d8febe860117697002d7fa64d02a0d173d0cf07d2b0965091'
THETA_MAX = F(1, 400)
MODEL = dict(
    actual_Wilson_covariant_response_remainder_bound=True,
    actual_full_and_conditional_gap_floor_improved=True,
    actual_conditional_flux_corrected_with_entire_nonlinear_residual=True,
    actual_dynamic_memory_bound_uniform_in_volume_and_block_factor=True,
    original_Wilson_Hamiltonian_modified=False,
    actual_vacuum_replaced_by_first_order_or_Gibbs_ansatz=False,
    nonlinear_residual_discarded=False,
    Haar_flux_alone_is_actual_weighted_flux=False,
    scalar_density_curvature_spectra_supply_missing_source_data=False,
    comparison_polynomials_are_numerical_Wilson_vacua=False,
    coupling_window='0 <= theta <= 1/400, unchanged',
    endpoint_bare_gap_lower='251/800',
    endpoint_relative_memory_upper='2105401/3189425625',
    weak_coupling_window_extended=False,
    nonzero_continuum_states_constructed=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def theta_value(theta):
    theta = p.nonnegative(theta)
    if theta > THETA_MAX:
        raise ValueError('Outside the proved window 0 <= theta <= 1/400')
    return theta


def remainder_budget(theta):
    t = theta_value(theta)
    preliminary = 918*t*t
    denominator = 1-144*t-2*preliminary
    refined = 310*t*t/denominator
    assert denominator >= F(25141, 40000)
    assert refined <= F(12400000, 25141)*t*t <= 500*t*t
    return dict(theta=t, preliminary_mixed_remainder=preliminary,
                bootstrap_denominator=denominator, refined_mixed_remainder=refined,
                mixed_remainder_upper=500*t*t,
                bootstrap_uses_existing_global_bound=True,
                quadratic_small_root_selected_without_input=False)


def actual_budget(theta, factor=2, damping=0):
    t, b, s = theta_value(theta), p.size(factor, 2), p.nonnegative(damping)
    eps, G, delta = 500*t*t, 42*t, 35*t
    beta = t+eps
    hrow = delta+beta
    d = F(1, 2)-2*hrow
    j = 8*t/9
    residual = 2*eps+2*hrow*j+F(8, 3)*t*G
    amplitude = j+residual/d
    assert amplitude == (F(4, 9)*t+1112*t*t)/d
    eta = amplitude*amplitude
    sigma2 = 4*beta*beta/d
    alpha = min(eta, sigma2/(d+s))
    assert d >= F(251, 800) and eta < F(1, 1500)
    assert eta <= 4*beta*beta/(d*d)
    return dict(theta=t, factor=b, damping=s, mixed_remainder_upper=eps,
                mixed_Hessian_row_upper=beta, full_Hessian_row_upper=hrow,
                bare_full_and_conditional_gap_lower=d,
                trial_flux_norm_over_sqrt_b=j,
                residual_gradient_over_sqrt_b=residual,
                residual_variance_over_b=residual*residual/d,
                correction_flux_energy_over_b=residual*residual/(d*d),
                completed_flux_norm_over_sqrt_b=amplitude,
                conditional_tensor_norm_upper=b*eta, relative_memory_upper=eta,
                source_energy_factor=sigma2, frequency_memory_upper=alpha,
                relative_resolvent_error=alpha/(1-alpha))


def endpoint_control():
    out = actual_budget(THETA_MAX)
    eta = out['relative_memory_upper']
    assert out['bare_full_and_conditional_gap_lower'] == F(251, 800)
    assert out['completed_flux_norm_over_sqrt_b'] == F(1451, 56475)
    assert eta == F(2105401, 3189425625)
    assert out['relative_resolvent_error'] == F(2105401, 3187320224)
    assert F(1, 144)/eta > F(21, 2)
    return dict(actual=out, old_gap_lower=F(3, 10), old_memory_upper=F(1, 144),
                old_relative_resolvent_error=F(1, 143),
                memory_upper_improvement_factor=F(1, 144)/eta,
                gap_floor_increase=F(11, 800),
                comparison_is_of_upper_certificates_not_actual_memory_values=True)


def covariant_block(tensor, e, f, direction):
    """Covariant derivative of a 3x3 T_e* tensor T_f* block."""
    return [[poly.add(D(tensor[i][j], direction),
                      *(poly.scale(tensor[k][j], -geometry.connection(direction, 3*e+i, 3*e+k))
                        for k in range(3)),
                      *(poly.scale(tensor[i][k], -geometry.connection(direction, 3*f+j, 3*f+k))
                        for k in range(3))) for j in range(3)] for i in range(3)]


@lru_cache(maxsize=None)
def plaquette_mixed(e, f):
    if type(e) is not int or type(f) is not int or e == f or not 0 <= min(e, f) <= max(e, f) < 4:
        raise ValueError('Distinct links of the four-link plaquette required')
    W = local.model.loop_polynomial()
    return [[D(D(W, 3*f+j), 3*e+i) for j in range(3)] for i in range(3)]


@lru_cache(maxsize=None)
def linear_resolvent_control(e=0, f=1):
    H = plaquette_mixed(e, f)
    rough = [[{} for _ in range(3)] for _ in range(3)]
    for k in range(12):
        second = covariant_block(covariant_block(H, e, f, k), e, f, k)
        for i in range(3):
            for j in range(3):
                rough[i][j] = poly.add(rough[i][j], second[i][j])
    wrong = 0
    for i in range(3):
        for j in range(3):
            assert poly.add(poly.scale(rough[i][j], F(-1, 3)), poly.scale(H[i][j], F(1, 3))) == H[i][j]
            raw = poly.add(*(D(D(H[i][j], k), k) for k in range(12)))
            wrong += raw != rough[i][j]
    assert wrong == 9
    return dict(links=[e, f], native_mixed_components_checked=9,
                leading_log_wave_coefficient=F(1, 3), plaquette_Casimir=F(3),
                mixed_tensor_Ricci_damping=F(1),
                components_changed_by_scalar_instead_of_rough_Laplacian=wrong)


@lru_cache(maxsize=None)
def third_derivative_control():
    W = local.model.loop_polynomial()
    qs = ((F(3, 5), F(4, 5), 0, 0), (0, F(5, 13), F(12, 13), 0),
          (F(1, 2),)*4, (0, 0, 0, 1))
    repeated, distinct, changed = 0, 0, 0
    for e, f in combinations(range(4), 2):
        H = plaquette_mixed(e, f)
        for g in range(4):
            for c in range(3):
                T = covariant_block(H, e, f, 3*g+c)
                for i in range(3):
                    for j in range(3):
                        if g == e:
                            expected = poly.scale(D(W, 3*f+j), -F(i == c, 4))
                            repeated += 1
                        elif g == f:
                            expected = poly.scale(D(W, 3*e+i), -F(j == c, 4))
                            repeated += 1
                        else:
                            expected = D(H[i][j], 3*g+c)
                            distinct += 1
                        assert T[i][j] == expected
                        changed += T[i][j] != D(H[i][j], 3*g+c)
                        assert abs(local.model.evaluate(T[i][j], qs)) <= F(1, 8)
    return dict(repeated_factor_components=repeated, distinct_factor_components=distinct,
                components_changed_by_omitting_connections=changed,
                covariant_third_directional_bound=F(1, 8),
                finite_samples_alone_prove_universal_norm_bound=False)


@lru_cache(maxsize=None)
def flux_polynomials(retained=(0, 1)):
    if not isinstance(retained, tuple) or not 1 <= len(retained) <= 2 or len(set(retained)) != len(retained):
        raise ValueError('One or two distinct retained plaquette links required')
    if any(type(e) is not int or not 0 <= e < 4 for e in retained):
        raise ValueError('Invalid retained plaquette link')
    directions = ((F(3, 5), F(4, 5), 0), (0, F(5, 13), F(12, 13)))
    W = local.model.loop_polynomial()
    source = poly.add(*(poly.scale(local.derivative(W, e, directions[i]), F(2, 3))
                        for i, e in enumerate(retained)))
    hidden = tuple(e for e in range(4) if e not in retained)
    eigenvalue = F(3*len(hidden), 4)
    phi = poly.scale(source, 1/eigenvalue)
    J = [D(phi, 3*e+i) for e in hidden for i in range(3)]
    divergence = poly.add(*(D(J[3*j+i], 3*e+i) for j, e in enumerate(hidden) for i in range(3)))
    assert poly.scale(divergence, -1) == source
    assert poly.scale(poly.add(*(D(D(source, 3*e+i), 3*e+i)
                                for e in hidden for i in range(3))), -1) == poly.scale(source, eigenvalue)
    averaged = source
    for e in hidden:
        averaged = local.integrate_link(averaged, e)
    assert not averaged
    return dict(retained=retained, hidden=hidden, directions=directions,
                source=source, potential=phi, flux=J, eigenvalue=eigenvalue)


def flux_control(retained):
    row = flux_polynomials(retained)
    return dict(retained=row['retained'], hidden=row['hidden'], conditional_free_Casimir=row['eigenvalue'],
                source_sha256=p.canonical_sha(p.encode(row['source'])),
                trial_potential_sha256=p.canonical_sha(p.encode(row['potential'])),
                exact_Haar_divergence=True, conditional_Haar_mean_zero=True,
                Haar_flux_alone_is_actual_weighted_flux=False)


@lru_cache(maxsize=None)
def weighted_residual_control():
    row = flux_polynomials((0, 1))
    W = local.model.loop_polynomial()
    # h=exp(2 W/5) is an engineered native parent, used only for derivative signs.
    chi = poly.scale(W, F(2, 5))
    gradient = [D(chi, i) for i in range(12)]
    H = geometry.hessian(chi, 4)
    raw = poly.add(*(poly.scale(local.derivative(chi, e, row['directions'][i]), 2)
                     for i, e in enumerate(row['retained'])))
    indices = [3*e+i for e in row['hidden'] for i in range(3)]
    J = row['flux']
    drift = poly.add(*(poly.scale(poly.mul(gradient[k], J[j]), 2) for j, k in enumerate(indices)))
    residual = poly.add(raw, poly.scale(row['source'], -1), drift)
    changed = 0
    for k in indices:
        first = D(poly.add(raw, poly.scale(row['source'], -1)), k)
        hess_term = poly.add(*(poly.scale(poly.mul(H[k][ell], J[j]), 2)
                              for j, ell in enumerate(indices)))
        covJ = [poly.add(D(J[j], k), *(poly.scale(J[v], -geometry.connection(k, ell, z))
                                     for v, z in enumerate(indices))) for j, ell in enumerate(indices)]
        third_term = poly.add(*(poly.scale(poly.mul(gradient[ell], covJ[j]), 2)
                               for j, ell in enumerate(indices)))
        assert D(residual, k) == poly.add(first, hess_term, third_term)
        changed += D(residual, k) != first
    assert changed == len(indices)
    return dict(kind='engineered_polynomial_log_density_not_Wilson_vacuum',
                conditional_gradient_components_checked=len(indices),
                components_wrong_if_weighted_drift_is_dropped=changed,
                residual_gradient_sha256=p.canonical_sha(p.encode([D(residual, k) for k in indices])),
                actual_conditional_center_is_a_Z_constant=True,
                numerical_conditional_exponential_integral_claimed=False)


@lru_cache(maxsize=None)
def incidence_control(factor):
    p.size(factor, 2)
    cell = blocking.blocked_geometry(4, factor)
    retained = {e for path in cell['paths'] for e in path}
    flux_rows, flux_cols, drift_rows, drift_cols = Counter(), Counter(), Counter(), Counter()
    histogram = Counter()
    for face in cell['fine']['faces']:
        r = [e for e, _ in face if e in retained]
        z = [e for e, _ in face if e not in retained]
        histogram[len(r)] += 1
        assert len(r) <= 2
        if not r:
            continue
        for e in r:
            for w in z:
                flux_cols[e] += F(2, 9*len(z))
                flux_rows[w] += F(2, 9*len(z))
                drift_cols[e] += F(1, 9)
                drift_rows[w] += F(1, 9)
    upper = dict(flux_row=max(flux_rows.values()), flux_column=max(flux_cols.values()),
                 drift_row=max(drift_rows.values()), drift_column=max(drift_cols.values()))
    assert upper['flux_row'] <= F(8, 9) and upper['flux_column'] <= F(8, 9)
    assert upper['drift_row'] <= F(8, 9) and upper['drift_column'] <= F(4, 3)
    return dict(factor=factor, fine_links=len(cell['fine']['edges']), retained_links=len(retained),
                coarse_links=len(cell['paths']), retained_face_histogram=dict(histogram),
                normalized_incidence_bounds=upper, copy_norm_squared=factor)


def approximate_flux_certificate(gradient, source, trial_flux, floor, young=1):
    N = a.matrix(gradient)
    R = a.matrix(source, len(N[0]))
    J = a.matrix(trial_flux, len(N), len(R[0]))
    d, t = p.positive(floor), p.positive(young)
    A = p.multiply(p.transpose(N), N)
    if not p.psd(p.add(A, a.scale(p.identity(len(A)), d), -1)):
        raise ValueError('The claimed centered conditional floor is invalid')
    residual = p.add(R, p.multiply(p.transpose(N), J), -1)
    correction = p.multiply(N, p.multiply(a.inverse(A), residual))
    complete = p.add(J, correction)
    assert p.multiply(p.transpose(N), complete) == R
    gram = lambda x:p.multiply(p.transpose(x), x)
    exact = p.multiply(p.transpose(R), p.multiply(a.inverse(A), R))
    residual_energy = p.multiply(p.transpose(residual), p.multiply(A, residual))
    upper = p.add(a.scale(gram(J), 1+t), a.scale(residual_energy, (1+1/t)/(d*d)))
    assert p.psd(p.add(gram(complete), exact, -1))
    assert p.psd(p.add(upper, gram(complete), -1))
    return dict(kind='finite_centered_operator_control_not_Wilson_vacuum',
                conditional_operator=A, source=R, trial_flux=J, residual=residual,
                correction_flux=correction, complete_flux=complete,
                trial_Gram=gram(J), corrected_Gram=gram(complete), residual_energy_Gram=residual_energy,
                exact_conditional_tensor=exact, full_residual_upper=upper)


def finite_flux_control():
    row = previous.matrix_control()
    N, R, A = row['gradient'], row['scores'], row['conditional_operator']
    optimum = p.multiply(N, p.multiply(a.inverse(A), R))
    trial = a.scale(optimum, F(2, 3))
    out = approximate_flux_certificate(N, R, trial, 1)
    assert not p.psd(p.add(out['trial_Gram'], out['exact_conditional_tensor'], -1))
    exact = approximate_flux_certificate([[1, 0], [0, 1]], [[1, 1], [0, 1]], [[1, 1], [0, 1]], 1)
    tensor = exact['exact_conditional_tensor']
    diagonal = [[tensor[i][i] if i == j else F(0) for j in range(2)] for i in range(2)]
    assert not p.psd(p.add(diagonal, tensor, -1))
    return dict(residual_required=out, cross_terms_required=exact)


def controls():
    return dict(endpoint=endpoint_control(),
                remainder_budgets=[remainder_budget(t) for t in (0, F(1, 8000), F(1, 1000), THETA_MAX)],
                actual_budgets=[actual_budget(t, b, s) for t in (0, F(1, 1000), THETA_MAX)
                                for b, s in ((2, 0), (3, F(1, 2)), (8, 10))],
                native_linear_resolvent=[linear_resolvent_control(0, 1), linear_resolvent_control(0, 2)],
                native_third_derivatives=third_derivative_control(),
                native_plaquette_fluxes=[flux_control(r) for r in ((0,), (0, 1), (0, 3))],
                weighted_residual=weighted_residual_control(),
                geometric_incidence=[incidence_control(b) for b in (2, 3, 4)],
                finite_fluxes=finite_flux_control(),
                inherited_frequency_counterexample=previous.previous.insufficient_zero_frequency_control())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'YM92' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM92 source base mismatch')
    if len(manifest['upstream_sha256']) != 427 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM92 source manifest shape mismatch')
    nsb = ROOT/'papers/nonlinear-spectral-blindness/certificates'
    if (nsb/'EXPECTED_NSB2.sha256').read_text().strip() != NSB2_HASH or p.canonical_sha(json.loads((nsb/'NSB2_RESULT.json').read_text())) != NSB2_HASH:
        raise ValueError('NSB2 parent evidence mismatch')
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
    return p.encode(dict(certificate_type='YM92_WILSON_FLUX_CORRECTOR', verdict='PASS',
                         runtime_policy='Python 3.12 only', model=MODEL,
                         input_sha256=source_checks(), controls=controls(),
                         claim_status='ACTUAL_WILSON_FLUX_AND_GAP_IMPROVEMENT__WEAK_COUPLING_AND_CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM92_WILSON_FLUX_CORRECTOR':
        raise ValueError('YM92 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('YM92 fresh result/pin mismatch')
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
    print('YM92 PASS', digest)
    print('Actual Wilson memory < 1/1500 and bare gap >= 251/800 at theta=1/400; 4D continuum open')


if __name__ == '__main__':
    main()
