"""YM80: bounded local phases with certified physical-energy budgets.

Python 3.12 only. Actual finite-regulator bounds use the unchanged SU(2)
vacuum and YM79's smooth sparse source. Gaussian controls test the scope
of the moment hypotheses; they are not Yang--Mills vacua. Default/--check
are read-only; --write changes only this stage's result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym79_local_witness_survival_certificates as previous

native, p = previous.native, previous.p
exact, positive, nonnegative, encode = previous.exact, previous.positive, previous.nonnegative, previous.encode
stage76 = previous.stage76
RESULT = HERE/'YM80_RESULT.json'
PIN = HERE/'EXPECTED_YM80.sha256'
SOURCES = HERE/'YM80_SOURCE_MANIFEST.json'
BASE_TREE = 'cde8043af9528f6d3bc83318b007d3bc18f209e6'
PARENT_HASH = '8f8e603ad24967292d5a6fb134b7fb52ac0b64feba08ee4bc84fa5507490c619'
DEFAULT_ETA, DEFAULT_B = F(7, 30), F(12000)
GRID = 10**18
MODEL = dict(actual_finite_phase_certificate=True,
             unchanged_actual_SU2_ground_state=True,
             physical_phase_energy_budget_certified=True,
             actual_new_bounded_smooth_phase_has_order_a_certificate=True,
             raw_linear_source_order_a_bound_proved=False,
             centered_bounded_phase_norm_at_most_one=True,
             generalized_regularizer_constants_certified=True,
             actual_UV_scale_spectral_packet_certified=True,
             uniform_positive_continuum_survival=False,
             normalization_fixes_UV_escape=False,
             divergent_phase_energy_is_fixed_physical_energy=False,
             OU_sharpness_control_is_actual_Yang_Mills=False,
             separated_time_bound_requires_actual_survival_lower=True,
             regular_local_continuum_Yang_Mills_constructed=False,
             full_Euclidean_covariance=False,
             positive_continuum_mass_gap=False,
             finite_controls_are_universal_proof=False)


def rounded_interval(lower, upper):
    """Outward rounding only; compact exact rationals survive JSON replay."""
    lower, upper = exact(lower), exact(upper)
    if lower > upper:
        raise ValueError('The interval must be ordered')
    return dict(lower=F((lower*GRID).__floor__(), GRID),
                upper=F((upper*GRID).__ceil__(), GRID))


def sqrt_interval(value):
    value = nonnegative(value)
    if not value:
        return dict(lower=F(0), upper=F(0))
    grid = GRID
    while value*grid*grid < 1:
        grid *= 10
    integer = isqrt(value.numerator*grid*grid//value.denominator)
    lower = F(integer, grid)
    return dict(lower=lower, upper=lower if lower*lower == value else F(integer+1, grid))


def exp_negative_interval(value):
    row = stage76.exp_negative_interval(nonnegative(value))
    row = rounded_interval(row['lower'], row['upper'])
    return dict(lower=max(F(0), row['lower']), upper=min(F(1), row['upper']))


@lru_cache(maxsize=None)
def regularizer_budget(eta=DEFAULT_ETA):
    """Analytic constants; use only when theta>=25 and eta**4<=theta.

    The latter condition is exactly epsilon<=1. The parameter eta itself
    need not be <=1. The written proof supplies the derivative estimates;
    this function replays their exact constants, not that universal proof.
    """
    eta = positive(eta)
    inherited = previous.clipped_regularizer_constants()
    K, Z = inherited['commutator_K'], inherited['commutator_Z_upper']
    C = 9*K/eta**2+Z/eta**4+8+24/eta
    inverse_sine = F(183, 5)*eta*eta
    beta = 1/(1+inverse_sine)
    root = sqrt_interval(beta)
    root_upper = F(29, 50) if eta == DEFAULT_ETA else root['upper']
    B = DEFAULT_B if eta == DEFAULT_ETA else sqrt_interval(C)['upper']
    assert root_upper*root_upper >= beta and B*B >= C
    if eta == DEFAULT_ETA:
        assert beta == F(1500, 4489)
        assert C == F(339994914968, 2401) < DEFAULT_B**2
    legacy_beta, legacy_root_upper, legacy_B = F(1280, 1463), F(15, 16), F(2**18)
    coefficient = beta*beta/(root_upper*B)
    legacy_coefficient = legacy_beta**2/(legacy_root_upper*legacy_B)
    return dict(eta=eta, inverse_sine_budget=inverse_sine, beta=beta,
                beta_sqrt_interval=root, beta_sqrt_upper=root_upper,
                commutator_K=K, commutator_Z_upper=Z,
                squared_second_moment_constant_upper=C,
                second_moment_constant=B, squared_constant_majorant=B*B,
                asymptotic_unit_coefficient_lower=coefficient,
                exact_asymptotic_unit_coefficient_interval=dict(
                    lower=beta*root['lower']/B, upper=beta*root['upper']/B),
                legacy_rational_asymptotic_unit_coefficient=legacy_coefficient,
                rational_coefficient_gain_over_legacy=coefficient/legacy_coefficient,
                epsilon_at_most_one_still_required=True,
                analytic_derivative_and_commutator_estimates_required=True,
                finite_arithmetic_proves_universal_derivative_estimates=False)


def phase_chain_moments(phase_squared, raw_energy, raw_second_moment,
                        carre_du_champ_second_moment):
    """Chain identity replay, without inferring a pointwise Gamma bound."""
    phase, q, m2, gamma2 = map(nonnegative, (phase_squared, raw_energy,
                                            raw_second_moment, carre_du_champ_second_moment))
    if gamma2 < q*q:
        raise ValueError('The second Gamma moment cannot be below its mean squared')
    return dict(phase_energy=phase*q,
                phase_second_moment=phase*phase*gamma2+phase*m2,
                real_imaginary_cross_term=F(0), identity_fixture_only=True,
                pointwise_Gamma_bound_inferred=False)


def phase_moment_budget(energy_budget, moment_scale, beta, time=1, beta_sqrt_upper=None):
    """Correlated moment certificate, after minimizing over u in [beta,1].

    The inputs encode Gamma(G)<=q, q_G=u*q, m3(G)<=r**2*q and e=lambda**2*q.
    The denominator at u=beta is not an upper bound on the actual graph
    norm; the *whole moment quotient* increases in u. This distinction
    is essential to the valid substitution of the lower endpoint.
    """
    e, r, beta, time = nonnegative(energy_budget), positive(moment_scale), positive(beta), positive(time)
    if beta > 1:
        raise ValueError('The retained energy fraction must be at most one')
    root_upper = sqrt_interval(beta)['upper'] if beta_sqrt_upper is None else positive(beta_sqrt_upper)
    if root_upper*root_upper < beta:
        raise ValueError('The directed square-root bound points the wrong way')
    denominator = 1+beta*time*e+beta*time*time*e*e+root_upper*r*time*time*e
    lower = beta*beta*time*time*e*e/denominator
    mass_lower = beta*beta*e/(beta*e+root_upper*r)
    assert 0 <= lower <= mass_lower <= 1
    return dict(time=time, energy_budget=e, moment_scale=r, beta=beta,
                beta_sqrt_upper=root_upper, phase_energy_lower=beta*e,
                phase_energy_upper=e, phase_second_moment_upper=e*e+r*e,
                correlated_denominator=denominator,
                reflected_norm_squared_lower=lower,
                centered_norm_squared_lower=mass_lower,
                centered_norm_squared_upper=F(1),
                whole_quotient_monotonicity_used=True,
                denominator_is_actual_graph_norm_upper=False,
                fixed_physical_budget_is_e=True)


def phase_energy_budget(theta, spacing, projected_L2_squared, phase_squared=1,
                        time=1, clock_scale=1, eta=DEFAULT_ETA):
    """Actual sparse smooth-source theorem, for even periodic side L>=4.

    The projected norm is that of YM78's disjoint selected xy cells; G is
    the analytic regularization of their angle sum. The unit-modulus
    phase and its centered vacuum state are in every finite-regulator
    operator domain. No third moment of the raw angle source is used.
    """
    theta, spacing, time, scale = map(positive, (theta, spacing, time, clock_scale))
    projected, phase = nonnegative(projected_L2_squared), nonnegative(phase_squared)
    regularizer = regularizer_budget(eta)
    if theta < 25:
        raise ValueError('The actual-vacuum estimate requires theta >=25')
    if regularizer['eta']**4 > theta:
        raise ValueError('The derivative estimate requires epsilon <=1')
    q = scale*projected
    B = regularizer['second_moment_constant']
    r, e = B*scale/spacing, phase*q
    row = phase_moment_budget(e, r, regularizer['beta'], time, regularizer['beta_sqrt_upper'])
    row.update(theta=theta, spacing=spacing, projected_L2_squared=projected,
               phase_parameter_squared=phase, sqrt_cE_cB=scale,
               regularizer=regularizer,
               regularizer_epsilon_fourth_power=regularizer['eta']**4/theta,
               exact_raw_physical_energy=q,
               smooth_energy_lower=regularizer['beta']*q, smooth_energy_upper=q,
               smooth_second_moment_upper=r*q,
               smooth_third_moment_upper=regularizer['squared_second_moment_constant_upper']*scale*scale*q/spacing**2,
               smooth_third_moment_majorant=r*r*q,
               normalized_a_lower_coefficient=row['reflected_norm_squared_lower']/spacing,
               fixed_energy_a_asymptotic_coefficient_lower=regularizer['asymptotic_unit_coefficient_lower']*e/scale,
               actual_finite_regulator_vacuum=True,
               uniformly_positive_continuum_lower_bound=False)
    return row


def optimized_phase_budget(theta, spacing, projected_L2_squared, physical_energy_budget,
                           time=1, clock_scale=1, eta=DEFAULT_ETA):
    projected, energy, scale = positive(projected_L2_squared), nonnegative(physical_energy_budget), positive(clock_scale)
    q = scale*projected
    row = phase_energy_budget(theta, spacing, projected, energy/q, time, scale, eta)
    row.update(chosen_phase_parameter_squared=energy/q,
               chosen_guaranteed_physical_energy_budget=energy,
               maximizes_certificate_under_guaranteed_energy_budget=True,
               actual_phase_energy_not_assumed_equal_to_budget=True)
    return row


def ultraviolet_packet_budget(moment_scale, beta):
    """Actual phase choice e=r; its packet energies diverge as r=Bc/a.

    Cauchy--Schwarz gives mass above beta*r/2 at least beta**2/8.
    The m2 bound 2*r**2 puts mass above 8*r/beta below beta**2/32.
    """
    r, beta = positive(moment_scale), positive(beta)
    if beta > 1:
        raise ValueError('The retained energy fraction must be at most one')
    return dict(energy_budget=r, phase_energy_lower=beta*r, phase_energy_upper=r,
                phase_second_moment_upper=2*r*r,
                spectral_window=dict(lower=beta*r/2, upper=8*r/beta),
                mass_above_lower_endpoint_lower=beta*beta/8,
                mass_above_upper_endpoint_upper=beta*beta/32,
                spectral_mass_lower=3*beta*beta/32,
                fixed_physical_energy_budget=False,
                packet_is_fixed_physical_energy_window=False,
                positive_continuum_survival_inferred=False)


def separated_time_budget(time, delay, survival_lower, energy_upper):
    """A positive input survival certificate gives a separated-time one.

    Jensen applies to the heat-filtered spectral measure. Its first
    moment is bounded by the unfiltered energy upper bound, not inferred
    from finite-regulator positivity alone.
    """
    time, delay = positive(time), nonnegative(delay)
    delta, energy = nonnegative(survival_lower), nonnegative(energy_upper)
    if delta > time*energy:
        raise ValueError('The proposed survival floor exceeds T times the energy bound')
    exponent = 2*delay*energy/delta if delta else None
    heat = exp_negative_interval(exponent) if delta else dict(lower=F(0), upper=F(0))
    lower = delta*heat['lower']
    return dict(time=time, delay=delay, survival_lower=delta, energy_upper=energy,
                Jensen_exponent_upper=exponent, exponential_interval=heat,
                separated_reflected_norm_squared_lower=lower,
                actual_positive_survival_premise_required=True,
                uniform_positive_input_survival_proved=False)


def _ou_interval_control(rate_lower, rate_upper, phase_energy, time):
    """Outward interval evaluation; rate and Poisson mean may be correlated."""
    low, high, energy, time = map(positive, (rate_lower, rate_upper, phase_energy, time))
    if low > high:
        raise ValueError('The OU rate interval must be ordered')
    mean_lower, mean_upper = energy/high, energy/low
    def correlation(at_time):
        decay_lower = exp_negative_interval(high*at_time)['lower']
        decay_upper = exp_negative_interval(low*at_time)['upper']
        exponent_lower = mean_lower*(1-decay_upper)
        exponent_upper = mean_upper*(1-decay_lower)
        return dict(lower=exp_negative_interval(exponent_upper)['lower'],
                    upper=exp_negative_interval(exponent_lower)['upper'])
    C1, C2 = correlation(time), correlation(2*time)
    atom = dict(lower=exp_negative_interval(mean_upper)['lower'],
                upper=exp_negative_interval(mean_lower)['upper'])
    mass = dict(lower=1-atom['upper'], upper=min(F(1), mean_upper, 1-atom['lower']))
    increment = dict(lower=2*(1-C1['upper']), upper=2*(1-C1['lower']))
    reflected = dict(lower=max(F(0), 1-2*C1['upper']+C2['lower'], increment['lower']**2/4),
                     upper=min(mass['upper'], increment['upper']/2, 1-2*C1['lower']+C2['upper']))
    centered_C1 = dict(lower=max(F(0), C1['lower']-atom['upper']),
                       upper=min(mass['upper'], C1['upper']-atom['lower']))
    assert 0 <= reflected['lower'] <= reflected['upper'] <= 1
    assert 0 <= centered_C1['lower'] <= centered_C1['upper'] <= 1
    return dict(time=time, rate_interval=dict(lower=low, upper=high),
                Poisson_mean_interval=dict(lower=mean_lower, upper=mean_upper),
                phase_energy=energy, phase_centered_mass=mass,
                phase_correlation_at_T=C1, phase_correlation_at_2T=C2,
                phase_centered_correlation_at_T=centered_C1,
                zero_energy_atom=atom, Weyl_increment=increment,
                Weyl_reflected_norm=reflected, OU_fixture_not_Yang_Mills_vacuum=True)


def ou_phase_control(rate, raw_energy=1, time=1, phase_parameter_squared=1):
    rate, q, time, phase = map(positive, (rate, raw_energy, time, phase_parameter_squared))
    energy = phase*q
    row = _ou_interval_control(rate, rate, energy, time)
    row.update(rate=rate, raw_energy=q, phase_parameter_squared=phase,
               Poisson_mean=energy/rate, phase_spectral_lattice_spacing=rate,
               phase_second_moment=energy*energy+rate*energy)
    return row


def ou_sharpness_control(moment_scale, raw_energy=1, energy_budget=1, time=1, beta=1):
    """OU with q_G=beta*q and exact rate r/sqrt(beta), directed throughout.

    Here raw_energy names the reference q in the hypotheses. The actual
    linear-source energy is beta*q. This example demonstrates insufficiency
    of moment data, and never asserts an actual Yang--Mills asymptotic.
    """
    r, q, e, time, beta = map(positive, (moment_scale, raw_energy, energy_budget, time, beta))
    if beta > 1:
        raise ValueError('The retained energy fraction must be at most one')
    root = sqrt_interval(beta)
    low, high = r/root['upper'], r/root['lower']
    row = _ou_interval_control(low, high, beta*e, time)
    row.update(moment_scale=r, reference_raw_energy=q, raw_energy=beta*q,
               raw_second_moment_interval=dict(lower=r*q*root['lower'], upper=r*q*root['upper']),
               raw_third_moment=r*r*q, energy_budget=e, beta=beta,
               phase_parameter_squared=e/q,
               phase_second_moment_interval=dict(lower=(beta*e)**2+low*beta*e,
                                                 upper=(beta*e)**2+high*beta*e),
               asymptotic_mass_scale_interval=dict(lower=beta*root['lower']*e/r,
                                                   upper=beta*root['upper']*e/r),
               available_moment_hypotheses_hold=True,
               proves_actual_Yang_Mills_phase_collapses=False)
    return row


def high_amplitude_control(rate, time=1):
    row = ou_phase_control(rate, raw_energy=1, time=time, phase_parameter_squared=rate)
    limit = exp_negative_interval(1)
    row.update(violates_fixed_energy_budget=True,
               limiting_centered_mass=dict(lower=1-limit['upper'], upper=1-limit['lower']),
               positive_time_centered_correlation_limit=F(0),
               limiting_covariance_continuous_at_zero=False,
               amplitude_renormalization_alone_fixes_UV_escape=False)
    return row


@lru_cache(maxsize=None)
def all_controls():
    regularizer = regularizer_budget()
    assert regularizer['rational_coefficient_gain_over_legacy'] == F(12039575625, 2337530036) > 5
    fixed = [phase_energy_budget(25, F(1, n), 1) for n in (4, 16, 256, 65536)]
    optimized = [optimized_phase_budget(25, F(1, 256), 2, energy)
                 for energy in (F(1, 4), F(1), F(4))]
    OU = [ou_sharpness_control(n, beta=regularizer['beta']) for n in (4, 16, 256, 4096)]
    high = [high_amplitude_control(n) for n in (4, 16, 256)]
    for row in fixed:
        assert row['reflected_norm_squared_lower'] > 0
    for row in OU:
        assert row['phase_centered_mass']['upper'] <= row['asymptotic_mass_scale_interval']['upper']
    return dict(regularizer_constants=regularizer,
                legacy_regularizer_directed_replay=regularizer_budget(F(1, 16)),
                actual_fixed_amplitude_phase_budgets=fixed,
                actual_optimized_energy_phase_budgets=optimized,
                actual_UV_scale_packets=[ultraviolet_packet_budget(row['moment_scale'], regularizer['beta']) for row in fixed],
                phase_chain_identity_controls=[phase_chain_moments(F(9, 64), 4, 8, 48),
                                               phase_chain_moments(F(1, 9), 18, 54, 4860)],
                OU_fixed_energy_sharpness_controls=OU,
                OU_high_amplitude_positive_time_failure=high,
                separated_time_transfers=[separated_time_budget(1, F(1, 4), row['reflected_norm_squared_lower'], row['phase_energy_upper']) for row in fixed],
                conditional_nonvanishing_separated_time_example=separated_time_budget(1, F(1, 4), F(1, 4), 1))


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 295 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM80 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM80 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM79 parent result/pin mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = native.digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = native.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = native.digest(SOURCES)
    return checked


def run():
    native.previous.runtime_check()
    return encode(dict(certificate_type='YM80_BOUNDED_PHASE_ENERGY_CERTIFICATES',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_BOUNDED_PHASE_A_ORDER_CERTIFICATE__UNIFORM_CONTINUUM_SURVIVAL_OPEN',
        proof_scope='Actual local phase energy/moment bounds with directed rational replay and labelled OU controls',
        open='Uniform positive fixed-energy local survival, regular continuum Yang--Mills, Euclidean covariance and mass gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM80_BOUNDED_PHASE_ENERGY_CERTIFICATES':
        raise ValueError('YM80 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM80 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM80 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], regularizer_eta='7/30',
                          physical_moment_constant=12000, uniform_continuum_survival=False)))


if __name__ == '__main__':
    main()
