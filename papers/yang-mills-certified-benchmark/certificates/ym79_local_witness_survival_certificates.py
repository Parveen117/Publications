"""YM79: actual local resolvent witnesses and bounded-phase survival bridges.

Python 3.12 only. The four-cube witness uses the actual untruncated vacuum
with certified uncertainty. Matrix and Ornstein--Uhlenbeck fixtures are
labelled comparisons, not interacting Yang--Mills vacua.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym78_stationary_derivative_energy_survival as previous

stage77, stage76, stage75 = previous.previous, previous.stage76, previous.stage75
stage74 = stage75.previous
native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative, encode = previous.positive, previous.nonnegative, previous.encode
model = previous.model
RESULT = HERE/'YM79_RESULT.json'
PIN = HERE/'EXPECTED_YM79.sha256'
SOURCES = HERE/'YM79_SOURCE_MANIFEST.json'
BASE_TREE = 'dcf8620a674b12edcba3a924b3f700ba87a17bcb'
PARENT_HASH = 'e12e0b51caea0528c63da7da99b2bb4a1bb8af7c2cead207a6724a90e1fa9f3e'
MODEL = dict(actual_four_cube_local_Wilson_resolvent_witness=True,
             actual_vacuum_error_certified_not_set_to_trial=True,
             actual_untruncated_neutral_hidden_space=True,
             all_coupling_fixed_regulator_source_residual_certificate=True,
             fixed_regulator_refinement_implies_uniform_continuum_survival=False,
             actual_bounded_Weyl_sources_exclude_IR_energy_mass=True,
             bounded_phase_reflected_norm_is_at_most_one=True,
             ordinary_increment_lower_bound_implies_bounded_phase_survival=True,
             actual_uniform_positive_increment_lower_bound_proved=False,
             actual_regularized_witness_improves_a4_to_a2_lower_bound=True,
             regularized_a2_lower_bound_is_uniformly_positive_as_a_to_zero=False,
             OU_IR_phase_rescue_is_actual_Yang_Mills_vacuum=False,
             survival_of_one_phase_state_implies_gap_for_all_rational_phases=False,
             one_time_Weyl_field_regularity_or_local_YM_identification_proved=False,
             full_Euclidean_covariance=False, positive_continuum_mass_gap=False,
             finite_controls_are_universal_proof=False)


def sqrt_interval(value):
    return stage74.sqrt_interval(value)


def norm_interval(vector):
    vector = tuple(map(exact, vector))
    return sqrt_interval(sum((x*x for x in vector), F(0)))


def vacuum_distance_budget(rayleigh, ground_lower, first_excited_lower):
    """Distance to the ground-state line; the chosen real sign needs certification."""
    rho, lower, excited = map(exact, (rayleigh, ground_lower, first_excited_lower))
    if not lower <= rho < excited:
        raise ValueError('The Rayleigh value must be between a ground floor and a strict excited floor')
    leakage = (rho-lower)/(excited-lower)
    overlap = sqrt_interval(1-leakage)
    distance = sqrt_interval(2*(1-overlap['lower']))['upper']
    return dict(rayleigh=rho, ground_energy_lower=lower, first_excited_energy_lower=excited,
                excited_weight_upper=leakage, ground_overlap_modulus_lower=overlap['lower'],
                distance_upper_after_certified_phase_alignment=distance,
                phase_alignment_must_be_verified=True)


@lru_cache(maxsize=None)
def native_source_moments():
    """Native one-holonomy Haar moments plus actual cube parity/private-edge checks."""
    graph = model.lattice(4)
    faces = [frozenset(edge for edge, _ in face) for face in graph['faces']]
    assert len(faces) == len(set(faces)) == 192
    selected = faces[0]
    private_edges = [len(face-selected) for face in faces[1:]]
    assert min(private_edges) >= 1
    # In V^2 W^(2k), distinct plaquettes leave an odd link. W^(2k) is even.
    pair_xors = {left ^ right for left, right in combinations(faces, 2)}
    assert all(pair_xors) and not pair_xors.intersection(faces)
    # The last check also eliminates the three-plaquette terms in <V^3>.
    coordinate, power = p.COORD[0], {(0, 0, 0, 0): F(1)}
    moments = {}
    for degree in range(1, 7):
        power = p.mul(power, coordinate)
        moments[degree] = native.phi(power, 1)
    assert [moments[k] for k in (1, 3, 5)] == [0, 0, 0]
    assert [moments[k] for k in (2, 4, 6)] == [F(1, 4), F(1, 8), F(5, 64)]
    # Integrating a private link in W_q^2 gives 1/4 independently of W_p.
    mixed2 = moments[4]+191*moments[2]*moments[2]
    mixed4 = moments[6]+191*moments[2]*moments[4]
    assert mixed2 == F(193, 16) and mixed4 == F(387, 64)
    return dict(side=4, plaquettes=192, selected_face_edges=sorted(selected),
                distinct_face_boundaries=192, odd_link_pair_checks=len(pair_xors),
                no_three_plaquette_parity_cancellation=True,
                other_faces_with_a_private_link=len(private_edges), minimum_private_links=min(private_edges),
                private_link_conditional_square_integral=F(1, 4),
                whole_plaquette_independence_assumed=False,
                W2=moments[2], W4=moments[4], W6=moments[6],
                V_W2=F(0), V_W4=F(0), V2=F(48), V3=F(0), V_LE_V=F(144),
                V2_W2=mixed2, V2_W4=mixed4)


def dual_witness_bound(time, coefficient, source_pairing_lower, response_norm_squared_upper):
    time, coefficient = positive(time), nonnegative(coefficient)
    pairing, norm_squared = exact(source_pairing_lower), nonnegative(response_norm_squared_upper)
    value = 2*time*coefficient*pairing-coefficient*coefficient*norm_squared
    return dict(time=time, trial_coefficient=coefficient, signed_variational_value=value,
                certified_resolvent_norm_squared_lower=max(F(0), value),
                certified_heat_filtered_norm_squared_lower=max(F(0), value))


@lru_cache(maxsize=None)
def actual_cube_witness():
    moments = native_source_moments()
    spectral = stage75.actual_cube_controls()
    theta, alpha, beta = F(1, 48), F(1, 144), F(9, 50)
    E_lower = spectral['ground_threshold']['energy']
    E_upper = spectral['retained_ground_energy_upper']
    excited = spectral['first_excited_threshold']['energy']
    normalization_squared = 1+alpha*alpha*moments['V2']
    rayleigh = (alpha*alpha*moments['V_LE_V']-2*theta*alpha*moments['V2'])/normalization_squared
    assert normalization_squared == F(433, 432) and rayleigh == -F(3, 433)
    vacuum = vacuum_distance_budget(rayleigh, E_lower, excited)
    # <h,1+alpha V>=(1-alpha E0/theta)<h,1> >0 by the actual eigen-equation.
    assert E_upper <= 0 and 1-alpha*E_upper/theta > 0
    assert vacuum['excited_weight_upper'] <= F(9999, 25000000)
    assert normalization_squared <= F(501, 500)**2
    assert vacuum['distance_upper_after_certified_phase_alignment'] <= F(1, 50)
    assert 3 <= 3-E_upper <= 3-E_lower <= F(301, 100)
    F_squared_upper = (3-E_lower)**2*moments['W4']+theta*theta*moments['V2_W4']
    response_squared_upper = (4-E_lower)**2*moments['W2']+theta*theta*moments['V2_W2']
    assert F_squared_upper <= F(107, 100)**2
    assert response_squared_upper <= F(403, 100)
    hand_pairing = (F(3, 4)-theta*alpha*moments['V2_W2'])*F(500, 501)-F(107, 5000)
    assert hand_pairing >= F(29, 40)
    hand = dual_witness_bound(1, beta, F(29, 40), F(403, 100))
    assert hand['signed_variational_value'] == F(32607, 250000) > F(13, 100)
    normalization_upper = sqrt_interval(normalization_squared)['upper']
    approximation_pairing = ((3-E_upper)*moments['W2']-theta*alpha*moments['V2_W2'])/normalization_upper
    directed_pairing = approximation_pairing-vacuum['distance_upper_after_certified_phase_alignment']*sqrt_interval(F_squared_upper)['upper']
    directed = dual_witness_bound(1, beta, directed_pairing, response_squared_upper)
    assert directed['signed_variational_value'] > hand['signed_variational_value']
    return dict(theta=theta, side=4, source='(W_p-mu(W_p)) h', trial='z=(9/50) W_p in native Haar Hilbert space',
                ground_trial='g=(1+V/144)/sqrt(433/432)', normalization_squared=normalization_squared,
                vacuum_approximation=vacuum, phase_verified_by_actual_eigen_equation=True,
                ground_energy_interval=dict(lower=E_lower, upper=E_upper),
                Haar_moments=moments, hand_source_pairing_lower=F(29, 40),
                hand_response_norm_squared_upper=F(403, 100), rational_hand_witness=hand,
                directed_source_pairing_lower=directed_pairing, sharper_directed_witness=directed,
                resolvent_time=F(1), time_units='inverse bare electric Hamiltonian units',
                hidden_space_untruncated=True, trial_identified_with_actual_vacuum=False,
                uniform_volume_or_continuum_survival_asserted=False)


def source_residual_budget(time, observable_norm_upper, vacuum_distance_upper,
                           ground_interval_width, source_norm_upper,
                           response_norm_lower, response_norm_upper, residual_norm_upper,
                           observable_approximation_error=0):
    """K=A-L=H+d. y solves (I+tK)y≈Og; response is Og-y."""
    time = positive(time)
    M, eps, width, source_norm, lower, upper, residual, observable_error = map(nonnegative,
        (observable_norm_upper, vacuum_distance_upper, ground_interval_width, source_norm_upper,
         response_norm_lower, response_norm_upper, residual_norm_upper, observable_approximation_error))
    if lower > upper:
        raise ValueError('The computed response norm interval must be ordered')
    shift_price = time*width/(1+time*width)
    error = M*eps+observable_error+shift_price*source_norm+residual
    norm_lower, norm_upper = max(F(0), lower-error), upper+error
    return dict(time=time, vacuum_source_error=M*eps, observable_approximation_error=observable_error,
                commuting_ground_shift_operator_error=shift_price,
                ground_shift_source_error=shift_price*source_norm,
                solve_residual_error=residual, total_response_error=error,
                resolvent_response_norm=dict(lower=norm_lower, upper=norm_upper),
                resolvent_norm_squared=dict(lower=norm_lower*norm_lower, upper=norm_upper*norm_upper),
                lower_A_bound_and_phase_alignment_required=True,
                inverse_uses_full_neutral_operator_not_a_hidden_truncation=True)


@lru_cache(maxsize=None)
def residual_matrix_controls():
    rows = []
    O = [[F(0), F(1)], [F(1), F(0)]]
    H, h = a.diag((0, 2)), [F(1), F(0)]
    for parameter, width, solve_error in ((F(1, 10), F(1, 10), F(1, 100)),
                                         (F(1, 100), F(1, 100), F(0)),
                                         (F(0), F(0), F(0))):
        g = [(1-parameter*parameter)/(1+parameter*parameter), 2*parameter/(1+parameter*parameter)]
        eps = norm_interval([x-y for x, y in zip(g, h)])['upper']
        w = a.apply(O, g)
        K = a.add(H, a.identity(2), width)
        operator = a.add(a.identity(2), K)
        exact_solve = a.apply(a.inverse(operator), w)
        y = [exact_solve[0]+solve_error, exact_solve[1]-solve_error]
        residual = [x-y for x, y in zip(w, a.apply(operator, y))]
        response = [x-y for x, y in zip(w, y)]
        response_norm = norm_interval(response)
        budget = source_residual_budget(1, 1, eps, width, 1, response_norm['lower'],
                                        response_norm['upper'], norm_interval(residual)['upper'])
        v = a.apply(O, h)
        exact_response = [x-y for x, y in zip(v, a.apply(a.inverse(a.add(a.identity(2), H)), v))]
        exact_squared = sum(x*x for x in exact_response)
        assert exact_squared == F(4, 9)
        assert budget['resolvent_norm_squared']['lower'] <= exact_squared <= budget['resolvent_norm_squared']['upper']
        rows.append(dict(parameter=parameter, ground_interval_width=width,
                         approximate_vacuum=g, exact_source_response_norm_squared=exact_squared,
                         residual=residual, budget=budget, matrix_fixture_not_Yang_Mills=True))
    return rows


def weyl_increment_budget(increment_lower, source_mass_upper=1):
    increment, mass = nonnegative(increment_lower), positive(source_mass_upper)
    if mass > 1 or increment > 2*mass:
        raise ValueError('A centered unit-modulus source has mass <=1 and increment <=2 mass')
    return dict(increment_lower=increment, centered_source_mass_upper=mass,
                reflected_norm_lower=increment*increment/(4*mass),
                weaker_universal_reflected_lower=increment*increment/4,
                energy_measure_below_epsilon_bound='epsilon times centered source mass',
                actual_positive_increment_premise_still_required=True)


def raw_phase_survival_budget(time, raw_energy_upper, raw_increment_variance_lower):
    """q must bound pointwise physical Gamma, not merely the mean form energy."""
    time, q, variance = positive(time), positive(raw_energy_upper), positive(raw_increment_variance_lower)
    if variance > 2*q*time:
        raise ValueError('The proposed increment floor exceeds the actual energy upper bound')
    phase_squared = variance/(2*q*q*time*time)
    increment = variance*variance/(4*q*q*time*time)
    reflected = variance**4/(64*q**4*time**4)
    assert reflected == increment*increment/4
    return dict(time=time, raw_energy_upper=q, raw_increment_variance_lower=variance,
                raw_fourth_moment_upper=12*q*q*time*time, chosen_phase_parameter_squared=phase_squared,
                Weyl_increment_lower=increment, Weyl_reflected_norm_lower=reflected,
                Weyl_energy_upper=phase_squared*q,
                pointwise_carre_du_champ_bound_required_not_just_mean_energy=True,
                requires_actual_raw_increment_lower_bound=True,
                Gaussian_joint_independence_of_martingales_assumed=False)


def ou_phase_control(rate, raw_energy=1, time=1, phase_parameter_squared=1):
    rate, q, time, phase_squared = map(positive, (rate, raw_energy, time, phase_parameter_squared))
    phase_energy = phase_squared*q
    poisson_mean = phase_energy/rate
    def correlation(at_time):
        exp = stage76.exp_negative_interval(rate*at_time)
        exponent_lower = poisson_mean*(1-exp['upper'])
        exponent_upper = poisson_mean*(1-exp['lower'])
        return dict(lower=stage76.exp_negative_interval(exponent_upper)['lower'],
                    upper=stage76.exp_negative_interval(exponent_lower)['upper'])
    C1, C2 = correlation(time), correlation(2*time)
    increment = dict(lower=2*(1-C1['upper']), upper=2*(1-C1['lower']))
    zero_atom = stage76.exp_negative_interval(poisson_mean)
    mass = dict(lower=1-zero_atom['upper'], upper=min(F(1), poisson_mean, 1-zero_atom['lower']))
    lower = max(F(0), 1-2*C1['upper']+C2['lower'], increment['lower']**2/4)
    upper = min(mass['upper'], increment['upper']/2, 1-2*C1['lower']+C2['upper'])
    assert 0 <= lower <= upper <= 1
    raw = stage77.finite_filter_control((rate,), (q/rate,), time, q)
    return dict(rate=rate, raw_energy=q, time=time, phase_parameter_squared=phase_squared,
                Poisson_mean=poisson_mean, phase_spectral_lattice_spacing=rate,
                phase_energy=phase_energy, phase_second_moment=phase_energy**2+rate*phase_energy,
                phase_centered_mass=mass, phase_correlation_at_T=C1, phase_correlation_at_2T=C2,
                Weyl_increment=increment, Weyl_reflected_norm=dict(lower=lower, upper=upper),
                raw_linear_reflected_norm=raw['filtered_norm_squared'],
                raw_increment_variance=raw['ordinary_increment_variance'],
                OU_fixture_not_Yang_Mills_vacuum=True)


def clipped_source_budget(theta, spacing, projected_L2_squared, L1_squared_upper,
                           time=1, clock_scale=1):
    """Actual regularized-angle witness with eta=1/16 and audited B=2^18."""
    theta, spacing, time, scale = map(positive, (theta, spacing, time, clock_scale))
    projected, l1 = nonnegative(projected_L2_squared), nonnegative(L1_squared_upper)
    if theta < 25:
        raise ValueError('The regularized actual-vacuum estimate requires theta >=25')
    if l1 < 4*spacing**3*projected:
        raise ValueError('The proposed L1 bound is inconsistent with the rectangular cell projection')
    eta, gamma, B = F(1, 16), F(9, 10), 2**18
    q = scale*projected
    m0 = F(45, 16)*l1/spacing**2
    smooth_m2 = B*scale*scale*projected/spacing
    denominator = m0+2*time*q+time*time*smooth_m2
    coefficient = time*gamma*q/denominator if denominator else F(0)
    lower = time*time*gamma*gamma*q*q/denominator if denominator else F(0)
    assert 2*time*coefficient*gamma*q-coefficient*coefficient*denominator == lower
    old_m2 = 10*scale*scale*l1/spacing**4
    old_denominator = m0+time*q+time*time*old_m2
    old_lower = time*time*q*q/old_denominator if old_denominator else F(0)
    root = sqrt_interval(theta)
    return dict(theta=theta, spacing=spacing, sqrt_cE_cB=scale, time=time,
                regularizer_eta=eta, regularizer_epsilon_fourth_power=eta**4/theta,
                regularizer_epsilon_squared=dict(lower=eta*eta/root['upper'], upper=eta*eta/root['lower']),
                exact_raw_physical_energy=q, mixed_raw_regularized_energy_lower=gamma*q,
                regularized_energy_upper=q, source_norm_squared_upper=m0,
                regularized_second_moment_upper=smooth_m2, regularized_second_moment_constant=B,
                witness_denominator_upper=denominator, optimized_witness_coefficient=coefficient,
                raw_source_resolvent_norm_squared_lower=lower,
                raw_source_heat_filtered_norm_squared_lower=lower,
                previous_a4_lower_bound=old_lower,
                current_lower_bound_is_larger=lower > old_lower,
                normalized_a2_lower_coefficient=lower/spacing**2,
                fixed_data_a2_asymptotic_coefficient=(time*time*gamma*gamma*q*q/(F(45, 16)*l1) if l1 else F(0)),
                continuum_positive_lower_bound_proved=False)


@lru_cache(maxsize=None)
def clipped_regularizer_constants():
    """Replay Bell arithmetic; derivative estimates are proved in the manuscript.

    A_k bounds derivatives of the regularized angle, not powers of its Hessian.
    The finite partition enumeration checks the constants used by the analytic
    commutator estimate; it is not a replacement for that estimate.
    """
    Cj = (0, 1, 1, 4, 24, 204)
    bounds, partitions = {}, {}
    for k in range(2, 6):
        terms = []
        for multiplicities in product(*(range(k//ell+1) for ell in range(1, k+1))):
            if sum(ell*m for ell, m in enumerate(multiplicities, 1)) != k:
                continue
            j = sum(multiplicities)
            divisor, factor = 1, Cj[j]
            for ell, m in enumerate(multiplicities, 1):
                divisor *= factorial(m)*factorial(ell)**m
                if ell >= 2:
                    factor *= (2**ell*factorial(ell))**m
            assert factorial(k) % divisor == 0
            coefficient = factorial(k)//divisor
            terms.append(dict(multiplicities=list(multiplicities), outer_derivative_order=j,
                              Bell_coefficient=coefficient, term=coefficient*factor))
        bounds[k], partitions[k] = sum(row['term'] for row in terms), terms
    assert bounds == {2: 9, 3: 76, 4: 984, 5: 17484}
    K = 8*bounds[2]**2+4*bounds[3]
    Z = 12*bounds[5]+24*bounds[3]**2+8*bounds[2]*bounds[4]
    total = 9*K*256+Z*65536+8+24*16
    assert (K, Z, total) == (952, 419280, 27480127880)
    assert total < 2**36
    return dict(outer_derivative_bounds=list(Cj[1:]),
                analytic_composite_derivative_bounds={str(k): v for k, v in bounds.items()},
                Bell_partition_terms={str(k): v for k, v in partitions.items()},
                commutator_K=K, commutator_Z_upper=Z, squared_second_moment_constant_upper=total,
                second_moment_constant=2**18, squared_constant_majorant=2**36,
                analytic_derivative_and_commutator_estimates_required=True,
                finite_enumeration_proves_universal_derivative_estimates=False)


@lru_cache(maxsize=None)
def all_controls():
    IR = [ou_phase_control(F(1, n)) for n in (2, 4, 32, 256)]
    UV = [ou_phase_control(n) for n in (2, 4, 32, 256)]
    for n, row in zip((2, 4, 32, 256), IR):
        assert row['raw_linear_reflected_norm']['upper'] <= F(1, n)
        assert row['Weyl_reflected_norm']['lower'] > F(1, 4)
    for n, row in zip((2, 4, 32, 256), UV):
        assert row['Weyl_reflected_norm']['upper'] <= F(1, n)
    # For E=1/n,n>=2, exponent n(1-e^-1/n)>=3/4 and e^(3/4)>=65/32.
    uniform_IR_hand_lower = F(33, 65)**2
    assert uniform_IR_hand_lower > F(1, 4)
    return dict(actual_local_four_cube_witness=actual_cube_witness(),
                all_coupling_residual_matrix_controls=residual_matrix_controls(),
                phase_increment_bridge=[weyl_increment_budget(value) for value in (F(1, 10), F(1, 2), 1)],
                conditional_actual_raw_to_phase_bridge=[raw_phase_survival_budget(1, 1, value) for value in (F(1, 2), 1, F(3, 2))],
                OU_infrared_phase_rescue=IR, OU_ultraviolet_phase_failure=UV,
                OU_all_n_ge_2_IR_hand_reflected_lower=uniform_IR_hand_lower,
                clipped_regularizer_constant_replay=clipped_regularizer_constants(),
                actual_regularized_source_budgets=[clipped_source_budget(256, F(1, n), 1, 1)
                    for n in (4, 16, 32, 256, 65536, 2**24)])


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 288 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM79 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM79 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM78 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM79_LOCAL_WITNESS_SURVIVAL_CERTIFICATES',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_LOCAL_FINITE_WITNESS_AND_PHASE_BRIDGE__UNIFORM_CONTINUUM_SURVIVAL_OPEN',
        proof_scope='Written actual-model proofs with native Haar geometry and labelled exact residual/OU controls',
        open='Uniform positive physical local survival, regular Yang--Mills fields, Euclidean covariance and mass gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM79_LOCAL_WITNESS_SURVIVAL_CERTIFICATES':
        raise ValueError('YM79 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM79 fresh certificate/pin mismatch')
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
    print('YM79 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], actual_local_witness_lower='32607/250000',
                          uniform_continuum_survival=False)))


if __name__ == '__main__':
    main()
