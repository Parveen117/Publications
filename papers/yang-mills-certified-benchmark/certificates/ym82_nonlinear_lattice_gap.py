"""YM82: nonlinear quantum plaquette gaps and explicitly defined parent models.

Python 3.12 only. The oscillation theorem concerns L_E-sum Phi_p(W_p)
with its actual vacuum. The e^(S/2) vacuum is exact only for the separately
defined parent Hamiltonian. Finite controls replay written estimates.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym81_local_third_moment_vacuum_response as previous

native, p, a, model = previous.native, previous.p, previous.a, previous.model
exact, positive, nonnegative, encode = previous.exact, previous.positive, previous.nonnegative, previous.encode
stage80 = previous.previous
RESULT = HERE/'YM82_RESULT.json'
PIN = HERE/'EXPECTED_YM82.sha256'
SOURCES = HERE/'YM82_SOURCE_MANIFEST.json'
BASE_TREE = '48e5eb6efd097809e649c57b5faa4db675392ae2'
PARENT_HASH = 'e1720eeb6f2c6aed083e64bcbc8a9ad3bf712e4be5377d06bf0b9b08b13218e7'
MODEL = dict(actual_bounded_nonlinear_plaquette_gap_theorem=True,
             actual_strong_internal_block_quantum_gap_theorem=True,
             strong_block_gap_uniform_in_fixed_internal_frequency=True,
             strong_block_gap_uniform_in_internal_amplitude=False,
             actual_quantum_ground_measure_used=True,
             unknown_quantum_vacuum_replaced_by_classical_Gibbs=False,
             weighted_half_oscillation_controls_volume_uniform_gap=True,
             rejected_sufficient_gate_proves_gaplessness=False,
             upstream_linear_Yang_Mills_Hamiltonian_redefined=False,
             separately_defined_parent_has_exact_exp_S_over_2_vacuum=True,
             parent_gap_is_uniform_in_oscillation_frequency=True,
             global_Bakry_curvature_positivity_is_necessary=False,
             parent_vacuum_identified_with_original_plaquette_vacuum=False,
             quartic_doublewell_gap_has_exponentially_small_upper_bound=True,
             doublewell_exact_vacuum_applies_to_overlapping_all_plaquettes=False,
             all_nonlinear_models_have_one_uniform_positive_gap=False,
             perturbative_remainders_certified_uniform_in_volume=False,
             spatial_continuum_Yang_Mills_constructed=False,
             positive_continuum_mass_gap=False,
             finite_controls_are_universal_proof=False)


def log_interval(value, terms=24):
    value = positive(value)
    native.width(terms, 1)
    reverse = value < 1
    if reverse:
        value = 1/value
    q = (value-1)/(value+1)
    lower = 2*sum((q**(2*k+1)/F(2*k+1) for k in range(terms)), F(0))
    error = 2*q**(2*terms+1)/(F(2*terms+1)*(1-q*q))
    return dict(lower=-lower-error, upper=-lower) if reverse else dict(lower=lower, upper=lower+error)


def quadratic_halfosc(theta, lam):
    theta, lam = exact(theta), exact(lam)
    t, ell = abs(theta), abs(lam)
    candidates = [(-F(1), -theta+lam), (F(1), theta+lam)]
    critical = -theta/(2*lam) if lam and t <= 2*ell else None
    if critical is not None:
        candidates.append((critical, theta*critical+lam*critical*critical))
    minimum, maximum = min(value for _, value in candidates), max(value for _, value in candidates)
    halfosc = (maximum-minimum)/2
    formula = t if t >= 2*ell else (t+2*ell)**2/(8*ell)
    assert halfosc == formula and halfosc >= t
    return dict(theta=theta, lambda_coefficient=lam, absolute_theta=t, absolute_lambda=ell,
                minimum=minimum, maximum=maximum, midpoint=(minimum+maximum)/2,
                halfosc=halfosc, critical_point=critical,
                halfosc_dominates_linear_amplitude=True,
                Lipschitz_upper=t+2*ell, mixed_native_derivative_coefficient=t+2*ell,
                constant_shift_does_not_change_gap=True)


def nonlinear_gap_budget(load, J=5, R=None):
    """S=max_e sum_(p containing e) osc(Phi_p)/2, with four-link faces."""
    load = nonnegative(load)
    native.width(J, 1)
    alpha = F(5, 2*J)+200*load*J
    if R is None:
        R = (1+1/alpha)/2 if alpha < 1 else None
    else:
        R = positive(R)
        if R <= 1:
            raise ValueError('The time weight ratio must exceed one')
    accepted = R is not None and R*alpha < 1
    logarithm = log_interval(R) if accepted else None
    gamma = {key: value/(10*J) for key, value in logarithm.items()} if logarithm else None
    return dict(weighted_halfosc_load=load, J=J, alpha_upper=alpha,
                time_weight_ratio=R, accepted=accepted,
                weighted_margin=1-R*alpha if R is not None else None,
                gamma_interval=gamma, sufficient_load_threshold=F(1, 2000),
                temporal_endpoint_load_coefficient=8, spatial_outgoing_load_coefficient=12,
                fine_total_load_coefficient=20, coarse_total_load_coefficient=200,
                full_quantum_Hilbert_space_gap=True if accepted else None,
                gauge_invariant_restriction_preserves_bound=True,
                actual_ground_measure_not_a_substituted_Gibbs_measure=True,
                finite_volume_uniform=True if accepted else None,
                rejection_proves_gaplessness=False)


def quadratic_gap_budget(theta, lam, J=5, R=None):
    interaction = quadratic_halfosc(theta, lam)
    row = nonlinear_gap_budget(4*interaction['halfosc'], J, R)
    row.update(interaction=interaction, plaquettes_per_link=4,
               sufficient_uniform_halfosc_threshold=F(1, 8000),
               potential='Phi(W)=theta W+lambda W^2')
    return row


def strong_block_gap_budget(B0, max_block_links, degree, arity, weak_coupling,
                            J=None, epsilon_lower=None, R=None):
    """Actual block Hamiltonians, each with osc(U_B)<=B0, plus weak bridges.

    The internal and bridge potentials are smooth and real. Each bridge has
    absolute value <=1, meets <=arity blocks and has finite carré du champ;
    each block meets <=degree bridges. No internal spectral gap is assumed.
    """
    B0, weak_coupling = nonnegative(B0), exact(weak_coupling)
    for size in (max_block_links, degree, arity):
        native.width(size, 1)
    heat = stage80.exp_negative_interval(10*B0)
    # The inherited grid may round a very small exponential lower bound to
    # zero. The elementary e<3 bound retains a positive rational certificate
    # for every finite B0, with no claim of amplitude-uniform constants.
    heat_lower = max(heat['lower'], F(1, 3**((10*B0).__ceil__())))
    free_ratio = F(19, 21)**max_block_links
    minorization = dict(lower=heat_lower*free_ratio, upper=heat['upper']*free_ratio)
    assert 0 < minorization['lower'] <= minorization['upper'] <= 1
    epsilon = minorization['lower'] if epsilon_lower is None else positive(epsilon_lower)
    if epsilon > minorization['lower']:
        raise ValueError('The supplied minorization floor exceeds its directed certificate')
    recommended_J = (4/epsilon).__ceil__()
    if J is None:
        J = recommended_J
    native.width(J, 1)
    temporal = 2/(epsilon*J)
    weak = 40*degree*(arity+1)*abs(weak_coupling)*J
    alpha = temporal+weak
    if R is None:
        R = (1+1/alpha)/2 if alpha < 1 else None
    else:
        R = positive(R)
        if R <= 1:
            raise ValueError('The time weight ratio must exceed one')
    accepted = R is not None and R*alpha < 1
    logarithm = log_interval(R) if accepted else None
    gamma = {key: value/(10*J) for key, value in logarithm.items()} if logarithm else None
    threshold = F(1, 80*degree*(arity+1)*J)
    return dict(internal_oscillation_bound=B0, max_block_links=max_block_links,
                bridge_degree=degree, bridge_arity=arity, weak_coupling=weak_coupling,
                free_kernel_ratio=free_ratio, minorization_interval=minorization,
                epsilon_lower=epsilon, J=J, recommended_J=recommended_J,
                temporal_alpha_upper=temporal, weak_alpha_upper=weak, alpha_upper=alpha,
                sufficient_weak_coupling_threshold=threshold,
                simple_sufficient_gate=J >= recommended_J and abs(weak_coupling) < threshold,
                time_weight_ratio=R, weighted_margin=1-R*alpha if R is not None else None,
                accepted=accepted, gamma_interval=gamma,
                full_quantum_Hilbert_space_gap=True if accepted else None,
                finite_volume_uniform=True if accepted else None,
                gauge_sector_restriction_requires_gauge_invariant_potentials=True,
                actual_block_ground_transform_used=True,
                internal_quantum_vacuum_assumed_Gibbs=False,
                local_spectral_gap_assumed=False,
                smooth_potentials_and_finite_bridge_Gamma_required=True,
                uniform_in_internal_amplitude=False,
                rejection_proves_gaplessness=False)


def strong_cosine_gap_budget(n, amplitude=F(1, 8), weak_coupling=F(1, 1228800)):
    """Selected four-link U_B=-amplitude*cos(nW), with weak other faces."""
    native.width(n, 1)
    amplitude = nonnegative(amplitude)
    extra = dict(J=128, epsilon_lower=F(1, 32)) if amplitude == F(1, 8) else {}
    row = strong_block_gap_budget(2*amplitude, 4, 12, 4, weak_coupling, **extra)
    row.update(n=n, amplitude=amplitude,
               potential='-A sum_selected cos(nW) - u sum_nonselected W',
               quantum_Hamiltonian='L_E plus the displayed multiplication potential',
               interacting_ground_state_assumed_exp_S_over_2=False,
               gauge_invariant_restriction_preserves_bound=True,
               gap_bound_independent_of_fixed_frequency=True,
               bridge_Gamma_finite=True, native_plaquette_Gamma_upper=F(1))
    return row


@lru_cache(maxsize=None)
def nonlinear_partition_geometry(side):
    native.width(side, 4)
    if side % 2:
        raise ValueError('The selected periodic plaquette partition requires even side length')
    graph = model.lattice(side)
    index = {edge: i for i, edge in enumerate(graph['edges'])}
    def shift(vertex, axis):
        out = list(vertex)
        out[axis] = (out[axis]+1) % side
        return tuple(out)
    owners, supports, selected_faces = {}, [], set()
    for x, y, z in product(range(0, side, 2), range(0, side, 2), range(side)):
        vertex = (x, y, z)
        support = frozenset((index[vertex, 0], index[shift(vertex, 0), 1],
                             index[shift(vertex, 1), 0], index[vertex, 1]))
        block = len(supports)
        supports.append(support)
        selected_faces.add(support)
        for edge in support:
            assert edge not in owners
            owners[edge] = block
    selected_count = len(supports)
    for edge in range(len(graph['edges'])):
        if edge not in owners:
            owners[edge] = len(supports)
            supports.append(frozenset((edge,)))
    degrees, row_mixed_counts = Counter(), Counter()
    cross_faces = 0
    for face in graph['faces']:
        support = frozenset(edge for edge, _ in face)
        if support in selected_faces:
            continue
        block_incidence = Counter(owners[edge] for edge in support)
        assert len(block_incidence) == 4 and set(block_incidence.values()) == {1}
        cross_faces += 1
        for block in block_incidence:
            degrees[block] += 1
            row_mixed_counts[block] += 3
    assert selected_count == side**3//4
    assert len(supports)-selected_count == 2*side**3
    assert set(degrees[block] for block in range(selected_count)) == {12}
    assert set(degrees[block] for block in range(selected_count, len(supports))) == {4}
    return dict(side=side, selected_four_link_blocks=selected_count,
                remaining_singleton_blocks=len(supports)-selected_count,
                total_links=len(graph['edges']), nonselected_faces=cross_faces,
                every_nonselected_face_meets_four_distinct_blocks=True,
                selected_block_cross_degree=12, singleton_cross_degree=4,
                mixed_derivative_bound_per_face_pair=F(1, 4),
                selected_block_mixed_row_per_abs_b=F(max(row_mixed_counts[b] for b in range(selected_count)), 4),
                singleton_mixed_row_per_abs_b=F(max(row_mixed_counts[b] for b in range(selected_count, len(supports))), 4),
                selected_support_potential_Hessian_lower_coefficient=F(-3),
                finite_geometry_replay_not_full_gap_proof=True)


def block_gap_matrix(conditional_gaps, cross_bounds, gap_target):
    gaps = tuple(map(positive, conditional_gaps))
    target = nonnegative(gap_target)
    cross = [list(map(nonnegative, row)) for row in cross_bounds]
    if not gaps or len(cross) != len(gaps) or any(len(row) != len(gaps) for row in cross):
        raise ValueError('Conditional and cross-block dimensions must match')
    if cross != a.transpose(cross) or any(cross[i][i] for i in range(len(gaps))):
        raise ValueError('The cross bounds must be symmetric with zero diagonal')
    matrix = a.add(a.diag(gaps), cross, -1)
    shifted = a.add(matrix, a.identity(len(gaps)), -target)
    accepted = a.psd(shifted)
    return dict(comparison_matrix=matrix, shifted_comparison_matrix=shifted,
                gap_target=target, accepted=accepted, positive_gap_certified=(accepted and target > 0),
                conditional_Poincare_and_cross_Hessian_premises_required=True,
                matrix_fixture_not_an_actual_Yang_Mills_vacuum=True)


def parent_gap_budget(n, A=F(1, 8), b=F(1, 192)):
    native.width(n, 1)
    A, b = nonnegative(A), nonnegative(b)
    oscillation = 2*A+24*b
    exponential = stage80.exp_negative_interval(oscillation)
    conditional = {key: F(3, 4)*value for key, value in exponential.items()}
    gap = {key: value-9*b for key, value in conditional.items()}
    elementary = F(3, 4)*(1-oscillation)-9*b
    gap['lower'] = max(gap['lower'], elementary)
    # For n>=4, a holonomy with W=pi/n exists and pi^2<10.
    Hessian_lower = A*(n*n-10)-3*b if n >= 4 else None
    return dict(n=n, A=A, b=b, log_density='S=A sum_selected cos(n W)+b sum_nonselected W',
                parent_ground_state='h_S=Z^(-1/2) exp(S/2)',
                parent_Hamiltonian='K_S=L_E-(1/2)L_E S+(1/4)Gamma(S)',
                conditional_log_density_oscillation_upper=oscillation,
                conditional_gap_interval=conditional, cross_Hessian_row_upper=9*b,
                parent_gap_interval=gap, elementary_parent_gap_lower=elementary,
                interval_encloses_comparison_rate_not_actual_gap=True,
                positive_gap_certified=(gap['lower'] > 0),
                gap_lower_independent_of_n=True,
                Hessian_S_direction_lower=Hessian_lower,
                product_Ricci_in_direction=F(1, 2),
                global_positive_Bakry_bound_fails=(Hessian_lower is not None and Hessian_lower > F(1, 2)),
                identified_with_original_plaquette_Hamiltonian=False)


@lru_cache(maxsize=None)
def native_nonlinear_controls():
    """Exact native one-holonomy calculus; a four-link loop has four Casimirs.

    The inherited full ordered four-link identities identify this radial
    subalgebra. Polynomial identities retain the radius until the unit
    carrier is imposed; Haar integration is native, not a discrete fixture.
    """
    inherited = model.plaquette_controls()
    w = p.COORD[0]
    radius = p.add(*(p.mul(q, q) for q in p.COORD))
    one = {(0, 0, 0, 0): F(1)}
    w2 = p.mul(w, w)
    w3, w4 = p.mul(w2, w), p.mul(w2, w2)
    def electric(poly):
        return p.scale(model.electric(poly, 1), 4)
    def energy(poly):
        return p.scale(model.energy(poly, 1), 4)
    Y = p.add(w2, p.scale(radius, -F(1, 4)))
    assert electric(w) == p.scale(w, 3)
    assert electric(Y) == p.scale(Y, 8)
    assert electric(w3) == p.add(p.scale(w3, 15), p.scale(p.mul(w, radius), -6))
    assert energy(w) == p.add(radius, p.scale(w2, -1))
    assert native.phi(Y, 1) == 0 and native.phi(p.mul(Y, Y), 1) == F(1, 16)
    assert native.phi(p.mul(w, Y), 1) == 0
    parameter_rows = []
    for t in (F(0), F(1, 2), F(4), F(16)):
        S = p.scale(w2, t)
        parent_potential = p.add(p.scale(electric(S), F(1, 2)), p.scale(energy(S), -F(1, 4)))
        F_t = p.add(p.scale(w4, t*t), p.scale(w2, 4*t-t*t))
        difference = p.add(parent_potential, p.scale(F_t, -1), p.scale(one, t))
        unit_remainder = p.mul(p.add(radius, p.scale(one, -1)),
                               p.add(p.scale(one, -t), p.scale(w2, -t*t)))
        assert difference == unit_remainder
        parameter_rows.append(dict(t=t, ground_energy=-t,
                                   identity_holds_modulo_native_unit_radius=True))
    return dict(inherited_full_four_link_controls=inherited,
                native_holonomy_reduction_uses_four_link_Casimirs=True,
                linear_mode_rate=F(3), centered_quadratic_mode_rate=F(8),
                linear_Haar_variance=F(1, 4), centered_quadratic_Haar_variance=F(1, 16),
                cross_mode_Haar_inner_product=F(0),
                cubic_electric_identity='L_E W^3=15 W^3-6 W on the unit carrier',
                exact_doublewell_parent_identities=parameter_rows)


def nonlinear_ward_response(theta=1, lam=1):
    theta, lam = exact(theta), exact(lam)
    controls = native_nonlinear_controls()
    tau = theta*theta/12+lam*lam/128
    R = 5*theta*theta/24+15*lam*lam/256
    assert R == F(3, 2)*theta*(theta/6)+4*lam*(lam/64)-tau/2
    return dict(theta=theta, lambda_coefficient=lam,
                exact_Ward_identity='R+tau/2=(3/2)mu(W Phi_prime)-(1/2)mu((1-W^2)Phi_second)',
                exact_quadratic_Ward_identity='R=(3/2)theta mu(W)+4 lambda mu(W^2)-lambda-tau/2',
                mean_potential_Haar=lam/4,
                first_linear_response=theta/6, first_centered_quadratic_response=lam/64,
                ground_energy_per_face_through_order_two=-lam/4-theta*theta/12-lam*lam/128,
                electric_energy_per_link_quadratic=tau,
                log_vacuum_Hessian_row_quadratic=R,
                mixed_theta_lambda_second_order_coefficient=F(0),
                source_modes=controls,
                response_coefficients_are_fixed_volume_perturbation=True,
                perturbative_remainder_uniform_in_volume=False,
                score_derivatives_are_physical_gauge_invariant_states=False)


def quartic_doublewell_budget(t):
    t = nonnegative(t)
    heat = stage80.exp_negative_interval(t/2)
    cap = F(55, 2048)
    assert cap > F(2, 81)
    assert F(7, 8)**2 >= F(3, 4)
    assert F(11, 32)**2 <= 1-F(15, 16)**2
    assert F(16, 5)**2 > 10
    return dict(t=t, W4_coefficient=t*t, W2_coefficient=4*t-t*t,
                exact_ground_state='h=Z^(-1/2) exp(t W^2/2)', exact_ground_energy=-t,
                gauge_invariant_trial='clip(2W,-1,1)',
                cap_intervals=[[F(-15, 16), F(-7, 8)], [F(7, 8), F(15, 16)]],
                rational_native_Haar_cap_lower=cap, simplified_cap_lower=F(2, 81),
                numerator_bound='4 exp(t/4)', denominator_bound='(2/81) exp(3t/4)',
                exponential_interval=heat, gauge_invariant_gap_upper=162*heat['upper'],
                analytic_gap_upper='162 exp(-t/2)', inherited_pi_squared_upper=F(10),
                exact_vacuum_requires_one_or_edge_disjoint_active_plaquettes=True,
                exact_vacuum_claimed_for_overlapping_all_plaquettes=False,
                added_W_even_sector_would_exclude_trial=True,
                zero_gap_at_any_fixed_finite_t_asserted=False)


@lru_cache(maxsize=None)
def all_controls():
    example = quadratic_gap_budget(F(1, 16000), F(1, 16000))
    assert example['interaction']['halfosc'] == F(9, 128000)
    assert example['weighted_halfosc_load'] == F(9, 32000)
    assert example['alpha_upper'] == F(25, 32)
    assert example['time_weight_ratio'] == F(57, 50)
    assert example['weighted_margin'] == F(7, 64)
    assert example['gamma_interval']['lower'] > F(7, 2675)
    parents = [parent_gap_budget(n) for n in (4, 8, 32, 256)]
    assert all(row['elementary_parent_gap_lower'] == F(27, 64) for row in parents)
    assert parents[0]['Hessian_S_direction_lower'] == F(47, 64)
    strong = [strong_cosine_gap_budget(n) for n in (4, 8, 32, 256)]
    assert all(row['alpha_upper'] == F(3, 4) for row in strong)
    assert all(row['time_weight_ratio'] == F(7, 6) for row in strong)
    assert all(row['gamma_interval']['lower'] > F(1, 8320) for row in strong)
    return dict(actual_quadratic_nonlinear_gap_example=example,
                actual_strong_block_cosine_gap_controls=strong,
                strong_block_barrier_controls=[strong_block_gap_budget(B0, 4, 12, 4, 0)
                                               for B0 in (0, F(1, 4), 1, 8)],
                quadratic_halfosc_controls=[quadratic_halfosc(theta, lam) for theta, lam in
                    ((0, 0), (1, 0), (0, 1), (1, 1), (1, -1), (-1, 1), (-1, -1), (2, 1), (3, 1))],
                sufficient_gate_boundaries=[nonlinear_gap_budget(load) for load in (0, F(1, 4000), F(1, 2000), F(1, 1000))],
                actual_native_block_partitions=[nonlinear_partition_geometry(side) for side in (4, 6, 8)],
                separate_parent_gap_controls=parents,
                conditional_block_matrix_fixture=block_gap_matrix((2, 5), ((0, 1), (1, 0)), 1),
                native_nonlinear_derivative_controls=native_nonlinear_controls(),
                nonlinear_Ward_and_Haar_response=nonlinear_ward_response(),
                exact_doublewell_counterfamily=[quartic_doublewell_budget(t) for t in (0, 4, 16, 32, 64)])


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 309 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM82 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM82 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM81 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM82_NONLINEAR_LATTICE_GAP', verdict='PASS',
        runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_NONLINEAR_FIXED_LATTICE_GAP_AND_SEPARATE_PARENT_THEOREMS__CONTINUUM_OPEN',
        proof_scope='Written quantum oscillation, actual strong-block and separate parent gap theorems with native polynomial and geometry controls',
        open='Uniform continuum survival, continuum field identification, Euclidean covariance and Yang--Mills mass gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM82_NONLINEAR_LATTICE_GAP':
        raise ValueError('YM82 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM82 fresh certificate/pin mismatch')
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
    print('YM82 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], nonlinear_gap_lower='>7/2675',
                          actual_strong_block_gap_lower='>1/8320',
                          separate_parent_gap_lower='27/64', continuum_mass_gap=False)))


if __name__ == '__main__':
    main()
