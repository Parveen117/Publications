"""YM81: actual vacuum response, local third moments and rational windows.

Python 3.12 only. Native geometric controls and exact Ward budgets concern
the unchanged SU(2) vacuum. Circulant examples are labelled fixtures. The
zero-momentum bound needed for uniform local survival remains a premise.
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
import ym80_bounded_phase_energy_certificates as previous

stage79 = previous.previous
native, p, a = previous.native, previous.p, stage79.a
model = stage79.model
exact, positive, nonnegative, encode = previous.exact, previous.positive, previous.nonnegative, previous.encode
RESULT = HERE/'YM81_RESULT.json'
PIN = HERE/'EXPECTED_YM81.sha256'
SOURCES = HERE/'YM81_SOURCE_MANIFEST.json'
BASE_TREE = 'fa8e2cc81dbe343aaef0403012f36ef2ad3102cc'
PARENT_HASH = '0f84160b93477194c7f18cec524f37d084b1cc7f46bbf67b1f58cc945a37f17f'
MODEL = dict(actual_vacuum_Hessian_row_Ward_identity=True,
             native_Ricci_connection_correction_required=True,
             vacuum_score_vectors_are_physical_gauge_invariant_states=False,
             actual_smooth_source_third_moment_has_nearest_neighbor_support=True,
             actual_physical_spatial_gradient_remainder_uniformly_bounded=True,
             zero_mean_smearing_removes_zero_mode_intercept=False,
             actual_uniform_scaled_zero_mode_upper_bound_proved=False,
             zero_mode_bound_conditionally_implies_local_phase_survival=True,
             actual_four_cube_rational_window_nonzero=True,
             rational_window_has_uniform_energy_and_inverse_energy_budgets=True,
             actual_smooth_phase_window_has_order_a3_lower_certificate=True,
             rational_filtering_alone_proves_uniform_continuum_survival=False,
             uniform_positive_continuum_survival=False,
             regular_local_continuum_Yang_Mills_constructed=False,
             full_Euclidean_covariance=False, positive_continuum_mass_gap=False,
             finite_controls_are_universal_proof=False)


def epsilon(i, j, k):
    if len({i, j, k}) < 3:
        return F(0)
    return F(1 if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1)


@lru_cache(maxsize=None)
def lie_ricci_control():
    connection = [[[epsilon(i, j, k)/2 for k in range(3)] for j in range(3)] for i in range(3)]
    def curvature(i, j, k, ell):
        return sum((connection[j][k][v]*connection[i][v][ell]
                    -connection[i][k][v]*connection[j][v][ell]
                    -epsilon(i, j, v)*connection[v][k][ell] for v in range(3)), F(0))
    ricci = [[sum((curvature(i, j, k, i) for i in range(3)), F(0))
              for k in range(3)] for j in range(3)]
    assert ricci == a.scale(a.identity(3), F(1, 2))
    generators = []
    for axis in range(3):
        matrix = [[F(0) for _ in range(4)] for _ in range(4)]
        matrix[0][axis+1], matrix[axis+1][0] = -F(1, 2), F(1, 2)
        for j in range(3):
            for k in range(3):
                matrix[k+1][j+1] = epsilon(axis, j, k)/2
        generators.append(matrix)
    casimir = [[-sum((a.multiply(g, g)[i][j] for g in generators), F(0))
                for j in range(4)] for i in range(4)]
    assert casimir == a.scale(a.identity(4), F(3, 4))
    gradient = [F(2, 3), -F(3, 5), F(7, 4)]
    Hessian = [[F(1), F(1, 3), F(-2)], [F(1, 3), F(2), F(3, 7)], [F(-2), F(3, 7), F(-1)]]
    ordered = [[Hessian[i][j]+sum(connection[i][j][k]*gradient[k] for k in range(3))
                for j in range(3)] for i in range(3)]
    correction = sum(x*x for row in ordered for x in row)-sum(x*x for row in Hessian for x in row)
    assert correction == sum(x*x for x in gradient)/2
    return dict(connection=connection, Ricci=ricci,
                fundamental_representation_matrices=generators,
                native_coordinate_action_matrices=[a.scale(g, -1) for g in generators],
                fundamental_Casimir=casimir, score_gradient_fixture=gradient,
                symmetric_Hessian_fixture=Hessian, ordered_derivative_fixture=ordered,
                ordered_minus_Hessian_squared_norm=correction,
                same_link_connection_correction_coefficient=F(1, 2),
                score_vectors_not_assumed_gauge_invariant=True)


def ward_response_budget(theta, spacing=1, clock_scale=1):
    theta, spacing, scale = map(positive, (theta, spacing, clock_scale))
    if theta < 25:
        raise ValueError('The actual weak-coupling row bound requires theta >=25')
    root = previous.sqrt_interval(theta)
    plaquette_lower = 1-F(9, 4)/root['lower']
    tau_lower, tau_upper = root['lower']/32, F(9, 4)*root['upper']
    refined_lower = F(3, 2)*theta*plaquette_lower-tau_upper/2
    bare_lower, bare_upper = max(F(3, 5)*theta, refined_lower), F(3, 2)*theta
    return dict(theta=theta, spacing=spacing, sqrt_cE_cB=scale,
                mean_plaquette_lower=plaquette_lower, mean_plaquette_upper=F(1),
                electric_energy_per_link=dict(lower=tau_lower, upper=tau_upper),
                Hessian_row_identity='R=(3/2) theta mu(W)-(1/2) tau',
                Hessian_row_definition='sum_f mu(||Hess_(f,e) log(h)||_HS^2), fixed link e',
                bare_Hessian_row=dict(lower=bare_lower, upper=bare_upper),
                universal_bare_Hessian_row=dict(lower=F(3, 5)*theta, upper=bare_upper),
                physical_Hessian_row=dict(lower=F(3, 5)*scale**2/spacing**2,
                                           upper=F(3, 2)*scale**2/spacing**2),
                actual_gauge_invariant_vacuum=True, score_vectors_are_gauge_invariant_states=False,
                positive_local_gap_inferred=False)


@lru_cache(maxsize=None)
def native_third_moment_geometry(side):
    native.width(side, 4)
    if side % 2:
        raise ValueError('The selected periodic plaquettes require even side length')
    graph = model.lattice(side)
    edge_index = {edge: index for index, edge in enumerate(graph['edges'])}
    shape = (side//2, side//2, side)
    def shifted(vertex, axis):
        out = list(vertex)
        out[axis] = (out[axis]+1) % side
        return tuple(out)
    selected, owner = {}, {}
    for cell in product(*(range(length) for length in shape)):
        vertex = (2*cell[0], 2*cell[1], cell[2])
        support = frozenset((edge_index[vertex, 0], edge_index[shifted(vertex, 0), 1],
                             edge_index[shifted(vertex, 1), 0], edge_index[vertex, 1]))
        selected[cell] = support
        for edge in support:
            assert edge not in owner
            owner[edge] = cell
    pair_counts, axis_counts = Counter(), Counter()
    max_blocks, bridge_incidence = 0, []
    for face in graph['faces']:
        incidence = Counter(owner[edge] for edge, _ in face if edge in owner)
        max_blocks = max(max_blocks, len(incidence))
        assert len(incidence) <= 2
        if len(incidence) == 2:
            first, second = sorted(incidence)
            assert sorted(incidence.values()) == [1, 1]
            axes = [axis for axis in range(3) if first[axis] != second[axis]]
            assert len(axes) == 1
            axis_counts['xyz'[axes[0]]] += 1
            pair_counts[first, second] += 1
            bridge_incidence.append([1, 1])
    expected = Counter()
    for cell in selected:
        for axis, multiplicity in enumerate((1, 1, 4)):
            target = list(cell)
            target[axis] = (target[axis]+1) % shape[axis]
            expected[tuple(sorted((cell, tuple(target))))] += multiplicity
    assert pair_counts == expected
    count = len(selected)
    assert dict(axis_counts) == dict(x=count, y=count, z=4*count)
    origin = (0, 0, 0)
    neighbors = [dict(cell=list(second if first == origin else first), bridge_faces=number)
                 for (first, second), number in sorted(pair_counts.items()) if origin in (first, second)]
    return dict(side=side, coarse_shape=list(shape), selected_plaquettes=count,
                occupied_links=len(owner), edge_supports_pairwise_disjoint=True,
                maximum_selected_blocks_per_potential_face=max_blocks,
                bridge_faces_by_axis=dict(axis_counts), origin_neighbor_bridges=neighbors,
                every_bridge_shares_one_link_with_each_block=True,
                parallel_horizontal_bonds_counted=(side == 4),
                native_mixed_W_derivative_upper=F(1, 16),
                coupling_over_theta_upper=dict(x=F(1, 8), y=F(1, 8), z=F(1, 2)),
                off_diagonal_identity='T3_bc=-2 theta mu(X_b X_c V)',
                finite_geometry_replay_not_universal_operator_proof=True)


def third_moment_symbol(diagonal, couplings, cosines):
    diagonal = exact(diagonal)
    couplings, cosines = tuple(map(exact, couplings)), tuple(map(exact, cosines))
    if len(couplings) != 3 or len(cosines) != 3 or any(abs(value) > 1 for value in cosines):
        raise ValueError('Three couplings and three cosines in [-1,1] are required')
    zero = diagonal+2*sum(couplings)
    value = diagonal+2*sum(K*cosine for K, cosine in zip(couplings, cosines))
    assert value == zero-2*sum(K*(1-cosine) for K, cosine in zip(couplings, cosines))
    return dict(symbol=value, zero_momentum_coefficient=zero,
                scalar_intercept_present_at_every_momentum=True,
                fixture_not_actual_vacuum_moments=True)


def coefficient_gradient_identity(shape, coefficients, diagonal, couplings):
    shape = tuple(shape)
    if len(shape) != 3:
        raise ValueError('A three-dimensional periodic shape is required')
    for length in shape:
        native.width(length, 2)
    cells = tuple(product(*(range(length) for length in shape)))
    coefficients = tuple(map(exact, coefficients))
    couplings, diagonal = tuple(map(exact, couplings)), exact(diagonal)
    if len(coefficients) != len(cells) or len(couplings) != 3:
        raise ValueError('Coefficient/coupling dimensions do not match the periodic grid')
    values = dict(zip(cells, coefficients))
    norm = sum(value*value for value in coefficients)
    differences, correlations = [], []
    for axis in range(3):
        diff, corr = F(0), F(0)
        for cell in cells:
            target = list(cell)
            target[axis] = (target[axis]+1) % shape[axis]
            diff += (values[tuple(target)]-values[cell])**2
            corr += values[tuple(target)]*values[cell]
        differences.append(diff)
        correlations.append(corr)
    zero = diagonal+2*sum(couplings)
    matrix_form = diagonal*norm+2*sum(K*corr for K, corr in zip(couplings, correlations))
    gradient_form = zero*norm-sum(K*diff for K, diff in zip(couplings, differences))
    assert matrix_form == gradient_form
    return dict(shape=list(shape), coefficient_sum=sum(coefficients), coefficient_squared_norm=norm,
                gradient_squared_by_axis=dict(zip('xyz', differences)),
                matrix_quadratic_form=matrix_form, gradient_quadratic_form=gradient_form,
                zero_momentum_coefficient=zero, scalar_intercept_contribution=zero*norm,
                diagonal_dominance_certifies_PSD=diagonal >= 2*sum(abs(K) for K in couplings),
                circulant_fixture_not_actual_Yang_Mills_vacuum=True)


def matrix_moment_compression(energies, columns):
    """Exact isometric-compression fixture; no interacting moments are supplied."""
    energies = tuple(map(positive, energies))
    columns = [list(map(exact, row)) for row in columns]
    if not energies or len(columns) != len(energies) or not columns[0]:
        raise ValueError('A nonempty energy list and rectangular source matrix are required')
    width = len(columns[0])
    if any(len(row) != width for row in columns):
        raise ValueError('The source matrix must be rectangular')
    moments = [a.multiply(a.transpose(columns), a.multiply(a.diag([E**order for E in energies]), columns))
               for order in (1, 2, 3)]
    u = moments[0][0][0]
    if u <= 0 or moments[0] != a.scale(a.identity(width), u):
        raise ValueError('The first source moment must equal u times the identity with u>0')
    deficit = a.add(a.scale(moments[2], u), a.multiply(moments[1], moments[1]), -1)
    assert a.psd(deficit)
    commutator = a.add(a.multiply(moments[1], moments[2]), a.multiply(moments[2], moments[1]), -1)
    return dict(energies=list(energies), source_columns=columns,
                M1=moments[0], M2=moments[1], M3=moments[2], u=u,
                compression_deficit=deficit, compression_deficit_PSD=True,
                commutator=commutator,
                matrix_fixture_not_actual_Yang_Mills_vacuum=True)


def physical_third_moment_budget(spacing, projected_L2_squared, gradient_L2_squared,
                                 scaled_zero_mode_upper, clock_scale=1):
    spacing, scale = positive(spacing), positive(clock_scale)
    projected, gradient = nonnegative(projected_L2_squared), nonnegative(gradient_L2_squared)
    C = positive(scaled_zero_mode_upper)
    q = scale*projected
    remainder = scale**3*gradient/2
    return dict(spacing=spacing, exact_raw_physical_energy=q,
                assumed_scaled_zero_mode_upper=C, spatial_gradient_remainder_upper=remainder,
                physical_third_moment_upper=C*q+remainder,
                spatial_gradient_constant=scale**3/2,
                actual_uniform_zero_mode_premise_certified=False,
                zero_mean_coefficients_remove_scalar_intercept=False)


def conditional_phase_budget(spacing, projected_L2_squared, gradient_L2_squared,
                             scaled_zero_mode_upper, phase_squared=1, time=1, clock_scale=1,
                             beta=None):
    phase, time = nonnegative(phase_squared), positive(time)
    regularizer = previous.regularizer_budget()
    beta = regularizer['beta'] if beta is None else positive(beta)
    if beta > 1:
        raise ValueError('The retained energy fraction must be at most one')
    raw = physical_third_moment_budget(spacing, projected_L2_squared, gradient_L2_squared,
                                      scaled_zero_mode_upper, clock_scale)
    q, m3 = raw['exact_raw_physical_energy'], raw['physical_third_moment_upper']
    r = previous.sqrt_interval(m3/q)['upper'] if q else F(1)
    result = previous.phase_moment_budget(phase*q, r, beta, time)
    result.update(raw_moment_budget=raw, conditional_phase_moment_scale=r,
                  retention_lower_bound_premise=beta,
                  retention_bound_follows_from_default_regularizer=(beta <= regularizer['beta']),
                  conditional_on_actual_uniform_scaled_zero_mode_bound=True,
                  premise_disproved_or_proved_by_this_certificate=False,
                  actual_unconditional_continuum_survival=False)
    return result


def window_scalar(value):
    x = nonnegative(value)
    denominator = 3*x*x+10*x+3
    value = 16*x/denominator
    b = x/(1+x)
    sos = 9*x*x*(x-4)**2+4*x*(x-F(13, 4))**2+F(71, 4)*x+9
    assert denominator-16*x == 3*(x-1)**2
    assert denominator**2-128*x**3 == sos >= 0
    assert 0 <= value <= 1 and x*value*value <= 2
    if x:
        assert value == 16/x/(3/x**2+10/x+3)
        assert value*value/x <= 2
        assert b**3/value <= (3*x+b*b)/16
        assert b**4/value <= 3*x/16
    return dict(x=x, filter=value, resolvent_response=b,
                resolvent_difference=2*(1/(1+x/3)-1/(1+3*x)),
                energy_SOS=sos, Holder_remainder=7*x*x+3*x,
                Cauchy_remainder=2*x**3+15*x*x+12*x+3,
                reciprocal_symmetry=True)


def rational_window_transfer(time, resolvent_lower, source_energy_upper, source_mass_upper=1):
    time, mass = positive(time), positive(source_mass_upper)
    R, energy = nonnegative(resolvent_lower), nonnegative(source_energy_upper)
    if R > mass or 4*R > time*energy:
        raise ValueError('The proposed resolvent floor violates its mass or energy bound')
    W = 16*R*R/(3*time*energy) if energy else F(0)
    Q = 256*R**3/(3*time*energy+R)**2 if R else F(0)
    basic = W*W/mass
    if max(Q, basic) > mass:
        raise ValueError('The proposed source bounds imply a filtered norm larger than the source mass')
    return dict(time=time, resolvent_norm_squared_lower=R,
                source_energy_upper=energy, source_mass_upper=mass,
                filter_expectation_lower=W, filtered_norm_squared_lower=Q,
                basic_Cauchy_filtered_norm_squared_lower=basic,
                filtered_energy_upper=2*mass/time,
                filtered_inverse_energy_upper=2*time*mass,
                positive_input_resolvent_premise_required=True)


def window_band_budget(time, filtered_norm_lower, source_mass_upper=1, ratio=32):
    time, mass, ratio = positive(time), positive(source_mass_upper), positive(ratio)
    lower = nonnegative(filtered_norm_lower)
    if ratio <= 1 or lower > mass:
        raise ValueError('A reciprocal window ratio >1 and a feasible filtered norm are required')
    tail = window_scalar(ratio)['filter']**2
    return dict(spectral_window=dict(lower=1/(ratio*time), upper=ratio/time),
                filter_squared_outside_upper=tail,
                original_source_spectral_mass_lower=max(F(0), lower-mass*tail),
                original_source_mass_upper=mass,
                bound_on_original_source_measure=True)


def phase_window_budget(theta, spacing, projected_L2_squared, phase_squared=1,
                        time=1, clock_scale=1, eta=previous.DEFAULT_ETA):
    phase = previous.phase_energy_budget(theta, spacing, projected_L2_squared,
                                         phase_squared, time, clock_scale, eta)
    e, r, beta, root = (phase[key] for key in ('energy_budget', 'moment_scale', 'beta', 'beta_sqrt_upper'))
    time, scale, spacing = positive(time), positive(clock_scale), positive(spacing)
    D = 1+2*beta*time*e+beta*time*time*e*e+root*r*time*time*e
    R = beta*beta*time*time*e*e/D
    Q = 256*time**4*e**4*beta**4/(D*(3*D+time*e*beta)**2)
    loose = rational_window_transfer(time, R, e)
    assert Q >= loose['filtered_norm_squared_lower'] and 0 <= Q <= 1
    B = phase['regularizer']['second_moment_constant']
    coefficient = 256*beta**4*e/(9*root**3*time*time*B**3*scale**3)
    return dict(actual_phase_budget=phase, resolvent_correlated_denominator=D,
                resolvent_norm_squared_lower=R, filtered_norm_squared_lower=Q,
                uncorrelated_transfer_lower=loose['filtered_norm_squared_lower'],
                correlated_retention_quotient_monotonicity_used=True,
                filtered_energy_upper=2/time, filtered_inverse_energy_upper=2*time,
                normalized_a3_lower_coefficient=Q/spacing**3,
                fixed_data_a3_asymptotic_coefficient_lower=coefficient,
                uniformly_positive_continuum_lower_bound=False)


@lru_cache(maxsize=None)
def actual_cube_window():
    anchor = stage79.actual_cube_witness()['rational_hand_witness']['signed_variational_value']
    transfer = rational_window_transfer(1, anchor, 1)
    lower = transfer['filtered_norm_squared_lower']
    assert lower == F(61632531115632, 1063322424390625) > F(36, 625)
    band = window_band_budget(1, lower)
    assert band['filter_squared_outside_upper'] == F(512, 3395)**2
    assert band['original_source_spectral_mass_lower'] > F(1, 32)
    return dict(side=4, theta=F(1, 48), source='(W_p-mu(W_p))h',
                inherited_resolvent_lower=anchor, rational_window_transfer=transfer,
                actual_original_source_band=band,
                bare_energy_units=True, physical_band='[kappa/32,32*kappa]',
                hidden_space_untruncated=True, uniform_continuum_survival_asserted=False)


@lru_cache(maxsize=None)
def all_controls():
    shape = (2, 2, 4)
    cells = tuple(product(*(range(length) for length in shape)))
    fixtures = []
    for signs in ((0, 0, 0), (1, 0, 0), (1, 1, 1)):
        coefficients = [F((-1)**sum(i*j for i, j in zip(signs, cell))) for cell in cells]
        row = coefficient_gradient_identity(shape, coefficients, 9, (F(-1), F(1, 2), F(-2)))
        symbol = third_moment_symbol(9, (-1, F(1, 2), -2), [(-1)**i for i in signs])
        assert row['matrix_quadratic_form'] == row['coefficient_squared_norm']*symbol['symbol']
        row['rational_Fourier_symbol'] = symbol
        fixtures.append(row)
    windows = [phase_window_budget(25, F(1, n), 1) for n in (4, 16, 256, 65536)]
    for row in windows:
        assert 0 < row['normalized_a3_lower_coefficient'] < row['fixed_data_a3_asymptotic_coefficient_lower']
    return dict(native_Lie_connection_and_Ricci=lie_ricci_control(),
                actual_Ward_response_budgets=[ward_response_budget(theta, F(1, 8)) for theta in (25, 256, 65536)],
                actual_native_third_moment_geometry=[native_third_moment_geometry(side) for side in (4, 6, 8)],
                exact_circulant_controls=fixtures,
                noncommuting_moment_compression=matrix_moment_compression((1, 4, 9, 16),
                    ((F(1, 2), F(7, 10)), (F(1, 4), -F(1, 20)),
                     (F(1, 6), F(1, 30)), (F(1, 8), -F(7, 40)))),
                conditional_uniform_phase_controls=[conditional_phase_budget(F(1, n), 1, 2, 3) for n in (4, 16, 256)],
                rational_window_scalar_controls=[window_scalar(x) for x in (0, F(1, 32), F(1, 3), 1, 3, 32)],
                actual_smooth_phase_window_budgets=windows,
                actual_four_cube_window=actual_cube_window())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 302 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM81 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM81 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM80 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM81_LOCAL_THIRD_MOMENT_VACUUM_RESPONSE',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_LOCAL_RESPONSE_AND_FINITE_WINDOW__UNIFORM_ZERO_MODE_AND_CONTINUUM_SURVIVAL_OPEN',
        proof_scope='Actual Ward and local third-moment identities, native geometry, conditional phase transfer and actual finite window',
        open='Uniform scaled zero-momentum control, nonzero local continuum survival, regular Yang--Mills fields, covariance and mass gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM81_LOCAL_THIRD_MOMENT_VACUUM_RESPONSE':
        raise ValueError('YM81 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM81 fresh certificate/pin mismatch')
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
    print('YM81 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], actual_cube_band_mass_lower='>1/32',
                          uniform_continuum_survival=False)))


if __name__ == '__main__':
    main()
