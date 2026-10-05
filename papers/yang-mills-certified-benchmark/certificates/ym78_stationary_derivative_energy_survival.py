"""YM78: stationary derivative fields, exact local energy and survival tests.

Python 3.12 only. Rational finite replays support the written actual-model
proofs; they do not establish a nonzero local Yang--Mills continuum limit.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym77_physical_smearing_continuum_survival as previous

stage76 = previous.previous
stage75 = stage76.previous
native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative, encode = previous.positive, previous.nonnegative, previous.encode
model = previous.model
RESULT = HERE/'YM78_RESULT.json'
PIN = HERE/'EXPECTED_YM78.sha256'
SOURCES = HERE/'YM78_SOURCE_MANIFEST.json'
BASE_TREE = '8f111b8b053947235927051427dac43b6b158a06'
PARENT_HASH = '6de6268090c349b97f9e43ddb56ceeeb7a1629505357233b20bf8c2a49d2d106'
MODEL = dict(actual_stationary_tempered_space_time_derivative_limit=True,
             derivative_reflection_uses_odd_time_parity=True,
             actual_sparse_local_smear_has_exact_physical_energy=True,
             actual_finite_regulator_reflected_norm_lower_bound=True,
             fixed_physical_frequency_contrast_is_a_survival_witness=True,
             fixed_frequency_survival_window_certified=False,
             nonzero_local_reflected_continuum_state=False,
             full_four_dimensional_Euclidean_covariance=False,
             local_Yang_Mills_field_identification=False,
             positive_continuum_mass_gap=False,
             local_sources_have_exactly_zero_charged_center_Schur_coupling=True,
             center_neutral_sector_equals_local_continuum_core=False,
             charged_resolvent_poles_may_be_ignored_in_full_inverse=False,
             actual_neutral_four_cube_cut_dimension_193=True,
             neutral_restriction_improves_previous_numeric_gap_interval=False,
             finite_replays_are_universal_proof=False,
             spectral_fixtures_are_interacting_Yang_Mills_vacua=False)


@lru_cache(maxsize=None)
def sparse_geometry(side):
    native.width(side, 4)
    if side % 2:
        raise ValueError('The sparse periodic construction requires even side length')
    graph, qs, _, _ = stage76.geometry_data(side)
    selected, used, nonpole_count = [], set(), 0
    for face in graph['faces']:
        root, first_axis = graph['edges'][face[0][0]]
        second_axis = graph['edges'][face[1][0]][1]
        if (first_axis, second_axis) != (0, 1) or root[0] % 2 or root[1] % 2:
            continue
        edges = {edge for edge, _ in face}
        assert used.isdisjoint(edges)
        used.update(edges)
        local_qs = tuple(qs[edge] for edge, _ in face)
        local_word = tuple((i, sign) for i, (_, sign) in enumerate(face))
        trace, gradient = stage76.face_gradient(local_qs, local_word)
        trace_energy = sum(map(stage76.norm2, gradient.values()), F(0))
        assert trace_energy == 1-trace*trace
        if abs(trace) != 1:
            assert trace_energy/(1-trace*trace) == 1
            nonpole_count += 1
        selected.append(root)
    assert len(selected) == side**3//4 and len(used) == side**3
    return dict(side=side, orientation=[0, 1], root_parity='even x, even y, arbitrary z',
                selected_plaquettes=len(selected), occupied_links=len(used),
                density_among_vertices=F(1, 4), rectangular_cell_sides_in_lattice_units=[2, 2, 1],
                nonpole_native_radial_diagonals_checked=nonpole_count,
                distinct_angle_gradients_have_disjoint_link_support=True,
                radial_Dirichlet_Gram_is_identity_almost_everywhere=True,
                interacting_vacuum_independence_assumed=False)


def polynomial_gradient_norm(poly):
    total = F(0)
    for axis in range(3):
        derivative = {}
        for powers, coefficient in poly.items():
            if powers[axis]:
                reduced = list(powers)
                reduced[axis] -= 1
                derivative[tuple(reduced)] = coefficient*powers[axis]
        total += previous.integrate_box(previous.polynomial_square(derivative), (0, 0, 0), (1, 1, 1))
    return total


def sparse_cell_projection(poly, cells_per_xy):
    """Exact 2a by 2a by a cell averages on [0,1]^3, with a=1/(2n)."""
    poly = previous.polynomial(poly)
    native.width(cells_per_xy)
    n, averages = cells_per_xy, []
    spacing = F(1, 2*n)
    sizes = (n, n, 2*n)
    volume = 4*spacing**3
    for index in product(*(range(size) for size in sizes)):
        lower = tuple(F(i, size) for i, size in zip(index, sizes))
        upper = tuple(F(i+1, size) for i, size in zip(index, sizes))
        averages.append(previous.integrate_box(poly, lower, upper)/volume)
    projected = volume*sum((value*value for value in averages), F(0))
    full = previous.integrate_box(previous.polynomial_square(poly), (0, 0, 0), (1, 1, 1))
    gradient = polynomial_gradient_norm(poly)
    floor = max(F(0), full-F(4, 9)*spacing*spacing*gradient)
    assert floor <= projected <= full
    return dict(cells_per_xy=n, spacing=spacing, cell_count=2*n**3, cell_volume=volume,
                cell_average_square_sum=sum((value*value for value in averages), F(0)),
                projected_L2_squared=projected, continuum_L2_squared=full,
                continuum_gradient_L2_squared=gradient,
                continuum_L1_squared_upper=full, L1_bound_uses_unit_volume_Cauchy_Schwarz=True,
                exact_projection_error=full-projected,
                conservative_Poincare_lower=floor, pi_squared_lower_used=F(9))


def smear_moment_budget(theta, spacing, L1_squared_upper, L2_squared,
                        gradient_L2_squared, clock_scale=1, exact_projected_L2_squared=None):
    theta, spacing, scale = positive(theta), positive(spacing), positive(clock_scale)
    if theta < 25:
        raise ValueError('The actual vacuum moment bounds require theta >= 25')
    l1, l2, gradient = map(nonnegative, (L1_squared_upper, L2_squared, gradient_L2_squared))
    lower = scale*max(F(0), l2-F(4, 9)*spacing*spacing*gradient)
    upper = scale*l2
    if exact_projected_L2_squared is not None:
        projected = nonnegative(exact_projected_L2_squared)
        if not lower <= scale*projected <= upper:
            raise ValueError('The proposed exact projection violates Jensen or Poincare')
        lower = upper = scale*projected
    m0, m2 = F(45, 16)*l1/spacing**2, 10*scale*scale*l1/spacing**4
    if lower*lower > m0*m2:
        raise ValueError('The declared source moments violate spectral Cauchy-Schwarz')
    return dict(theta=theta, spacing=spacing, sqrt_cE_cB=scale,
                smear_coefficient_fourth_power=16*theta*spacing**8,
                physical_first_moment=dict(lower=lower, upper=upper),
                physical_source_norm_squared_upper=m0, physical_second_moment_upper=m2,
                pi_squared_upper_used=F(10), actual_sparse_exact_energy_identity=True,
                bounds_are_uniform_in_volume=True,
                positive_continuum_reflected_norm_inferred=False)


def reflected_norm_lower(time, q_lower, m0_upper, m2_upper, q_upper=None):
    time = positive(time)
    q, m0, m2 = map(nonnegative, (q_lower, m0_upper, m2_upper))
    upper = q if q_upper is None else nonnegative(q_upper)
    if q > upper or q*q > m0*m2:
        raise ValueError('Inconsistent first, zeroth or second spectral moments')
    denominator = m0+time*upper+time*time*m2
    lower = time*time*q*q/denominator if denominator else F(0)
    return dict(filter_time=time, source_energy_lower=q, source_energy_upper=upper,
                source_norm_upper=m0, source_second_moment_upper=m2,
                certified_reflected_norm_lower=lower, exact_denominator=denominator,
                fixed_regulator_lower_bound_only=True)


def scalar_filter_inequalities(value):
    x = nonnegative(value)
    denominator = 1+x+x*x/2
    quadratic_minorant = x*(1+x/2)/denominator
    square_remainder = (1+x+x*x)*(1+x/2)**2-denominator**2
    assert square_remainder == x*x*(1+x)/4 >= 0
    quartic = x**4-2*x**3+x*x-2*x+4
    if x <= 2:
        decomposition = x*x*(x-1)**2+2*(2-x)
    else:
        y = x-2
        decomposition = y**4+6*y**3+13*y*y+10*y+4
    assert quartic == decomposition >= 0
    kernel = 6*x*x/((x*x+1)*(x*x+4))
    assert 0 <= kernel <= F(2, 3)
    # e^x >= 1+x+x^2/2 proves the first comparison, e^x >=1+x the second.
    assert quadratic_minorant**2 >= x*x/(1+x+x*x)
    assert x*x/(1+x)**2 >= x*kernel/3
    return dict(x=x, rational_exponential_minorant=quadratic_minorant,
                quadratic_comparison_remainder=square_remainder,
                frequency_contrast_kernel=kernel, contrast_quartic=quartic,
                nonnegative_quartic_decomposition=decomposition,
                reflected_norm_scalar_lower=x*x/(1+x+x*x),
                contrast_reflected_norm_scalar_lower=x*kernel/3)


def derivative_moment_budget(energy_bound, temporal_L2_squared, order):
    energy_bound, time_norm = nonnegative(energy_bound), nonnegative(temporal_L2_squared)
    native.width(order)
    A = energy_bound*time_norm
    return dict(spatial_energy_bound=energy_bound, temporal_test_L2_squared=time_norm,
                mgf_quadratic_coefficient=A, variance_upper=2*A,
                half_moment_order=order, even_moment_upper=2*factorial(order)*(4*A)**order,
                derivative_reflection_is_odd=True,
                uniform_pointwise_Dirichlet_bound_required=True)


def spectral_measure(atoms, weights):
    atoms, weights = tuple(map(nonnegative, atoms)), tuple(map(nonnegative, weights))
    if not atoms or len(atoms) != len(weights):
        raise ValueError('A finite nonnegative measure with matching dimensions is required')
    return atoms, weights


def derivative_kernel(atoms, weights, frequency):
    atoms, weights = spectral_measure(atoms, weights)
    s = positive(frequency)
    return sum((2*w*E*s*s/(E*E+s*s) for E, w in zip(atoms, weights)), F(0))


def derivative_spectral_control(atoms, weights, frequency=1):
    atoms, weights = spectral_measure(atoms, weights)
    s = positive(frequency)
    m0 = sum(weights)
    q = sum((E*w for E, w in zip(atoms, weights)), F(0))
    m2 = sum((E*E*w for E, w in zip(atoms, weights)), F(0))
    K, K2 = derivative_kernel(atoms, weights, s), derivative_kernel(atoms, weights, 2*s)
    contrast = K2-K
    kernel_sum = sum((E*w*scalar_filter_inequalities(E/s)['frequency_contrast_kernel']
                      for E, w in zip(atoms, weights)), F(0))
    assert contrast == kernel_sum and 0 <= contrast <= F(2, 3)*q
    lower = reflected_norm_lower(1/s, q, m0, m2)['certified_reflected_norm_lower']
    direct_lower = direct_upper = F(0)
    for E, w in zip(atoms, weights):
        exp = stage76.exp_negative_interval(E/s)
        direct_lower += w*(1-exp['upper'])**2
        direct_upper += w*min(F(1), (E/s)**2, (1-exp['lower'])**2)
    certified = max(lower, contrast/(3*s), direct_lower)
    assert 0 <= certified <= direct_upper
    return dict(atoms=atoms, weights=weights, frequency=s, source_norm_squared=m0,
                source_energy=q, source_second_moment=m2, derivative_kernel=K,
                derivative_kernel_at_double_frequency=K2, frequency_contrast=contrast,
                contrast_reflected_norm_lower=contrast/(3*s), moment_reflected_norm_lower=lower,
                reflected_norm_interval=dict(lower=certified, upper=direct_upper),
                algebraic_fixture_not_interacting_vacuum=True)


def dyadic_contrast_control(atoms, weights, base_frequency=1, levels=6):
    atoms, weights = spectral_measure(atoms, weights)
    base = positive(base_frequency)
    native.width(levels)
    terms = []
    for exponent in range(-levels, levels):
        frequency = base*F(2)**exponent
        contrast = derivative_kernel(atoms, weights, 2*frequency)-derivative_kernel(atoms, weights, frequency)
        assert contrast >= 0
        terms.append(dict(exponent=exponent, frequency=frequency, contrast=contrast))
    q = sum((E*w for E, w in zip(atoms, weights)), F(0))
    low = derivative_kernel(atoms, weights, base*F(2)**(-levels))
    high = derivative_kernel(atoms, weights, base*F(2)**levels)
    total = sum((row['contrast'] for row in terms), F(0))
    assert total == high-low and total <= 2*q
    return dict(base_frequency=base, levels=levels, terms=terms,
                exact_sum=total, endpoint_difference=high-low,
                all_frequency_sum=2*q, missing_frequency_tails=low+2*q-high,
                fixed_frequency_window_survival_inferred=False)


@lru_cache(maxsize=None)
def neutral_cube_controls():
    native_geometry = previous.center_geometry(4)
    theta, ground_z, excited_z = F(1, 48), -F(712209, 10**8), F(48376192, 10**8)
    rows = []
    for energy in (ground_z, excited_z):
        full = stage75.cube_schur_inertia(theta, energy)
        transverse = stage75.scalar_inertia([full['transverse_pivot']])
        neutral = {key: full['uniform_inertia'][key]+191*transverse[key] for key in transverse}
        assert sum(neutral.values()) == 193
        rows.append(dict(energy=energy, neutral_inertia=neutral,
                         uniform_pivots=full['uniform_pivots'], transverse_pivot=full['transverse_pivot']))
    assert rows[0]['neutral_inertia'] == dict(negative=0, zero=0, positive=193)
    assert rows[1]['neutral_inertia'] == dict(negative=1, zero=0, positive=192)
    actual = stage75.actual_cube_controls()
    graph = model.lattice(4)
    index = {edge: i for i, edge in enumerate(graph['edges'])}
    counts = {}
    for axis in range(3):
        other = [j for j in range(3) if j != axis]
        for transverse in product(range(4), repeat=2):
            word = []
            for k in range(4):
                vertex = [0, 0, 0]
                vertex[axis] = k
                for j, coordinate in zip(other, transverse):
                    vertex[j] = coordinate
                word.append(index[tuple(vertex), axis])
            charge = tuple(sum(graph['edges'][e][1] == j and graph['edges'][e][0][j] == 3
                               for e in word) % 2 for j in range(3))
            counts[charge] = counts.get(charge, 0)+1
    assert sorted(counts.values()) == [16, 16, 16]
    return dict(theta=theta, full_R3_dimension=241, neutral_R3_dimension=193,
                removed_charged_winding_dimension=48,
                winding_character_multiplicities=[dict(charge=list(k), count=v) for k, v in sorted(counts.items())],
                actual_neutral_hidden_free_floor=F(9, 2),
                free_floor_is_attained_by_contractible_six_edge_rectangle=True,
                interacting_hidden_floor=F(1, 2), exact_neutral_thresholds=rows,
                certified_neutral_gap=actual['certified_bare_gap'],
                gap_units='bare electric Hamiltonian units',
                previous_full_gap_same_interval=True,
                native_center_seam_check=native_geometry['all_plaquettes_have_even_seam_parity'],
                inverse_must_be_taken_on_neutral_hidden_complement=True,
                neutral_sector_identified_with_local_continuum_core=False)


@lru_cache(maxsize=None)
def all_controls():
    pi = stage76.spatial.pi_interval()
    assert 9 < pi.lo*pi.lo <= pi.hi*pi.hi < 10
    polys = ({(0, 0, 0): 1}, {(1, 0, 0): 1},
             {(2, 0, 0): 1, (0, 1, 1): -1, (0, 0, 0): F(1, 3)})
    cells = [sparse_cell_projection(poly, n) for poly in polys for n in (2, 3, 4)]
    budgets = []
    for row in cells:
        budget = smear_moment_budget(256, row['spacing'], row['continuum_L1_squared_upper'],
            row['continuum_L2_squared'], row['continuum_gradient_L2_squared'],
            exact_projected_L2_squared=row['projected_L2_squared'])
        q = budget['physical_first_moment']['lower']
        budgets.append(dict(moments=budget, reflected_lower=reflected_norm_lower(1, q,
            budget['physical_source_norm_squared_upper'], budget['physical_second_moment_upper'])))
    ultraviolet = [derivative_spectral_control((n,), (F(1, n),)) for n in (2, 4, 32, 256)]
    infrared = [derivative_spectral_control((F(1, n),), (n,)) for n in (2, 4, 32, 256)]
    return dict(actual_sparse_geometry=[sparse_geometry(n) for n in (4, 6)],
                exact_rectangular_cell_projections=cells, actual_moment_and_reflection_budgets=budgets,
                scalar_positive_inequalities=[scalar_filter_inequalities(x) for x in (0, F(1, 100), F(1, 2), 1, 2, 3, 100)],
                derivative_test_moments=[derivative_moment_budget(2, F(1, 3), order) for order in (1, 2, 3)],
                ultraviolet_energy_escape=ultraviolet, infrared_energy_escape=infrared,
                surviving_fixed_energy_fixture=derivative_spectral_control((1,), (1,)),
                exact_dyadic_telescoping=[dyadic_contrast_control((F(1, 8), 1, 32), (2, F(1, 3), F(1, 64)), levels=n)
                                          for n in (2, 4, 8)],
                actual_neutral_sector=neutral_cube_controls())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 281 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM78 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM78 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM77 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM78_STATIONARY_DERIVATIVE_ENERGY_SURVIVAL',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_DERIVATIVE_FIELD_AND_EXACT_ENERGY__LOCAL_REFLECTED_SURVIVAL_OPEN',
        proof_scope='Written actual-model theorems with exact finite native geometry and rational spectral controls',
        open='Nonzero local reflected limit, full Euclidean covariance, Yang--Mills identification and positive local gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM78_STATIONARY_DERIVATIVE_ENERGY_SURVIVAL':
        raise ValueError('YM78 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM78 fresh certificate/pin mismatch')
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
    print('YM78 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], neutral_R3_dimension=193,
                          stationary_derivative_field=True, nonzero_local_reflected_limit=False)))


if __name__ == '__main__':
    main()
