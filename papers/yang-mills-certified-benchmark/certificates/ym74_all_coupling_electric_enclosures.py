"""YM74: all-coupling electric cutoffs, vacuum flatness and finite-volume gaps.

Python 3.12 only. The scalar budgets apply to the actual finite SU(2) model;
labelled rational matrix fixtures test algebra and are not its vacuum.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from math import comb, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym73_physical_mass_memory_metric as previous

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative = previous.positive, previous.nonnegative
matrix, sub, symmetric = previous.matrix, previous.sub, previous.symmetric
model = previous.adapted.model
encode = previous.previous.spatial.encode
GRID = previous.GRID
RESULT = HERE/'YM74_RESULT.json'
PIN = HERE/'EXPECTED_YM74.sha256'
SOURCES = HERE/'YM74_SOURCE_PINS.json'
BASE_TREE = 'ad47015dc1abbba36fb338f5a5f531074e3f54a1'

MODEL = dict(actual_all_coupling_vacuum_energy_and_flatness_bound=True,
             gauge_noninvariant_trial_used_only_for_ground_energy=True,
             all_coupling_hidden_floor=True, whole_low_energy_capture=True,
             fixed_regulator_arbitrary_precision=True,
             electric_cutoff_is_gauge_invariant_and_finite_rank=True,
             hidden_space_is_not_truncated=True,
             actual_periodic_four_cube_241_dimensional_cut=True,
             winding_four_cycles_retained=True,
             actual_finite_cube_gap_beyond_YM62_window=True,
             uniform_weak_coupling_positive_gap=False,
             continuum_R4_QFT=False, physical_mass_prediction=False,
             finite_fixture_is_Yang_Mills_vacuum=False,
             coverage_implies_nonvacuum_low_energy_states_exist=False,
             plaquette_flatness_controls_continuum_curvature=False)


def ceil_fraction(value):
    value = exact(value)
    return -((-value.numerator)//value.denominator)


def sqrt_interval(value):
    value = nonnegative(value)
    integer = isqrt(value.numerator*GRID*GRID//value.denominator)
    lower = F(integer, GRID)
    upper = lower if lower*lower == value else F(integer+1, GRID)
    return dict(lower=lower, upper=upper)


def ceil_sqrt(value):
    value = nonnegative(value)
    integer = isqrt(value.numerator//value.denominator)
    return integer if integer*integer == value else integer+1


def product_trial(index):
    """Exact normalized moments of f_n(q)=(1+Re q)^n; no fitted state."""
    native.width(index, 0)
    mean = F(2*index, 2*index+3)
    electric = F(3*index*index, 4*(4*index+1))
    return dict(index=index, link_scalar_mean=mean, link_electric_energy=electric,
                plaquette_mean=mean**4, plaquette_deficiency=1-mean**4,
                beta_alpha=F(4*index+3, 2), beta_beta=F(3, 2))


def vacuum_budget(theta, plaquettes=1, index=None):
    """Periodic cubic model (number of links equals number of plaquettes)."""
    theta = nonnegative(theta)
    native.width(plaquettes)
    if index is None:
        index = ceil_sqrt(32*theta)
    trial = product_trial(index)
    density = trial['link_electric_energy']+theta*trial['plaquette_deficiency']
    deficit = plaquettes*min(theta, density)
    flat = min(F(1), density/theta) if theta else F(1)
    asymptotic = (F(3, 2)*sqrt_interval(2*theta)['upper']+F(3, 16))
    return dict(theta=theta, plaquettes=plaquettes, trial=trial,
                trial_energy_upper=plaquettes*(density-theta),
                ground_energy_lower=-plaquettes*theta,
                ground_energy_upper=deficit-plaquettes*theta,
                actual_shifted_vacuum_deficit_upper=deficit,
                actual_total_electric_energy_upper=deficit,
                actual_link_electric_energy_upper=deficit/plaquettes,
                actual_plaquette_deficiency_upper=flat,
                actual_holonomy_Frobenius_squared_upper=4*flat,
                optimized_density_asymptotic_upper=asymptotic,
                optimized_index=index == ceil_sqrt(32*theta),
                volume_uniform_per_plaquette=True,
                applies_to_actual_vacuum_not_trial_expectation=True)


def electric_cutoff_budget(theta, plaquettes, cutoff, physical_band=0, clock=1,
                           next_electric_lower=None):
    """Bound Q_R H Q_R and the whole physical band.

    An optional next_electric_lower must be independently certified; the
    default cutoff itself is automatically a valid (possibly nonsharp) floor.
    """
    cutoff, physical_band, clock = nonnegative(cutoff), nonnegative(physical_band), positive(clock)
    next_floor = cutoff if next_electric_lower is None else positive(next_electric_lower)
    if next_floor < cutoff or not next_floor:
        raise ValueError('A positive next-electric lower bound >= cutoff is required')
    cell = vacuum_budget(theta, plaquettes)
    deficit = cell['actual_shifted_vacuum_deficit_upper']
    hidden = next_floor-deficit
    leakage = min(F(1), (deficit+physical_band/clock)/next_floor)
    return dict(cutoff=cutoff, next_electric_lower=next_floor,
                vacuum_deficit_upper=deficit, hidden_bare_floor=hidden,
                hidden_physical_floor=clock*hidden,
                positive_hidden_floor_certified=hidden > 0,
                entire_low_energy_subspace_leakage_squared_upper=leakage,
                vacuum_leakage_squared_upper=min(F(1), deficit/next_floor),
                physical_band=physical_band, clock=clock,
                nonvacuum_low_energy_states_asserted=False)


def combined_cutoff(theta, plaquettes, clock, growth, reference_energy=1,
                    physical_band=1, precision=None):
    theta, clock, reference_energy = nonnegative(theta), positive(clock), positive(reference_energy)
    growth, physical_band = positive(growth), nonnegative(physical_band)
    native.width(plaquettes)
    if growth < 1:
        raise ValueError('Growth must be at least one')
    deficit = vacuum_budget(theta, plaquettes)['actual_shifted_vacuum_deficit_upper']
    capture = max(F(3), growth*(deficit+reference_energy/clock))
    chosen = F(ceil_fraction(capture))
    precision_cut = None
    if precision is not None:
        precision = positive(precision)
        b = theta*plaquettes
        precision_cut = F(4+ceil_fraction(2*b+2*clock*b*b/precision))
        chosen = max(chosen, precision_cut)
    actual = electric_cutoff_budget(theta, plaquettes, chosen, physical_band, clock)
    coarse_leakage = min(F(1), max(F(1), physical_band/reference_energy)/growth)
    assert actual['hidden_physical_floor'] >= growth*reference_energy
    assert actual['entire_low_energy_subspace_leakage_squared_upper'] <= coarse_leakage
    width = None
    if precision is not None:
        width = 2*clock*(theta*plaquettes)**2/(chosen-3-2*theta*plaquettes)
        assert width <= precision
    return dict(cutoff=chosen, capture_cutoff=capture, precision_cutoff=precision_cut,
                growth=growth, reference_energy=reference_energy,
                guaranteed_hidden_physical_floor=growth*reference_energy,
                guaranteed_band_leakage_squared_upper=coarse_leakage,
                physical_gap_total_width_upper=width, requested_precision=precision,
                actual_budget=actual)


def finite_shell(cutoff):
    cutoff = nonnegative(cutoff)
    ceiling = ceil_sqrt(cutoff)
    return dict(retained_electric_cutoff=cutoff,
                exact_coupling_shell_upper=cutoff+2*ceiling+3,
                ceil_sqrt_cutoff=ceiling,
                four_fundamental_edge_multiplications=True)


def representation_rank_bound(cutoff, edges):
    cutoff = nonnegative(cutoff)
    native.width(edges)
    twice_spin = isqrt((4*cutoff+1).numerator//(4*cutoff+1).denominator)-1
    twice_spin = max(0, twice_spin)
    dimension = (twice_spin+1)*(twice_spin+2)*(2*twice_spin+3)//6
    return dict(maximum_twice_edge_spin=twice_spin,
                unrestricted_single_edge_dimension=dimension,
                physical_retained_rank_upper=dimension**edges,
                gauge_projection_can_only_reduce_rank=True)


def ritz_root_interval(retained_energy, hidden_floor, coupling_squared):
    """Outward f(c,d,q)=(c+d-sqrt((c-d)^2+4q))/2, q=||B||^2 bound."""
    c, d, q = exact(retained_energy), exact(hidden_floor), nonnegative(coupling_squared)
    radius = sqrt_interval((c-d)**2+4*q)
    return dict(lower=(c+d-radius['upper'])/2,
                upper=(c+d-radius['lower'])/2)


def ritz_eigenvalue_interval(retained_lower, retained_upper, hidden_floor, coupling_squared):
    lo, hi = exact(retained_lower), exact(retained_upper)
    if hi < lo:
        raise ValueError('Reversed retained eigenvalue interval')
    lower = ritz_root_interval(lo, hidden_floor, coupling_squared)['lower']
    return dict(lower=lower, upper=hi)


def ritz_gap_interval(ground, excited, hidden_floor, coupling_squared, clock=1):
    clock = positive(clock)
    if len(ground) != 2 or len(excited) != 2:
        raise ValueError('Two endpoint pairs required')
    ground, excited = tuple(map(exact, ground)), tuple(map(exact, excited))
    if ground[0] > ground[1] or excited[0] > excited[1]:
        raise ValueError('Reversed retained eigenvalue interval')
    if excited[1] < ground[0]:
        raise ValueError('First excited Ritz value cannot precede the ground Ritz value')
    first = ritz_eigenvalue_interval(*ground, hidden_floor, coupling_squared)
    second = ritz_eigenvalue_interval(*excited, hidden_floor, coupling_squared)
    lower, upper = max(F(0), second['lower']-first['upper']), second['upper']-first['lower']
    if upper < lower:
        raise ValueError('Inconsistent spectral enclosures')
    return dict(lower=clock*lower, upper=clock*upper, width=clock*(upper-lower),
                ground_eigenvalue=first, first_excited_eigenvalue=second,
                clock=clock, positive_gap_certified=lower > 0)


def ritz_width_budget(potential_bound, cutoff, clock=1, arithmetic_width=0):
    b, cutoff, clock = nonnegative(potential_bound), nonnegative(cutoff), positive(clock)
    arithmetic_width = nonnegative(arithmetic_width)
    margin = cutoff-3-2*b
    if margin <= 0:
        raise ValueError('This universal width bound requires cutoff > 3+2b')
    return dict(physical_total_width_upper=clock*(2*b*b/margin+arithmetic_width),
                retained_Ritz_upper=3+b, electric_separation_margin=margin,
                ground_plus_excited_arithmetic_width=arithmetic_width)


def coupling_parameter_budget(plaquettes, theta_radius, clock=1):
    """Weyl parameter error at a FIXED physical clock, with ||V||<=Np."""
    native.width(plaquettes)
    theta_radius, clock = nonnegative(theta_radius), positive(clock)
    one_energy = clock*plaquettes*theta_radius
    return dict(physical_energy_error_radius=one_energy,
                physical_gap_error_radius=2*one_energy,
                total_gap_interval_width_increase_upper=4*one_energy,
                physical_clock_held_fixed=True)


def two_by_two_eigenvalues(rows):
    rows = symmetric(rows)
    if len(rows) != 2:
        raise ValueError('Two by two matrix required')
    trace = rows[0][0]+rows[1][1]
    radius = sqrt_interval((rows[0][0]-rows[1][1])**2+4*rows[0][1]**2)
    return (dict(lower=(trace-radius['upper'])/2, upper=(trace-radius['lower'])/2),
            dict(lower=(trace+radius['lower'])/2, upper=(trace+radius['upper'])/2))


def four_cube_gap(theta, clock=1):
    theta, clock = exact(theta), positive(clock)
    ground_radius = sqrt_interval(9+192*theta*theta)
    ground = ((3-ground_radius['upper'])/2, (3-ground_radius['lower'])/2)
    hidden = F(9, 2)-192*abs(theta)
    q = F(383, 4)*theta*theta
    result = ritz_gap_interval(ground, (3, 3), hidden, q, clock)
    return dict(theta=theta, periodic_side=4, physical_retained_dimension=241,
                retained_ground_interval=dict(lower=ground[0], upper=ground[1]),
                retained_first_excited_energy=F(3), next_free_electric_floor=F(9, 2),
                unshifted_hidden_floor=hidden, exact_coupling_squared=q,
                gap=result, actual_finite_cube=True,
                evaluated_interacting_vacuum=False, continuum_claim=False)


@lru_cache(maxsize=None)
def geometry_controls():
    cell = model.lattice(4)
    adjacent = {vertex: [] for vertex in cell['vertices']}
    lookup = {}
    for edge, (left, right) in enumerate(cell['ends']):
        assert sum(left) % 2 != sum(right) % 2
        adjacent[left].append(right)
        adjacent[right].append(left)
        lookup[frozenset((left, right))] = edge
    cycles = set()
    for start in cell['vertices']:
        for v1 in adjacent[start]:
            for v2 in adjacent[v1]:
                if v2 == start:
                    continue
                for v3 in adjacent[v2]:
                    if v3 in (start, v1) or start not in adjacent[v3]:
                        continue
                    word = ((start, v1), (v1, v2), (v2, v3), (v3, start))
                    cycles.add(frozenset(lookup[frozenset(pair)] for pair in word))
    plaquettes = {frozenset(edge for edge, sign in face) for face in cell['faces']}
    winding = cycles-plaquettes
    assert len(cycles) == 240 and len(plaquettes) == 192 and len(winding) == 48
    assert all(len(left.intersection(right)) <= 1 for left, right in combinations(cycles, 2))
    pair_xors = {left ^ right for left, right in combinations(plaquettes, 2)}
    winding_xors = {left ^ right for left, right in combinations(winding, 2)}
    mixed_xors = {left ^ right for left in plaquettes for right in winding}
    assert len(pair_xors) == comb(192, 2) == 18336
    assert not pair_xors.intersection(cycles | winding_xors | mixed_xors)
    assert not cycles.intersection(mixed_xors | winding_xors)
    assert all(len({cell['edges'][edge][1] for edge in loop}) == 1 for loop in winding)
    native_plaquette = model.plaquette_controls()
    x = p.COORD[0]
    second = native.phi(p.mul(x, x), 1)
    fourth = native.phi(p.mul(p.mul(x, x), p.mul(x, x)), 1)
    assert second == F(1, 4) and fourth == F(1, 8)
    assert native_plaquette['variance'] == second
    # The parity audit leaves only paired Haar moments in P V^2 P.
    # Subtracting (P V P)^2 removes the vacuum-star return contribution.
    plaquette_diagonal = 4*(fourth+191*second*second)-4*second*second
    plaquette_off_diagonal = 8*second*second-4*second*second
    plaquette_identity = plaquette_diagonal-plaquette_off_diagonal
    winding_diagonal = 4*192*second*second
    largest = max(plaquette_identity+192*plaquette_off_diagonal, winding_diagonal)
    assert largest == F(383, 4)
    # Product of independent Haar links is Haar; W^4 reduces to this one-link moment.
    return dict(vertices=64, links=192, plaquettes=192, simple_four_cycles=240,
                winding_four_cycles=48, physical_R3_dimension=241,
                periodic_four_cube_graph_bipartite=True,
                distinct_plaquette_pair_xors=len(pair_xors),
                maximum_distinct_four_cycle_edge_overlap=1,
                all_four_cycle_and_mixed_pair_parity_collisions_absent=True,
                Haar_W_squared=second, Haar_W_fourth=fourth,
                native_plaquette=native_plaquette,
                next_electric_floor=F(9, 2),
                BBstar_plaquette_identity_coefficient=plaquette_identity,
                BBstar_plaquette_all_ones_coefficient=plaquette_off_diagonal,
                BBstar_winding_identity_coefficient=winding_diagonal,
                BBstar_largest_eigenvalue_over_theta_squared=largest)


@lru_cache(maxsize=None)
def trial_controls():
    rows = []
    # Independently integrate f_n^2 and its electric form using native Haar moments.
    x = p.COORD[0]
    coordinate_powers = [{(0, 0, 0, 0): F(1)}]
    for index in range(16):
        coordinate_powers.append(p.mul(coordinate_powers[-1], x))
    moments = [native.phi(polynomial, 1) for polynomial in coordinate_powers]
    for n in range(8):
        norm = sum(F(comb(2*n, k))*moments[k] for k in range(2*n+1))
        mean = sum(F(comb(2*n, k))*moments[k+1] for k in range(2*n+1))/norm
        energy = F(0) if not n else F(n*n, 4)*sum(
            F(comb(2*n-2, k))*(moments[k]-moments[k+2]) for k in range(2*n-1))/norm
        row = product_trial(n)
        assert mean == row['link_scalar_mean'] and energy == row['link_electric_energy']
        rows.append(dict(index=n, independent_Haar_norm=norm,
                         independent_Haar_mean=mean, independent_Haar_electric=energy))
    regimes = []
    for theta in map(F, (1, 16, 256, 4096, 65536, 10**6)):
        row = vacuum_budget(theta, 192)
        density = row['trial']['link_electric_energy']+theta*row['trial']['plaquette_deficiency']
        assert density <= row['optimized_density_asymptotic_upper']
        assert row['ground_energy_lower'] <= row['ground_energy_upper'] <= 0
        assert vacuum_budget(theta, 3*16**3)['actual_plaquette_deficiency_upper'] == row['actual_plaquette_deficiency_upper']
        regimes.append(row)
    assert all(regimes[i]['actual_plaquette_deficiency_upper'] > regimes[i+1]['actual_plaquette_deficiency_upper']
               for i in range(len(regimes)-1))
    return dict(independent_native_integrations=rows, actual_vacuum_weak_coupling_budgets=regimes,
                nonphysical_trial_not_used_for_excited_Ritz_eigenvalues=True)


@lru_cache(maxsize=None)
def matrix_controls():
    rows = []
    # Rational Householder conjugates have independently known FULL eigenvalues.
    for vector, eigenvalues in (((1, 2, 3, 4), (-2, 1, 5, 9)),
                                ((2, -1, 1, 3), (0, 1, 4, 11)),
                                ((1, 1, -2, 1), (-5, -2, 3, 7))):
        vector = list(map(F, vector))
        denominator = sum(value*value for value in vector)
        outer = [[x*y for y in vector] for x in vector]
        U = a.add(a.identity(4), outer, -2/denominator)
        assert a.multiply(a.transpose(U), U) == a.identity(4)
        H = a.multiply(a.multiply(U, a.diag(eigenvalues)), a.transpose(U))
        C, D, B = sub(H, (0, 1)), sub(H, (2, 3)), sub(H, (0, 1), (2, 3))
        BB = a.multiply(B, a.transpose(B))
        c0, c1 = two_by_two_eigenvalues(C)
        hidden = two_by_two_eigenvalues(D)[0]['lower']
        q = two_by_two_eigenvalues(BB)[1]['upper']
        assert a.psd(a.add(D, a.identity(2), -hidden))
        assert a.psd(a.add(a.scale(a.identity(2), q), BB, -1))
        interval = ritz_gap_interval((c0['lower'], c0['upper']),
                                     (c1['lower'], c1['upper']), hidden, q)
        for j, retained in enumerate((c0, c1)):
            row = ritz_eigenvalue_interval(retained['lower'], retained['upper'], hidden, q)
            assert row['lower'] <= eigenvalues[j] <= row['upper']
        assert interval['lower'] <= eigenvalues[1]-eigenvalues[0] <= interval['upper']
        commutator = a.add(a.multiply(C, BB), a.multiply(BB, C), -1)
        assert any(value for line in commutator for value in line)
        rows.append(dict(exact_full_spectrum=list(eigenvalues), retained_matrix=C,
                         hidden_matrix=D, coupling_matrix=B, coupling_squared_upper=q,
                         nonzero_retained_coupling_commutator=commutator,
                         exact_full_gap=F(eigenvalues[1]-eigenvalues[0]), certified_gap=interval))
    # A retained block entirely ABOVE an independent hidden level needs f=min(c,d).
    invisible = ritz_gap_interval((3, 3), (5, 5), 1, 0)
    assert invisible['lower'] == 0 and invisible['upper'] == 4
    return dict(rational_fixtures_not_Yang_Mills_vacua=rows,
                invisible_hidden_level=invisible)


@lru_cache(maxsize=None)
def cutoff_controls():
    rows = []
    for theta, plaquettes, clock, growth in ((F(1, 8000), 192, F(1), F(4)),
            (F(16), 192, F(1, 4), F(16)),
            (F(256), 1536, F(1, 100), F(64)),
            (F(10**6), 12288, F(1000), F(256))):
        row = combined_cutoff(theta, plaquettes, clock, growth,
                              reference_energy=2, physical_band=3, precision=F(1, 100))
        assert row['physical_gap_total_width_upper'] <= F(1, 100)
        row['theta'], row['plaquettes'], row['clock'] = theta, plaquettes, clock
        rows.append(row)
    shells = [finite_shell(value) for value in (F(0), F(3), F(17, 2), F(64))]
    for row in shells:
        assert (row['ceil_sqrt_cutoff'])**2 >= row['retained_electric_cutoff']
    # This is a rank bound, not enumeration of the displayed enormous spaces.
    ranks = [representation_rank_bound(F(3), edges) for edges in (1, 4, 12)]
    fixed = []
    for cutoff in (F(16), F(64), F(256), F(1024)):
        fixed.append(dict(cutoff=cutoff, **ritz_width_budget(F(2), cutoff, F(3))))
    assert all(fixed[i]['physical_total_width_upper'] > fixed[i+1]['physical_total_width_upper']
               for i in range(len(fixed)-1))
    parameter = coupling_parameter_budget(192, F(1, 10**10), F(3, 2))
    assert parameter['physical_gap_error_radius'] == F(576, 10**10)
    return dict(all_coupling_cutoff_examples=rows, finite_electric_coupling_shells=shells,
                finite_rank_bounds=ranks, fixed_regulator_width_convergence=fixed,
                finite_coupling_parameter_error=parameter,
                no_exponential_size_matrix_was_numerically_diagonalized=True)


@lru_cache(maxsize=None)
def actual_cube_controls():
    geometry_controls()
    rows = [four_cube_gap(F(1, denominator)) for denominator in (8000, 128, 64, 48, 40)]
    chosen = rows[3]['gap']
    assert chosen['lower'] > F(4904, 10000)
    assert chosen['upper'] < F(307874, 100000)
    assert rows[-1]['gap']['lower'] == 0
    assert four_cube_gap(0)['gap']['lower'] == four_cube_gap(0)['gap']['upper'] == 3
    return dict(actual_periodic_cube_gap_intervals=rows,
                theta_one_over_48_positive_bare_gap=True,
                actual_volume_independent_gap_at_this_theta_asserted=False)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 253 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM74 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM74 base tree mismatch')
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
    return encode(dict(certificate_type='YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        native_geometry=geometry_controls(), actual_vacuum_trials=trial_controls(),
        rational_spectral_controls=matrix_controls(), electric_cutoff_controls=cutoff_controls(),
        actual_finite_cube=actual_cube_controls(), model=MODEL,
        claim_status='ACTUAL_ALL_COUPLING_CUTOFF_AND_VACUUM_BOUNDS__CONTINUUM_MASS_NOT_PROVED',
        proof_scope='Written operator proofs with native exact Haar, graph and rational spectral controls',
        open='Uniform weak-coupling positive physical edge, nontrivial local field limits and 4D Yang--Mills construction'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM74_ALL_COUPLING_ELECTRIC_ENCLOSURES':
        raise ValueError('YM74 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM74 fresh certificate/pin mismatch')
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
    print('YM74 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], all_coupling_hidden_floor=True,
                          fixed_regulator_arbitrary_precision=True,
                          actual_four_cube_theta_one_over_48_gap_lower='>0.4904',
                          continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
