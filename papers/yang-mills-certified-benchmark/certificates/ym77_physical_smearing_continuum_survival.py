"""YM77: physical smearing, anchored compactness and explicit survival tests.

Python 3.12 only. Actual center-sector trial states and smearing estimates
are distinguished from finite spectral controls and unproved local survival.
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
import ym76_local_defect_weak_coupling_dynamics as previous

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative, encode = previous.positive, previous.nonnegative, previous.encode
model, norm2 = previous.model, previous.norm2
RESULT = HERE/'YM77_RESULT.json'
PIN = HERE/'EXPECTED_YM77.sha256'
SOURCES = HERE/'YM77_SOURCE_MANIFEST.json'
BASE_TREE = '7fac8100c57a0bec7553a80793dd2af914aa59df'
PARENT_HASH = '1e96e5526edb84fd7d2a45fd0dc06877098cda21669e196186804c2ec9b1183d'
SUBSETS = tuple(s for size in (1, 2, 3) for s in combinations(range(3), size))
MODEL = dict(actual_anchored_physical_space_time_compactness=True,
             anchored_reanchored_temporal_semigroup_constructed=True,
             nontrivial_temporal_state=False,
             actual_uniform_physical_smearing_Dirichlet_budget=True,
             nonzero_anchored_limit=False,
             full_time_Euclidean_covariance=False,
             continuum_Yang_Mills=False, uniform_positive_gap=False,
             local_field_identification=False,
             heat_filter_survival_lower_bound_proved=False,
             energy_upper_bound_alone_implies_survival=False,
             ordinary_increment_survival_implies_filter_survival=False,
             actual_all_coupling_center_sector_gap_upper_bound=True,
             actual_explicit_weak_coupling_gap_closing_trajectory=True,
             center_sector_bound_applies_to_center_neutral_gap=False,
             winding_states_are_local_continuum_fields=False,
             weak_coupling_alone_selects_Yang_Mills_continuum=False,
             finite_controls_are_universal_analytic_proof=False,
             finite_spectral_fixtures_are_Yang_Mills_vacua=False)


def trace_channels(quaternions):
    qs = tuple(previous.quaternion(q) for q in quaternions)
    if len(qs) != 3:
        raise ValueError('Three exact SU(2) holonomies are required')
    traces = {}
    for subset in SUBSETS:
        out = model.y.ONE
        for i in subset:
            out = model.y.qmul(out, qs[i])
        traces[subset] = out[0]
    A = sum((traces[(i,)]**2 for i in range(3)), F(0))
    B = sum((traces[s]**2 for s in SUBSETS if len(s) == 2), F(0))
    C = traces[(0, 1, 2)]**2
    vectors = [q[1:] for q in qs]
    gram = [[sum(x*y for x, y in zip(left, right)) for right in vectors] for left in vectors]
    determinant = (gram[0][0]*(gram[1][1]*gram[2][2]-gram[1][2]*gram[2][1])
                   -gram[0][1]*(gram[1][0]*gram[2][2]-gram[1][2]*gram[2][0])
                   +gram[0][2]*(gram[1][0]*gram[2][1]-gram[1][1]*gram[2][0]))
    triple = (sum(qs[i][0]*traces[tuple(j for j in range(3) if j != i)] for i in range(3))
              -2*qs[0][0]*qs[1][0]*qs[2][0]-traces[(0, 1, 2)])
    assert determinant == triple*triple
    assert A+B+C > F(1, 5)
    return dict(singleton_squares=A, pair_squares=B, triple_square=C,
                total_trace_squares=A+B+C, imaginary_Gram=gram,
                imaginary_Gram_determinant=determinant, signed_triple_product=triple,
                traces=[dict(subset=list(s), value=traces[s]) for s in SUBSETS])


def word_gradient(quaternions, word):
    """Prefix/suffix replay for an edge-simple closed walk, including inverses."""
    qs = tuple(previous.quaternion(q) for q in quaternions)
    model.incidence((word,), len(qs))
    if not word:
        raise ValueError('A nonempty word with distinct edges is required')
    factors = [qs[e] if sign == 1 else model.y.qdagger(qs[e]) for e, sign in word]
    prefix, suffix = [model.y.ONE], [model.y.ONE]*(len(word)+1)
    for factor in factors:
        prefix.append(model.y.qmul(prefix[-1], factor))
    for i in range(len(word)-1, -1, -1):
        suffix[i] = model.y.qmul(factors[i], suffix[i+1])
    gradients = {}
    for i, (edge, sign) in enumerate(word):
        values = []
        for axis in (1, 2, 3):
            raw = (model.y.qmul(model.UNITS[axis], qs[edge]) if sign == 1
                   else model.y.qmul(model.y.qdagger(qs[edge]), model.UNITS[axis]))
            derivative = tuple(-sign*x/2 for x in raw)
            values.append(model.y.qmul(model.y.qmul(prefix[i], derivative), suffix[i+1])[0])
        gradients[edge] = tuple(values)
    return prefix[-1][0], gradients


@lru_cache(maxsize=None)
def center_geometry(side):
    graph, qs, _, _ = previous.geometry_data(side)
    indices = {edge: i for i, edge in enumerate(graph['edges'])}
    loops = []
    for axis in range(3):
        words = []
        for k in range(side):
            vertex = tuple(k if j == axis else 0 for j in range(3))
            words.append((indices[vertex, axis], 1))
        loops.append(tuple(words))
    assert all(set(e for e, _ in loops[i]).isdisjoint(e for e, _ in loops[j])
               for i, j in combinations(range(3), 2))
    seams = [{i for i, (vertex, direction) in enumerate(graph['edges'])
              if direction == axis and vertex[axis] == side-1} for axis in range(3)]
    for seam in seams:
        assert all(sum(e in seam for e, _ in face) % 2 == 0 for face in graph['faces'])
    holonomies = [model.loop_value(qs, loop) for loop in loops]
    trace_control = trace_channels(holonomies)
    rows = []
    for subset in SUBSETS:
        word = tuple(edge for i in subset for edge in loops[i])
        value, gradients = word_gradient(qs, word)
        energy = sum(map(norm2, gradients.values()), F(0))
        assert energy == F(len(word), 4)*(1-value*value)
        parities = tuple(sum(e in seam for e, _ in word) % 2 for seam in seams)
        assert parities == tuple(int(i in subset) for i in range(3))
        for axis, seam in enumerate(seams):
            transformed = tuple(tuple((-1 if e in seam else 1)*x for x in q) for e, q in enumerate(qs))
            assert model.loop_value(transformed, word)[0] == (-1)**parities[axis]*value
        rows.append(dict(subset=list(subset), length=len(word), center_parities=parities,
                         trace=value, exact_Dirichlet_density=energy))
    characters = [[(-1)**sum(bits[i] for i in subset) for bits in product((0, 1), repeat=3)]
                  for subset in SUBSETS]
    assert a.multiply(characters, a.transpose(characters)) == a.scale(a.identity(7), 8)
    return dict(side=side, rooted_axial_loop_lengths=[len(loop) for loop in loops],
                center_seam_sizes=list(map(len, seams)), all_plaquettes_have_even_seam_parity=True,
                edge_disjoint_rooted_cycles=True, channel_controls=rows,
                seven_character_Gram=a.scale(a.identity(7), 8), trace_algebra=trace_control,
                sum_channel_lengths=12*side, finite_replay_not_universal_proof=True)


def center_state_budget(side, clock=1, threshold=1, time=1):
    native.width(side, 4)
    clock, threshold, time = positive(clock), positive(threshold), nonnegative(time)
    return dict(side=side, physical_clock=clock, actual_state='sum of seven nonempty winding subset traces times h',
                state_norm_squared=dict(lower=F(1, 5), upper=F(7)),
                bare_form_energy_upper=F(3*side), physical_form_energy_upper=3*clock*side,
                full_bare_gap_upper=F(15*side), full_physical_gap_upper=15*clock*side,
                threshold=threshold, normalized_mass_above_threshold_upper=min(F(1), 15*clock*side/threshold),
                time=time, unnormalized_heat_loss_upper=min(F(7), 3*time*clock*side),
                includes_nontrivial_global_center_sectors=True,
                bound_for_center_neutral_or_local_QFT_gap=False)


def closing_trajectory(index):
    native.width(index, 2)
    n = index
    spacing, side, coupling_squared = F(1, n), n*n, F(1, n**4)
    theta, clock = F(n**8), F(1, n**3)
    assert theta == 1/coupling_squared**2 and clock == coupling_squared/spacing
    return dict(index=n, spacing=spacing, side=side, physical_extent=side*spacing,
                bare_coupling_squared=coupling_squared, theta=theta, physical_clock=clock,
                actual_center_state=center_state_budget(side, clock),
                actual_plaquette_deficiency_upper=F(9, 4*n**4),
                actual_local_UV_band=dict(lower=F(n, 45), upper=F(80*n), mass_lower=F(1, 3600)),
                selected_asymptotically_free_trajectory=False,
                local_continuum_fields_or_survival_asserted=False,
                large_lattice_geometry_allocated=False)


def polynomial(value):
    out = {}
    for monomial, coefficient in value.items():
        if len(monomial) != 3:
            raise ValueError('A three-variable polynomial is required')
        for power in monomial:
            native.width(power, 0)
        coefficient = exact(coefficient)
        if coefficient:
            out[tuple(monomial)] = coefficient
    return out


def polynomial_square(poly):
    out = {}
    for left, x in poly.items():
        for right, y in poly.items():
            powers = tuple(i+j for i, j in zip(left, right))
            out[powers] = out.get(powers, F(0))+x*y
    return out


def integrate_box(poly, lower, upper):
    lower, upper = tuple(map(exact, lower)), tuple(map(exact, upper))
    if len(lower) != 3 or len(upper) != 3 or any(lo > hi for lo, hi in zip(lower, upper)):
        raise ValueError('An ordered three-dimensional exact box is required')
    total = F(0)
    for powers, coefficient in poly.items():
        term = coefficient
        for power, lo, hi in zip(powers, lower, upper):
            term *= (hi**(power+1)-lo**(power+1))/(power+1)
        total += term
    return total


def cell_jensen(poly, cells_per_axis):
    """Exact averages on [0,1]^3; a cubical fixture for the written L2 theorem."""
    poly = polynomial(poly)
    native.width(cells_per_axis)
    n, averages = cells_per_axis, []
    spacing = F(1, n)
    for index in product(range(n), repeat=3):
        lower = tuple(F(i, n) for i in index)
        upper = tuple(F(i+1, n) for i in index)
        averages.append(integrate_box(poly, lower, upper)/spacing**3)
    discrete = spacing**3*sum((value*value for value in averages), F(0))
    continuous = integrate_box(polynomial_square(poly), (F(0),)*3, (F(1),)*3)
    assert discrete <= continuous
    return dict(cells_per_axis=n, spacing=spacing, cell_count=n**3,
                cell_average_square_sum=sum((value*value for value in averages), F(0)),
                volume_weighted_square_sum=discrete, continuum_L2_squared=continuous,
                exact_Jensen_defect=continuous-discrete)


def physical_smear_budget(theta, spacing, average_square_sum, continuum_L2_squared,
                          sqrt_electric_magnetic=1):
    theta, spacing, scale = positive(theta), positive(spacing), positive(sqrt_electric_magnetic)
    average_square_sum, continuum_L2_squared = nonnegative(average_square_sum), nonnegative(continuum_L2_squared)
    discrete = spacing**3*average_square_sum
    if discrete > continuum_L2_squared:
        raise ValueError('The supplied cell averages violate the required Jensen bound')
    root = previous.stage74.sqrt_interval(theta)
    Cf = 4*scale*continuum_L2_squared
    return dict(theta=theta, spacing=spacing, sqrt_cE_cB=scale,
                smear_coefficient_fourth_power=theta*spacing**8,
                smear_coefficient_squared=dict(lower=root['lower']*spacing**4,
                                                upper=root['upper']*spacing**4),
                exact_physical_discrete_Dirichlet_upper=4*scale*discrete,
                uniform_physical_carre_du_champ_upper=Cf,
                physical_form_energy_upper=Cf,
                nonzero_variance_lower_bound_provided=False,
                actual_angle_incidence_bound_used=True)


def increment_moment_budget(Cf, time_difference, order):
    Cf, dt = nonnegative(Cf), nonnegative(time_difference)
    native.width(order)
    return dict(Cf=Cf, time_difference=dt, half_moment_order=order,
                even_moment_upper=2*factorial(order)*(4*Cf*dt)**order,
                second_moment_upper=2*Cf*dt,
                requires_uniform_pointwise_carre_du_champ_not_just_mean_energy=True,
                anchored_time_increment=True)


def filter_budget(time, Cf):
    time, Cf = positive(time), positive(Cf)
    return dict(filter_time=time, source_energy_upper=Cf, filtered_norm_squared_upper=time*Cf,
                filtered_energy_upper=Cf, filtered_inverse_energy_upper=time*time*Cf,
                uniform_positive_filtered_norm_inferred=False)


def filter_survival_budget(time, Cf, norm_squared_floor):
    row = filter_budget(time, Cf)
    floor = positive(norm_squared_floor)
    if floor > row['filtered_norm_squared_upper']:
        raise ValueError('The proposed survival floor exceeds the proved norm ceiling')
    lower, upper = floor/(4*time*time*Cf), 4*Cf/floor
    assert lower < upper
    return dict(**row, independently_required_norm_squared_floor=floor,
                physical_band=dict(lower=lower, upper=upper),
                unnormalized_band_mass_lower=floor/2,
                # Discarded mass <= floor/2 <= actual filtered mass/2.
                normalized_band_mass_lower=F(1, 2),
                normalized_physical_mean_energy=dict(lower=floor/(time*time*Cf), upper=Cf/floor),
                low_band_mass_upper=floor/4, high_band_mass_upper=floor/4,
                survival_premise_is_not_discharged=True)


def finite_filter_control(atoms, weights, time=1, Cf=None):
    atoms, weights = tuple(map(nonnegative, atoms)), tuple(map(nonnegative, weights))
    time = positive(time)
    if not atoms or len(atoms) != len(weights):
        raise ValueError('A finite positive spectral measure with matching dimensions is required')
    source_energy = sum((x*w for x, w in zip(atoms, weights)), F(0))
    Cf = positive(source_energy if Cf is None else Cf)
    if source_energy > Cf:
        raise ValueError('The spectral source violates the declared energy budget')
    filtered, increments = [], []
    for energy, weight in zip(atoms, weights):
        exp = previous.exp_negative_interval(time*energy)
        lo = (1-exp['upper'])**2
        hi = min(F(1), (time*energy)**2, (1-exp['lower'])**2)
        filtered.append(dict(lower=weight*lo, upper=weight*hi))
        increments.append(dict(lower=2*weight*(1-exp['upper']),
                               upper=2*weight*min(F(1), time*energy, 1-exp['lower'])))
    mass = {key: sum((row[key] for row in filtered), F(0)) for key in ('lower', 'upper')}
    energy = {key: sum((x*row[key] for x, row in zip(atoms, filtered)), F(0)) for key in ('lower', 'upper')}
    inverse_energy = {key: sum((row[key]/x for x, row in zip(atoms, filtered) if x), F(0))
                      for key in ('lower', 'upper')}
    increment = {key: sum((row[key] for row in increments), F(0)) for key in ('lower', 'upper')}
    assert 0 <= mass['lower'] <= mass['upper'] <= time*Cf
    assert 0 <= energy['lower'] <= energy['upper'] <= Cf
    assert 0 <= inverse_energy['lower'] <= inverse_energy['upper'] <= time*time*Cf
    assert 0 <= increment['lower'] <= increment['upper'] <= 2*time*Cf
    return dict(**filter_budget(time, Cf), atoms=atoms, weights=weights,
                source_norm_squared=sum(weights), source_energy=source_energy,
                filtered_atom_masses=filtered, filtered_norm_squared=mass,
                filtered_energy=energy, filtered_inverse_energy=inverse_energy,
                ordinary_increment_variance=increment,
                algebraic_fixture_not_Yang_Mills_vacuum=True)


def finite_survival_control(atoms, weights, time, Cf, norm_squared_floor):
    time, Cf, norm_squared_floor = positive(time), positive(Cf), positive(norm_squared_floor)
    control = finite_filter_control(atoms, weights, time, Cf)
    budget = filter_survival_budget(time, Cf, norm_squared_floor)
    if control['filtered_norm_squared']['lower'] < norm_squared_floor:
        raise ValueError('The exact directed control does not certify the proposed survival floor')
    band = budget['physical_band']
    masses = control['filtered_atom_masses']
    direct = sum((row['lower'] for atom, row in zip(control['atoms'], masses)
                  if band['lower'] <= atom <= band['upper']), F(0))
    low = sum((row['upper'] for atom, row in zip(control['atoms'], masses) if atom < band['lower']), F(0))
    high = sum((row['upper'] for atom, row in zip(control['atoms'], masses) if atom > band['upper']), F(0))
    assert low <= budget['low_band_mass_upper'] and high <= budget['high_band_mass_upper']
    return dict(control=control, conditional_budget=budget,
                independently_certified_band_mass_lower=max(direct, norm_squared_floor/2),
                low_mass_upper=low, high_mass_upper=high)


@lru_cache(maxsize=None)
def all_controls():
    palette = (model.y.ONE, model.UNITS[1], model.UNITS[2], model.UNITS[3],
               (F(3, 5), F(4, 5), 0, 0), (0, F(3, 5), F(4, 5), 0))
    traces = [trace_channels(qs)['total_trace_squares'] for qs in product(palette, repeat=3)]
    assert F(32, 135) > F(1, 5)
    polynomials = ({(0, 0, 0): 1}, {(1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): -1},
                   {(2, 0, 0): 3, (0, 1, 1): -2, (0, 0, 0): F(1, 3)})
    cells = [cell_jensen(poly, n) for poly in polynomials for n in (2, 3, 4)]
    smearing = [physical_smear_budget(26, row['spacing'], row['cell_average_square_sum'],
                                     row['continuum_L2_squared']) for row in cells]
    ultraviolet = [finite_filter_control((n,), (F(1, n),)) for n in (2, 4, 32, 64, 1024)]
    infrared = [finite_filter_control((F(1, n),), (n,)) for n in (2, 4, 32, 64, 1024)]
    for n, row in zip((2, 4, 32, 64, 1024), ultraviolet):
        assert row['filtered_norm_squared']['upper'] <= F(1, n)
    for n, row in zip((2, 4, 32, 64, 1024), infrared):
        assert row['filtered_norm_squared']['upper'] <= F(1, n)
        assert row['ordinary_increment_variance']['lower'] >= 2-F(1, n)
    return dict(actual_center_geometry=[center_geometry(n) for n in (4, 5)],
                quaternion_trace_controls=dict(exact_triples_checked=len(traces),
                    minimum_replay_trace_sum=min(traces), universal_written_lower_bound=F(1, 5),
                    contradiction_determinant_lower=F(32, 135), contradiction_determinant_upper=F(1, 5)),
                actual_gap_closing_trajectories=[closing_trajectory(n) for n in (2, 4, 32, 64, 256)],
                exact_cell_Jensen_controls=cells, actual_physical_smearing=smearing,
                anchored_increment_moments=[increment_moment_budget(4, F(1, n), order)
                    for n in (1, 4, 16) for order in (1, 2, 3)],
                conditional_survival=filter_survival_budget(1, 1, F(1, 4)),
                certified_finite_survival=finite_survival_control((F(1, 100), 1, 100),
                    (1, F(1, 2), F(1, 1000)), 1, 1, F(1, 8)),
                ultraviolet_filter_collapse=ultraviolet, infrared_filter_collapse=infrared)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 274 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM77 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM77 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM76 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM77_PHYSICAL_SMEARING_CONTINUUM_SURVIVAL',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_PHYSICAL_COMPACTNESS_AND_SURVIVAL_CRITERION__NONZERO_YM_LIMIT_OPEN',
        proof_scope='Written actual-model estimates and explicit center states with exact finite replay controls',
        open='Nonzero local anchored limit, surviving filtered states, Euclidean covariance and positive local gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM77_PHYSICAL_SMEARING_CONTINUUM_SURVIVAL':
        raise ValueError('YM77 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM77 fresh certificate/pin mismatch')
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
    print('YM77 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], actual_center_gap_upper='15*kappa*L',
                          anchored_compactness=True, nonzero_local_continuum_limit=False)))


if __name__ == '__main__':
    main()
