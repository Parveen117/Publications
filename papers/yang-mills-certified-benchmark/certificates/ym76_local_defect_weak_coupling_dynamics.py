"""YM76: actual local defects, weak-coupling vacuum scales and microscopic dynamics.

Python 3.12 only. Exact finite geometry replays accompany the written universal
proofs; beta and spectral fixtures are not interacting Yang--Mills vacua.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym75_rh_source_resolved_transfer as previous

stage74 = previous.previous
native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative, encode = previous.positive, previous.nonnegative, previous.encode
model, GRID = stage74.model, stage74.GRID
spatial = stage74.previous.previous.spatial
RESULT = HERE/'YM76_RESULT.json'
PIN = HERE/'EXPECTED_YM76.sha256'
SOURCES = HERE/'YM76_SOURCE_MANIFEST.json'
BASE_TREE = 'f781c1781874a45dded79cffa8dbad4fc05696dd'
PARENT_HASH = '5a4afd2e02ec5eaa45d98c3936fc4e97170d82b8ed4c61746ad3ed313c4ebc37'
MODEL = dict(actual_SU2_local_defect_factorization=True,
             volume_uniform_two_sided_weak_coupling_vacuum_scales=True,
             actual_centered_plaquette_angle_state_nonzero=True,
             actual_microscopic_spectral_band_weight=True,
             actual_collective_radial_Dirichlet_Gram_upper_bound=True,
             finite_geometry_replay_is_universal_proof=False,
             radial_beta_fixtures_are_interacting_vacua=False,
             microscopic_clock_is_physical_clock=False,
             physical_ultraviolet_weight_escape=True,
             physical_limit_must_reproduce_lattice_heat_correlations=True,
             unfiltered_probe_has_norm_preserving_strongly_continuous_physical_limit=False,
             zero_atom_allowed_in_limiting_source_spectrum=True,
             local_probe_excludes_other_slow_modes=False,
             uniform_weak_coupling_positive_gap=False,
             continuum_R4_QFT=False, physical_mass_prediction=False)


def norm2(vector):
    return sum((value*value for value in vector), F(0))


def quaternion(value):
    out = tuple(map(exact, value))
    if len(out) != 4 or norm2(out) != 1:
        raise ValueError('An exact unit quaternion is required')
    return out


def face_gradient(quaternions, face):
    """Native left generators D_a q=-(e_a q)/2, including inverse links."""
    qs = tuple(quaternion(q) for q in quaternions)
    model.incidence((face,), len(qs))
    if len(face) != 4:
        raise ValueError('An elementary four-link face is required')
    factors = [qs[e] if sign == 1 else model.y.qdagger(qs[e]) for e, sign in face]
    prefix, suffix = [model.y.ONE], [model.y.ONE]*(len(face)+1)
    for factor in factors:
        prefix.append(model.y.qmul(prefix[-1], factor))
    for i in range(len(face)-1, -1, -1):
        suffix[i] = model.y.qmul(factors[i], suffix[i+1])
    gradient = {}
    for i, (edge, sign) in enumerate(face):
        values = []
        for axis in (1, 2, 3):
            raw = (model.y.qmul(model.UNITS[axis], qs[edge]) if sign == 1
                   else model.y.qmul(model.y.qdagger(qs[edge]), model.UNITS[axis]))
            derivative = tuple(-sign*value/2 for value in raw)
            values.append(model.y.qmul(model.y.qmul(prefix[i], derivative), suffix[i+1])[0])
        gradient[edge] = tuple(values)
    return prefix[-1][0], gradient


@lru_cache(maxsize=None)
def geometry_data(side, kind='pattern'):
    graph = model.lattice(side)
    palette = (model.y.ONE, (F(3, 5), F(4, 5), F(0), F(0)),
               (F(5, 13), F(0), F(12, 13), F(0)),
               (F(8, 17), F(0), F(0), F(15, 17)), model.UNITS[1])
    if kind == 'identity':
        qs = tuple(model.y.ONE for _ in graph['edges'])
    elif kind == 'center':
        qs = tuple(tuple((-1 if i % 7 == 0 else 1)*x for x in model.y.ONE)
                   for i in range(len(graph['edges'])))
    elif kind == 'pattern':
        qs = tuple(palette[(i+sum(vertex)+2*axis) % len(palette)]
                   for i, (vertex, axis) in enumerate(graph['edges']))
    else:
        raise ValueError('Unknown exact geometry fixture')
    # Validate once here; each local gradient only needs its four quaternions.
    qs = tuple(quaternion(q) for q in qs)
    values, gradients = [], []
    for face in graph['faces']:
        local = tuple(qs[e] for e, _ in face)
        local_face = tuple((i, sign) for i, (_, sign) in enumerate(face))
        value, local_gradient = face_gradient(local, local_face)
        values.append(value)
        gradients.append({edge: local_gradient[i] for i, (edge, _) in enumerate(face)})
    return graph, qs, tuple(values), tuple(gradients)


def geometry_replay(side, kind='pattern'):
    """Finite exact replay of R=8X-Gamma(V), not a proof by sampling."""
    graph, _, values, gradients = geometry_data(side, kind)
    count = len(values)
    incidence, _ = model.incidence(graph['faces'], len(graph['edges']))
    assert set(incidence) == {4} and count == len(graph['edges']) == 3*side**3
    at_edge = [[] for _ in graph['edges']]
    for value, gradient in zip(values, gradients):
        assert sum(map(norm2, gradient.values()), F(0)) == 1-value*value
        for edge, vector in gradient.items():
            assert norm2(vector) == (1-value*value)/4
            at_edge[edge].append(vector)
    Gamma = sum((norm2(tuple(sum(v[k] for v in vectors) for k in range(3)))
                 for vectors in at_edge), F(0))
    mismatch = sum((norm2(tuple(x-y for x, y in zip(left, right)))
                    for vectors in at_edge for left, right in combinations(vectors, 2)), F(0))
    X = sum((1-value for value in values), F(0))
    on_site = 4*sum(((1-value)**2 for value in values), F(0))
    residual = 8*X-Gamma
    assert residual == on_site+mismatch >= 0
    square_replays = []
    V = count-X
    for t in (F(1, 16), F(1, 2), F(2)):
        theta = 3*t+8*t*t
        Q_potential = -3*t*V+t*t*Gamma
        assert theta*X == Q_potential+3*t*count+t*t*residual
        square_replays.append(dict(t=t, theta=theta, scalar_identity_error=F(0)))
    return dict(side=side, configuration=kind, links=count, plaquettes=count,
                plaquettes_per_link=4, compared_incident_pairs=6*count,
                X=X, Gamma_V=Gamma, local_square_deficiency=on_site,
                local_gradient_mismatch=mismatch, exact_nonnegative_defect=residual,
                square_completion_replays=square_replays,
                finite_replay_not_universal_proof=True)


def collective_replay(side, kind='pattern'):
    """Pointwise radial coefficients c_p=d_p sqrt(1-W_p^2), with rational d_p.

    This replays the pointwise Gram inequality; it computes no interacting
    vacuum covariance or expectations. No irrational gradients are required.
    """
    graph, _, values, gradients = geometry_data(side, kind)
    rows = []
    patterns = (tuple(F(1) for _ in values),
                tuple(F((i % 7)-3, (i % 3)+1) for i in range(len(values))))
    for coefficients in patterns:
        at_edge = [[] for _ in graph['edges']]
        for coefficient, gradient in zip(coefficients, gradients):
            for edge, vector in gradient.items():
                at_edge[edge].append(tuple(coefficient*x for x in vector))
        energy = sum((norm2(tuple(sum(v[k] for v in vectors) for k in range(3)))
                      for vectors in at_edge), F(0))
        norm = sum((c*c*(1-w*w) for c, w in zip(coefficients, values)), F(0))
        remainder = sum((norm2(tuple(x-y for x, y in zip(left, right)))
                         for vectors in at_edge for left, right in combinations(vectors, 2)), F(0))
        assert 4*norm-energy == remainder >= 0
        rows.append(dict(radial_coefficient_norm_squared=norm,
                         collective_Dirichlet_energy=energy,
                         four_I_minus_K_quadratic_form=remainder))
    # Disjoint support is the geometric reason for exact orthogonality.
    chosen, used = [], set()
    for i, gradient in enumerate(gradients):
        if not used.intersection(gradient):
            chosen.append(i)
            used.update(gradient)
    assert all(set(gradients[i]).isdisjoint(gradients[j]) for i, j in combinations(chosen, 2))
    return dict(side=side, configuration=kind, incidence_operator_bound=4,
                exact_congruence_controls=rows, edge_disjoint_plaquette_count=len(chosen),
                disjoint_radial_Dirichlet_Gram_is_identity=True,
                independence_of_vacuum_plaquettes_assumed=False)


def weak_budget(theta):
    theta = positive(theta)
    if theta < 25:
        raise ValueError('The simple weak-coupling constants require theta >= 25')
    root = stage74.sqrt_interval(theta)
    lo, hi = root['lower'], root['upper']
    vacuum = stage74.vacuum_budget(theta)
    density = vacuum['actual_shifted_vacuum_deficit_upper']
    t_root = stage74.sqrt_interval(9+32*theta)
    t = dict(lower=(t_root['lower']-3)/16, upper=(t_root['upper']-3)/16)
    pi = spatial.pi_interval()
    pi_lo, pi_hi = pi.lo, pi.hi
    assert pi_hi*pi_hi < 10
    assert density <= F(9, 4)*lo and density/theta <= F(1, 2)
    assert 3*t['lower'] <= density
    return dict(theta=theta, sqrt_theta=root, square_completion_parameter=t,
                vacuum_excess_density=dict(lower=3*t['lower'], upper=density),
                actual_link_electric_energy=dict(lower=lo/32, upper=F(9, 4)*hi),
                actual_plaquette_deficiency=dict(lower=1/(32*hi), upper=F(9, 4)/lo),
                actual_trace_mean_lower=1-density/theta,
                actual_trace_gradient_mean_lower=1/(16*hi),
                actual_local_defect_density_upper=(density-3*t['lower'])/t['lower']**2,
                actual_trace_variance_lower=1/(2*pi_hi**4*(1+4*density)**2),
                actual_angle_variance=dict(lower=1/(40*hi), upper=9*pi_hi*pi_hi/(8*lo)),
                sharper_angle_variance_lower=1/(16*density+4),
                sharper_angle_variance_upper=pi_hi*pi_hi*density/(2*theta),
                pi_interval=dict(lower=pi_lo, upper=pi_hi),
                microscopic_first_moment=dict(lower=F(4, 45), upper=F(40)),
                microscopic_second_over_first_upper=F(40),
                explicit_rescaled_state_norm_squared=dict(lower=F(1, 40), upper=9*pi_hi*pi_hi/8),
                explicit_rescaled_state_first_moment=F(1),
                explicit_rescaled_state_second_moment_upper=F(40),
                bounds_are_actual_vacuum_and_uniform_in_volume=True)


def radial_beta(index):
    """Haar f_n=(1+W)^n fixture; evaluates singular terms through beta moments."""
    native.width(index, 0)
    alpha, beta = F(4*index+3, 2), F(3, 2)
    inverse_y = (alpha+beta-1)/(alpha-1)
    inverse_one_minus_y = (alpha+beta-1)/(beta-1)
    cot_squared = (inverse_y+inverse_one_minus_y)/4-1
    mixed = 1-inverse_y/2
    radial_energy = index*index*(inverse_y-1)
    J_squared = radial_energy+cot_squared-2*index*mixed
    assert J_squared == radial_energy+1
    return dict(index=index, beta_alpha=alpha, beta_beta=beta,
                radial_electric_energy=radial_energy, cotangent_squared_mean=cot_squared,
                cosine_over_one_plus_cosine_mean=mixed, J_norm_squared=J_squared,
                exact_radial_cancellation=J_squared-radial_energy,
                interacting_vacuum=False)


def spectral_band_budget(first_moment_floor=F(4, 45), second_over_first=40):
    lower, ratio = positive(first_moment_floor), positive(second_over_first)
    if lower > ratio:
        raise ValueError('The moments of a probability measure require a <= b')
    return dict(first_moment_floor=lower, second_over_first_upper=ratio,
                band_lower=lower/4, band_upper=2*ratio,
                band_mass_lower=lower/(8*ratio),
                zero_atom_excluded=False, full_gap_inferred=False)


def finite_spectral_control(atoms, weights, first_moment_floor=F(4, 45), second_over_first=40):
    atoms, weights = tuple(map(nonnegative, atoms)), tuple(map(nonnegative, weights))
    if not atoms or len(atoms) != len(weights) or sum(weights) != 1:
        raise ValueError('A normalized finite nonnegative spectral measure is required')
    row = spectral_band_budget(first_moment_floor, second_over_first)
    first = sum((x*w for x, w in zip(atoms, weights)), F(0))
    second = sum((x*x*w for x, w in zip(atoms, weights)), F(0))
    if first < row['first_moment_floor'] or second > row['second_over_first_upper']*first:
        raise ValueError('The source moments do not satisfy the claimed premises')
    mass = sum((w for x, w in zip(atoms, weights) if row['band_lower'] <= x <= row['band_upper']), F(0))
    assert mass >= row['band_mass_lower']
    return dict(**row, atoms=atoms, weights=weights, first_moment=first,
                second_moment=second, actual_band_mass=mass,
                algebraic_fixture_not_gauge_vacuum=True)


def exp_negative_interval(value):
    """Directed rational Taylor/repeated-squaring enclosure of exp(-value)."""
    value = nonnegative(value)
    if not value:
        return dict(lower=F(1), upper=F(1))
    scaled, squarings = value, 0
    while scaled > 1:
        scaled /= 2
        squarings += 1
    total, term = F(1), F(1)
    for k in range(1, 49):
        term *= scaled/k
        total += term
    omitted = term*scaled/49
    upper_exp = total+omitted/(1-scaled/50)
    lower, upper = 1/upper_exp, 1/total
    for _ in range(squarings):
        lower = F((lower*lower*GRID).__floor__(), GRID)
        upper = F((upper*upper*GRID).__ceil__(), GRID)
    return dict(lower=lower, upper=min(F(1), upper))


def heat_budget(time):
    time = nonnegative(time)
    band = spectral_band_budget()
    lower = exp_negative_interval(40*time)['lower']
    upper = 1-band['band_mass_lower']*(1-exp_negative_interval(time/45)['upper'])
    assert 0 <= lower <= upper <= 1
    return dict(microscopic_time=time, heat_correlation_lower=lower,
                heat_correlation_upper=upper, continuous_at_zero=True,
                decay_to_zero_not_proved=True)


def physical_clock_budget(time, spacing, sqrt_electric_magnetic=1):
    time, spacing, scale = nonnegative(time), positive(spacing), positive(sqrt_electric_magnetic)
    # H_phys=(sqrt(c_E c_B)/a) (H/sqrt(theta)).
    row = heat_budget(time*scale/spacing)
    return dict(physical_time=time, spacing=spacing, sqrt_cE_cB=scale,
                microscopic_time=row['microscopic_time'],
                physical_heat_correlation_upper=row['heat_correlation_upper'],
                positive_time_lim_sup_upper=F(3599, 3600) if time else F(1),
                physical_energy_band=dict(lower=scale/(45*spacing), upper=80*scale/spacing),
                escaping_spectral_mass_lower=F(1, 3600),
                excludes_norm_preserving_strongly_continuous_unfiltered_probe_limit=True,
                limit_required_to_reproduce_lattice_heat_correlations=True,
                excludes_smeared_fields_or_remaining_infrared_component=False)


@lru_cache(maxsize=None)
def all_controls():
    geometry = [geometry_replay(side, kind) for side in (4, 5)
                for kind in ('identity', 'center', 'pattern')]
    weak = [weak_budget(theta) for theta in (25, 26, 256, 4096, 65536, 10**6)]
    radial = [radial_beta(n) for n in (0, 1, 2, 8, 64)]
    collapsing = [finite_spectral_control((F(1, n), 1), (F(1, 2), F(1, 2)))
                  for n in (2, 16, 256, 65536)]
    zero = finite_spectral_control((0, 1), (F(1, 2), F(1, 2)))
    # A first moment alone permits complete escape from every fixed finite band.
    no_second = [dict(parameter=n, atoms=[0, n], weights=[1-F(1, n), F(1, n)],
                      first_moment=F(1), second_moment=F(n)) for n in (128, 1024, 8192)]
    pi = spatial.pi_interval()
    free_variance = dict(lower=pi.lo*pi.lo/12-F(1, 2), upper=pi.hi*pi.hi/12-F(1, 2))
    assert free_variance['lower'] > F(1, 4)
    return dict(actual_local_geometry=geometry,
                collective_incidence=[collective_replay(side) for side in (4, 5)],
                weak_coupling_budgets=weak, radial_beta_identities=radial,
                free_Haar_angle=dict(variance=free_variance, Dirichlet_energy=F(1), second_energy_moment=F(4)),
                actual_microscopic_band=spectral_band_budget(),
                microscopic_heat=[heat_budget(s) for s in (0, F(1, 100), F(1, 4), 1, 10)],
                physical_probe_budgets=[physical_clock_budget(1, F(1, n)) for n in (1, 16, 256, 4096)],
                closing_gap_source_fixtures=collapsing, limiting_zero_atom_fixture=zero,
                first_moment_without_second_moment_failure=no_second)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 267 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM76 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM76 base tree mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM75 parent result/pin mismatch')
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
    return encode(dict(certificate_type='YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        parent_result_sha256=PARENT_HASH, controls=all_controls(), model=MODEL,
        claim_status='ACTUAL_WEAK_COUPLING_LOCAL_DYNAMICS_PROVED__CONTINUUM_GAP_OPEN',
        proof_scope='Written actual-model proofs with exact finite geometry, beta and spectral controls',
        open='Uniform positive physical infrared gap and nontrivial local four-dimensional field construction'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM76_LOCAL_DEFECT_WEAK_COUPLING_DYNAMICS':
        raise ValueError('YM76 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM76 fresh certificate/pin mismatch')
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
    print('YM76 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], microscopic_band=['1/45', '80'],
                          microscopic_band_mass_lower='1/3600', continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
