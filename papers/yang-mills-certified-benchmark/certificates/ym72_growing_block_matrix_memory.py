"""YM72: growing-block matrix memory, moment budgets and visible relaxation.

Python 3.12 only. Directed budgets accompany written proofs for the actual
model. Finite fixtures check algebra and claim boundaries, not a discretized
vacuum. Default/--check are read-only; --write updates this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym71_interacting_memory_enclosures as previous

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
leading, spatial = previous.previous, previous.spatial
positive, nonnegative = leading.positive, leading.nonnegative
up, sqrt_upper, exponential_upper = previous.up, previous.sqrt_upper, previous.exponential_upper
GAP = F(1, 125)
MAX_RHO = F(3, 10**6)
RESULT = HERE/'YM72_RESULT.json'
PIN = HERE/'EXPECTED_YM72.sha256'
SOURCES = HERE/'YM72_SOURCE_PINS.json'


@lru_cache(maxsize=None, typed=True)
def envelope(rho):
    rho = nonnegative(rho)
    if rho > MAX_RHO:
        raise ValueError('The universal envelope requires rho <= 3/1000000')
    K = exponential_upper(4992*spatial.pi_interval().hi*rho)
    root_K = sqrt_upper(K)
    inverse_free_norm = sqrt_upper(F(24, 7))
    delta = F(512, 3)*rho
    forcing = delta+768*rho/GAP
    residual_over_rate = up(F(16, 9)*((240+F(8, 9)+8*sqrt_upper(80))*rho+7680*rho*rho))
    dn, fn, en = (up(value*inverse_free_norm) for value in (delta, forcing, residual_over_rate))
    norm = up(root_K+dn)
    lower_error = up(F(13, 9)*(K-1+2*root_K*dn+dn*dn+(2*fn+en)*norm))
    upper_error = up(lower_error+F(312, 7)/GAP*(forcing+residual_over_rate)**2)
    boundary_K = exponential_upper(384*spatial.pi_interval().hi*rho)
    return dict(rho=rho, gap_lower=GAP, union_density_ratio=K, boundary_density_ratio=boundary_K,
                lift_error=delta, forcing_error_per_theta=forcing,
                eigen_residual_over_minimum_rate=residual_over_rate,
                relative_lower_error=lower_error, relative_upper_error=upper_error,
                normalized_commutator_upper=up(2*upper_error*upper_error))


def parameters(theta, factor, readouts=1):
    theta = exact(theta)
    native.width(factor, 3)
    native.width(readouts)
    rho = abs(theta)*factor*readouts
    return dict(theta=theta, block_factor=factor, readouts=readouts, **envelope(rho))


def matrix_response_interval(theta, factor, readouts, damping, clock=F(1)):
    damping, clock = nonnegative(damping), positive(clock)
    cell = parameters(theta, factor, readouts)
    coefficient = leading.atom_response(leading.leading_atoms(factor), damping/clock)[0][0]
    center = clock*exact(theta)**2*coefficient
    return dict(lower=max(F(0), 1-cell['relative_lower_error'])*center,
                upper=(1+cell['relative_upper_error'])*center, leading=center,
                meaning='Scalar endpoints multiply the identity in readout coefficient space')


def scaled_response_interval(theta, factor, readouts, spectral):
    theta, spectral = exact(theta), nonnegative(spectral)
    if not theta:
        raise ValueError('Division by theta squared requires nonzero coupling')
    cell = parameters(theta, factor, readouts)
    center = leading.atom_response(leading.leading_atoms(factor), factor*spectral)[0][0]
    return dict(lower=max(F(0), 1-cell['relative_lower_error'])*center,
                upper=(1+cell['relative_upper_error'])*center, leading=center)


def limit_profile(spectral):
    spectral = nonnegative(spectral)
    return 1/(12*(spectral+3))


def finite_block_correction(spectral):
    spectral = nonnegative(spectral)
    return (spectral-1)/(24*(spectral+3)**2)


def shape_budget(theta, factor, readouts=1):
    cell = parameters(theta, factor, readouts)
    lower_free, upper_free = F(6*factor+3, 6*factor+8), 1+F(1, 2*factor)
    lower = max(F(0), 1-cell['relative_lower_error'])*lower_free
    upper = (1+cell['relative_upper_error'])*upper_free
    return dict(relative_lower=lower, relative_upper=upper,
                first_correction_error=up(cell['relative_upper_error']*upper_free+F(43, 36*factor*factor)))


def normalized_moment_interval(theta, factor, readouts, order):
    """Loewner interval for r^n M_n/theta^2, from positive response stencils."""
    theta = exact(theta)
    if not theta:
        raise ValueError('Nonzero coupling required for normalized moments')
    native.width(order, 0)
    cell = parameters(theta, factor, readouts)
    zero = scaled_response_interval(theta, factor, readouts, 0)
    gap = GAP/factor
    lower, upper = F(0), zero['upper']/gap**order
    if not order:
        lower = zero['lower']
    last_step = None
    if order:
        for exponent in range(-28, 17):
            step = gap*F(2)**exponent
            rows = [scaled_response_interval(theta, factor, readouts, k*step) for k in range(order+1)]
            candidate = previous.stencil_interval([(v['lower'], v['upper']) for v in rows], step, gap)
            lo, hi = max(lower, candidate['lower']), min(upper, candidate['upper'])
            if hi-lo < upper-lower:
                last_step = step
            lower, upper = lo, hi
    if lower > upper:
        raise ArithmeticError('Moment enclosures do not overlap')
    center = factor**order*leading.atom_moments(leading.leading_atoms(factor), order+1)[order][0][0]
    return dict(order=order, lower=lower, upper=upper, leading=center,
                limiting_moment=F(1, 12*3**(order+1)), last_improving_scaled_step=last_step,
                controlling_product=abs(theta)*readouts*factor**(order+1),
                certified_error_times_r_power=cell['relative_upper_error']*factor**order)


def dynamics_budget(theta, factor, readouts=1):
    cell = parameters(theta, factor, readouts)
    K = cell['boundary_density_ratio']
    residual = up(192*abs(exact(theta))*sqrt_upper(readouts*K))
    return dict(centered_Gram_lower=1/(4*K), centered_Gram_upper=K/4,
                normalized_generator_residual=residual,
                all_nonnegative_times_correlation_error=up(residual/3),
                spectral_weight_below_three_halves=up(F(4, 9)*residual*residual),
                full_physical_gap_upper_after_block_clock=F(100, 19*factor),
                full_physical_gap_lower_after_block_clock=GAP/factor)


def matrix_from_scalar_budget(theta, factor, readouts):
    """Direct multi-readout constants before the universal rho estimate."""
    base = previous.parameters(theta, factor)
    native.width(readouts)
    root_N = sqrt_upper(readouts)
    K = exponential_upper(4*spatial.pi_interval().hi*base['wave_score_upper']
                          *readouts*base['union_support_upper'])
    G2 = base['free_forcing_norm_squared']
    Gstar = sqrt_upper(K*G2)
    delta, force, residual = [up(root_N*base[key]) for key in
                            ('lift_error', 'forcing_error_per_theta', 'eigen_equation_residual')]
    H = up(Gstar+delta)
    lo_rate, hi_rate, gap = base['minimum_free_rate'], base['maximum_free_rate'], base['gap_lower']
    comparison = up((K-1)*G2+2*Gstar*delta+delta*delta)
    error = up(hi_rate/lo_rate*(comparison
               +2*force*H+H*residual/lo_rate)/G2)
    squared = up(hi_rate/gap*(force+residual/lo_rate)**2/G2)
    return dict(lower_error=error, upper_error=up(error+squared))


def star_geometry(factor, origin, side):
    """Local fine-link geometry, without allocating the ambient lattice."""
    native.width(factor, 3)
    native.width(side)
    if len(origin) != 3 or any(isinstance(v, bool) or not isinstance(v, int) for v in origin):
        raise ValueError('Three integer origin coordinates required')
    if any(v % factor for v in origin):
        raise ValueError('The loop origin must be a coarse vertex')
    if side % factor or side < 4*factor:
        raise ValueError('Admissible periodic fine side required')
    def vertex(v):
        return tuple(x % side for x in v)
    def move(v, axis, amount=1):
        w = list(v)
        w[axis] += amount
        return vertex(w)
    origin = vertex(origin)
    boundary = set()
    for i in range(factor):
        boundary.update(((move(origin, 0, i), 0),
                         (move(move(origin, 0, factor), 1, i), 1),
                         (move(move(origin, 1, factor), 0, i), 0),
                         (move(origin, 1, i), 1)))
    labels = set()
    for v, axis in boundary:
        for other in range(3):
            if other == axis:
                continue
            axes = tuple(sorted((axis, other)))
            labels.update(((v, axes), (move(v, other, -1), axes)))
    faces = []
    for v, (axis, other) in sorted(labels):
        faces.append(frozenset(((v, axis), (move(v, axis), other),
                                (move(v, other), axis), (v, other))))
    support = boundary.union(*faces)
    retained = lambda edge: all(edge[0][j] % factor == 0 for j in range(3) if j != edge[1])
    assert all(retained(edge) for edge in boundary)
    assert all(any(not retained(edge) for edge in face) for face in faces)
    return dict(boundary=frozenset(boundary), faces=tuple(faces), support=frozenset(support))


@lru_cache(maxsize=None)
def geometry_controls():
    rows = []
    for factor in (3, 4, 6):
        for count in (1, 2, 5):
            side = (3*count+1)*factor
            stars = [star_geometry(factor, (3*factor*i, 0, 0), side) for i in range(count)]
            assert all(left['support'].isdisjoint(right['support'])
                       for i, left in enumerate(stars) for right in stars[i+1:])
            assert all(len(star['boundary']) == 4*factor and len(star['faces']) == 16*factor-4
                       for star in stars)
            assert all(len(star['support']) <= 52*factor-16 for star in stars)
            rows.append(dict(block_factor=factor, readouts=count, fine_side=side,
                             each_star_links=len(stars[0]['support']), disjoint=True))
    # Independent comparison with the inherited complete lattice enumeration.
    geometry = previous.adapted.previous.blocked_geometry(4, 3)
    word = geometry['coarse']['faces'][0]
    active = {e for link, sign in word for e in geometry['paths'][link]}
    faces = [{e for e, sign in face} for face in geometry['fine']['faces']
             if active.intersection(e for e, sign in face)]
    actual = {geometry['fine']['edges'][e] for e in active.union(*faces)}
    assert actual == star_geometry(3, (0, 0, 0), 12)['support']
    return dict(growing_families=rows, native_full_lattice_comparison=True,
                theorem_for_arbitrary_sizes_uses_the_written_separation_argument=True)


def spectral_fixture(damping, order=None):
    damping = nonnegative(damping)
    if order is not None:
        native.width(order, 0)
    I = a.identity(2)
    U = previous.adapted.matrix(((1, 1), (1, 1)))
    V = a.diag((1, 0))
    epsilon = F(1, 10000)
    scalar = lambda rate: 1/(rate+damping) if order is None else F(rate)**(-(order+1))
    return a.add(a.add(a.scale(I, scalar(3)/12), U, epsilon*scalar(4)), V, epsilon*scalar(5))


def matrix_controls():
    error = F(36, 10000)
    normalized = []
    for s in (F(0), F(1, 7), F(5), F(1000)):
        normalized.append(a.scale(spectral_fixture(s), 1/limit_profile(s)))
        assert a.psd(a.add(normalized[-1], a.identity(2), -1))
        assert a.psd(a.add(a.scale(a.identity(2), 1+error), normalized[-1], -1))
    commutator = a.add(a.multiply(normalized[0], normalized[2]),
                      a.multiply(normalized[2], normalized[0]), -1)
    assert commutator[0][1] and abs(commutator[0][1]) <= 2*error*error
    stencil_checks = 0
    from math import comb, factorial
    for order in range(5):
        moment = spectral_fixture(0, order)
        for step in (F(1, 10), F(1, 2), F(2)):
            total = a.scale(a.identity(2), 0)
            scalar_intervals = []
            for k in range(order+1):
                total = a.add(total, spectral_fixture(k*step), (-1)**k*comb(order, k))
                phi = limit_profile(k*step)
                scalar_intervals.append((phi, (1+error)*phi))
            total = a.scale(total, 1/(factorial(order)*step**order))
            correction = previous.product(1+k*step/3 for k in range(1, order+1))
            assert a.psd(a.add(moment, total, -1))
            assert a.psd(a.add(a.scale(total, correction), moment, -1))
            interval = previous.stencil_interval(scalar_intervals, step, 3)
            assert a.psd(a.add(moment, a.identity(2), -interval['lower']))
            assert a.psd(a.add(a.scale(a.identity(2), interval['upper']), moment, -1))
            stencil_checks += 1
    return dict(fixture='Positive noncommuting matrix spectral measure; not the Yang--Mills vacuum',
                relative_enclosure_error=error, nonzero_normalized_commutator=commutator,
                matrix_stieltjes_stencil_checks=stencil_checks)


def profile_controls():
    rows = []
    for factor in (3, 4, 16, 64, 256):
        for z in (F(0), F(1), F(7, 3), F(1000)):
            value = leading.atom_response(leading.leading_atoms(factor), factor*z)[0][0]
            phi, correction = limit_profile(z), finite_block_correction(z)
            assert F(6*factor+3, 6*factor+8)*phi <= value <= (1+F(1, 2*factor))*phi
            remainder = value-phi-correction/factor
            c = z+3
            signed_sum = sum((delta*(aa*delta-bb*c)/(factor*factor*c*c*(c+delta/factor))
                             for aa, bb, delta in ((F(0), F(1, 16), F(0)),
                                (F(1, 16), -F(1, 32), F(3, 2)),
                                (F(1, 48), -F(1, 96), F(7, 2)),
                                (F(0), F(1, 48), F(4)))), F(0))
            assert remainder == signed_sum
            assert abs(remainder) <= phi*F(43, 36*factor*factor)
            rows.append(dict(block_factor=factor, spectral=z, exact_profile=value,
                             limit=phi, first_correction=correction, remainder=remainder))
    return dict(exact_free_profile_checks=rows, normalized_total_free_residue_limit=F(1, 12),
                limiting_rate=3, full_interacting_generator_has_only_this_rate=False)


def counterexample_controls():
    moments = []
    for target_order in (1, 2, 3):
        for factor in (4, 16, 64):
            delta, weight = GAP/factor, F(1, factor**target_order)
            atoms = ((F(3), (1-weight)/12), (delta, weight*delta/36))
            for z in (F(0), F(1, 7), F(1000)):
                value = leading.atom_response(atoms, z)[0][0]
                assert (1-weight)*limit_profile(z) <= value <= limit_profile(z)
            observed = leading.atom_moments(atoms, target_order+1)[target_order][0][0]
            baseline = F(1, 12*3**(target_order+1))
            excess = 1/(36*GAP**target_order)
            assert observed == (1-weight)*baseline+excess
            moments.append(dict(block_factor=factor, target_order=target_order,
                                relative_response_error=weight, first_uncontrolled_moment=observed,
                                nonvanishing_moment_excess=excess, hidden_floor=delta))
    relaxation = []
    for k in (2, 4, 6):
        delta = F(1, 2**k)
        tilt = delta
        soft = 2*tilt/(1+tilt*tilt)
        visible = (1-tilt*tilt)/(1+tilt*tilt)
        assert soft*soft+visible*visible == 1
        residual = (3-delta)*soft
        # t = 2^k log(2), making both heat factors exact rationals.
        target = F(1, 2**(3*2**k))
        correlation = soft*soft/2+visible*visible*target
        assert abs(correlation-target) <= residual*(1-target)/3
        assert soft*soft <= F(4, 9)*residual*residual
        relaxation.append(dict(full_gap=delta, time_in_log_two_units=2**k,
                               observed_correlation=correlation, rate_three_correlation=target,
                               generator_residual=residual, low_energy_weight=soft*soft,
                               zero_tilt_gives_exact_rate_three_at_the_same_full_gap=True))
    return dict(fixtures_are_positive_spectral_models_not_Yang_Mills_solutions=True,
                higher_moment_budget_is_necessary_for_an_inference_from_these_bounds=moments,
                visible_rate_does_not_determine_full_gap=relaxation)


def budget_controls():
    assert spatial.parameters(F(1, 10**6))['gamma'].lo > GAP
    endpoint = envelope(F(1, 10**8))
    assert endpoint['relative_lower_error'] < F(5396, 10**6)
    assert endpoint['relative_upper_error'] < F(10608, 10**6)
    # These large sizes evaluate proved scalar constants, not ambient matrices.
    comparisons = []
    for factor, count in ((3, 1), (4, 2), (16, 4), (64, 16), (10**5, 10**3)):
        theta = F(1, 10**8*factor*count)
        bound = parameters(theta, factor, count)
        direct = matrix_from_scalar_budget(theta, factor, count)
        assert direct['lower_error'] <= bound['relative_lower_error']
        assert direct['upper_error'] <= bound['relative_upper_error']
        comparisons.append(dict(block_factor=factor, readouts=count, theta=theta,
                                direct_errors=direct, universal_error=bound['relative_upper_error']))
    shapes = []
    for factor, count in ((16, 4), (64, 16), (256, 64)):
        theta = F(1, 10**8*count*factor*factor)
        shapes.append(dict(theta=theta, block_factor=factor, readouts=count,
                           **shape_budget(theta, factor, count), dynamics=dynamics_budget(theta, factor, count)))
    moments = []
    for factor in (4, 16, 64):
        count = factor
        theta = F(1, 10**24*count*factor**4)
        for order in (1, 2):
            moments.append(dict(theta=theta, block_factor=factor, readouts=count,
                                **normalized_moment_interval(theta, factor, count, order)))
    assert dynamics_budget(F(1, 3*10**8), 3, 1)['all_nonnegative_times_correlation_error'] < F(2134, 10**10)
    return dict(universal_product_window=dict(max_rho=F(1, 10**8), lower=F(99, 100),
                upper=F(51, 50), every_integer_block_factor_at_least_three=True,
                no_upper_limit_on_readout_count=True, constants=endpoint),
                direct_vs_universal= comparisons, growing_shapes=shapes,
                actual_matrix_moment_intervals=moments)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
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
    return spatial.encode(dict(
        certificate_type='YM72_GROWING_BLOCK_MATRIX_MEMORY_AND_VISIBLE_RELAXATION',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        geometry=geometry_controls(), matrices=matrix_controls(), profiles=profile_controls(),
        counterexamples=counterexample_controls(), budgets=budget_controls(),
        model=dict(actual_YM62_vacuum_and_YM64_spatial_cut_retained=True,
                   matrix_memory_for_separated_Wilson_readouts=True,
                   growing_blocks_and_readout_counts_under_product_budget=True,
                   uniform_in_admissible_finite_volume=True,
                   joint_normalized_one_pole_memory_limit=True,
                   higher_moment_order_has_its_own_scaling_budget=True,
                   actual_matrix_moments_enclosed=True,
                   actual_normalized_loop_time_correlations_enclosed=True,
                   full_gap_vanishes_on_the_declared_block_clock=True,
                   fixed_nonzero_coupling_growing_block_accuracy=False,
                   unnormalized_memory_has_nonzero_limit_on_block_clock=False,
                   full_spatial_memory_distance_decay=False,
                   infinite_order_uniformity=False,
                   larger_bare_coupling_window=False,
                   continuum_R4_QFT=False, physical_mass_prediction=False),
        claim_status='JOINT_COUPLING_BLOCK_READOUT_CONTROL__VISIBLE_RATE_IS_NOT_THE_FULL_GAP',
        proof_scope='Written matrix and scaling proofs with directed constants and exact native/algebraic controls',
        open='Persistent interaction, slow-mode retention, spatial memory decay and the physical continuum trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM72 fresh certificate/pin mismatch')
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
    print('YM72 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], product_budget='abs(theta)*r*N <= 1/100000000',
                          relative_matrix_bounds=['99/100', '51/50'],
                          fixed_coupling_growing_block_control=False, continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
