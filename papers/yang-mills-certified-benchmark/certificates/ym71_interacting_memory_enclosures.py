"""YM71: volume-uniform enclosures of the actual interacting Wilson memory.

Python 3.12 only. Written native form proofs supply the inequalities;
directed rational calculations enclose their constants and moment inputs.
No vacuum discretization or continuum construction is claimed.
Default/--check are read-only; --write changes this result and hash only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from math import comb, factorial, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym70_certified_spectral_memory as previous
import ym67_pointwise_coarse_vacuum as pointwise

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
adapted, spatial, model = previous.previous, pointwise.spatial, previous.model
RESULT = HERE/'YM71_RESULT.json'
PIN = HERE/'EXPECTED_YM71.sha256'
SOURCES = HERE/'YM71_SOURCE_PINS.json'
GRID = 10**36


def up(value):
    value = exact(value)
    return F(-((-value.numerator*GRID)//value.denominator), GRID)


def sqrt_upper(value):
    value = previous.nonnegative(value)
    if not value:
        return F(0)
    floor = isqrt(value.numerator*GRID*GRID//value.denominator)
    return F(floor if F(floor, GRID)**2 == value else floor+1, GRID)


def exponential_upper(value):
    value = previous.nonnegative(value)
    if not value:
        return F(1)
    halvings = 0
    while value > 1:
        value /= 2
        halvings += 1
    result = up(pointwise.exp_interval(value).hi)
    for _ in range(halvings):
        result = up(result*result)
    return result


@lru_cache(maxsize=None, typed=True)
def parameters(theta, factor=3):
    theta = exact(theta)
    native.width(factor, 3)
    eta = abs(theta)
    gap = spatial.parameters(theta)['gamma'].lo
    B = 24*eta
    m, ell, support = 4*factor, 4*factor+3, 52*factor-16
    A = F(16*factor, 3)
    lo_rate, hi_rate = F(3*factor), F(3*factor+4)
    G2 = F(2*factor+1, 24)
    K = exponential_upper(4*spatial.pi_interval().hi*B*support)
    Gbar = sqrt_upper(K*G2)
    delta = up(F(4, 3)*B*A)
    # Every use of an irrational quantity is rounded upwards explicitly.
    residual = up(A*(2*B*ell+F(8, 3)*eta+8*B*B*ell
                     +8*sqrt_upper(2*eta*B*ell)))
    forcing_per_theta = up(delta+8*m*B/gap)
    hidden_norm = up(Gbar+delta)
    comparison = up((K-1)*G2+2*Gbar*delta+delta*delta)
    linear_error = up((hi_rate/lo_rate)*(comparison+2*forcing_per_theta*hidden_norm
                      +hidden_norm*residual/lo_rate)/G2)
    squared_error = up(max(F(1), hi_rate/gap)*(forcing_per_theta+residual/lo_rate)**2/G2)
    return dict(theta=theta, eta=eta, block_factor=factor, gap_lower=gap,
                wave_score_upper=B, fine_boundary_links=m, component_support_upper=ell,
                union_support_upper=support, absolute_word_weight=A,
                minimum_free_rate=lo_rate, maximum_free_rate=hi_rate,
                free_forcing_norm_squared=G2, marginal_density_ratio_upper=K,
                free_forcing_weighted_norm_upper=Gbar, lift_error=delta,
                eigen_equation_residual=residual, forcing_error_per_theta=forcing_per_theta,
                lifted_forcing_norm_upper=hidden_norm, inner_product_comparison=comparison,
                relative_lower_error=linear_error, relative_upper_error=up(linear_error+squared_error),
                residual_squared_error=squared_error)


def response_interval(theta, factor, damping, clock=F(1)):
    damping, clock = previous.nonnegative(damping), previous.positive(clock)
    cell = parameters(theta, factor)
    coefficient = previous.atom_response(previous.leading_atoms(factor), damping/clock)[0][0]
    center = clock*cell['eta']**2*coefficient
    lower = max(F(0), (1-cell['relative_lower_error'])*center)
    upper = (1+cell['relative_upper_error'])*center
    return dict(lower=lower, upper=upper, leading=center)


def stencil_interval(intervals, step, gap):
    """Positive inverse moments from signed finite differences of responses."""
    step, gap = previous.positive(step), previous.positive(gap)
    if not intervals:
        raise ValueError('At least one response interval required')
    cells = []
    for lower, upper in intervals:
        lower, upper = previous.nonnegative(lower), previous.nonnegative(upper)
        if lower > upper:
            raise ValueError('Reversed response interval')
        cells.append((lower, upper))
    order = len(cells)-1
    left = right = F(0)
    for k, (lower, upper) in enumerate(cells):
        weight = (-1)**k*comb(order, k)
        left += weight*(lower if weight >= 0 else upper)
        right += weight*(upper if weight >= 0 else lower)
    denominator = factorial(order)*step**order
    left, right = max(F(0), left/denominator), right/denominator
    if right < 0:
        raise ValueError('Intervals exclude a positive Stieltjes difference')
    correction = F(1)
    for k in range(1, order+1):
        correction *= 1+k*step/gap
    return dict(lower=left, upper=correction*right,
                divided_difference_lower=left, divided_difference_upper=right,
                support_correction=correction)


def direct_moment_radius(cell, order):
    """Separate inverse-power residual argument; no derivative of an inequality."""
    native.width(order, 0)
    power = order+1
    eta, gap, rate = cell['eta'], cell['gap_lower'], cell['minimum_free_rate']
    norm, residual, force = (cell['lifted_forcing_norm_upper'],
                             cell['eigen_equation_residual'], cell['forcing_error_per_theta'])
    telescope = sum((rate**(-(power-k+1))*gap**(-k) for k in range(1, power+1)), F(0))
    return eta*eta*(cell['inner_product_comparison']/rate**power
                   +norm*residual*telescope+(2*norm*force+force*force)*gap**(-power))


def moment_interval(theta, factor, order, clock=F(1)):
    native.width(order, 0)
    clock = previous.positive(clock)
    cell = parameters(theta, factor)
    leading = cell['eta']**2*previous.atom_moments(previous.leading_atoms(factor), order+1)[order][0][0]
    radius = direct_moment_radius(cell, order)
    zero = response_interval(theta, factor, 0)
    lower = max(F(0), leading-radius)
    upper = min(leading+radius, zero['upper']*cell['gap_lower']**(-order))
    best_step = None
    if order == 0:
        lower, upper = max(lower, zero['lower']), min(upper, zero['upper'])
    else:
        # Each candidate is a rigorous enclosure; intersect, do not fit to truth.
        for exponent in range(-28, 9):
            step = cell['minimum_free_rate']*F(2)**exponent
            rows = [response_interval(theta, factor, k*step) for k in range(order+1)]
            candidate = stencil_interval([(row['lower'], row['upper']) for row in rows],
                                         step, cell['gap_lower'])
            candidate_lower, candidate_upper = max(lower, candidate['lower']), min(upper, candidate['upper'])
            if candidate_upper-candidate_lower < upper-lower:
                best_step = step
            lower, upper = candidate_lower, candidate_upper
    if lower > upper:
        raise ArithmeticError('Independent moment enclosures do not overlap')
    scale = clock**(1-order)
    return dict(order=order, clock=clock, lower=scale*lower, upper=scale*upper,
                leading=scale*leading, direct_radius=scale*radius, last_improving_unscaled_step=best_step)


def integrate_link(poly, edge):
    """Exact native Haar integration in one link; retain all other coordinates."""
    out = {}
    for powers, coefficient in poly.items():
        local = adapted.y.phi({powers[4*edge:4*edge+4]: F(1)})
        key = powers[:4*edge]+powers[4*edge+4:]
        out[key] = out.get(key, F(0))+coefficient*local
    return p.clean(out)


@lru_cache(maxsize=None)
def native_controls():
    g = model.loop_polynomial()
    constant = model.source.constant(1, 4)
    rows = []
    for t in (F(1, 7), F(-1, 5)):
        h = p.add(constant, p.scale(g, t))
        for edge in range(4):
            coupling = p.add(*(p.mul(native.derivative(h, edge, axis),
                                     native.derivative(g, edge, axis)) for axis in (1, 2, 3)))
            numerator = p.add(p.mul(h, g), p.scale(coupling, -F(8, 3)))
            assert integrate_link(p.mul(h, numerator), edge) == {}
            raw = integrate_link(p.mul(p.mul(h, h), g), edge)
            assert raw  # The free function does not stay hidden for this weight.
            for axis in (1, 2, 3):
                dh = native.derivative(h, edge, axis)
                assert model.electric(dh, 4) == native.derivative(model.electric(h, 4), edge, axis)
            rows.append(dict(weight_parameter=t, hidden_link=edge,
                             weighted_lift_has_exact_zero_conditional_mean=True,
                             raw_free_component_has_nonzero_conditional_mean=True))
    geometry_rows = []
    for factor in (3, 4, 5):
        geometry = adapted.previous.blocked_geometry(4, factor)
        retained = {edge for path in geometry['paths'] for edge in path}
        word = geometry['coarse']['faces'][0]
        boundary = {edge for link, sign in word for edge in geometry['paths'][link]}
        faces = [{e for e, sign in face} for face in geometry['fine']['faces']
                 if boundary.intersection(e for e, sign in face)]
        assert all(face-retained for face in faces)
        support = boundary.union(*faces)
        assert len(support) <= 52*factor-16
        word_weight = sum((F(len(face & boundary), 8)+F(5*len(face & boundary), 24)
                           for face in faces), F(0))
        assert word_weight == F(16*factor, 3)
        geometry_rows.append(dict(block_factor=factor, incident_plaquettes=len(faces),
                                  actual_union_links=len(support), union_bound=52*factor-16,
                                  every_component_has_an_unretained_fundamental_link=True,
                                  absolute_word_weight=word_weight))
    return dict(weighted_fixture='h=1+t W_p, a positive gauge-invariant test weight, not the interacting vacuum',
                conditional_lifts=rows, actual_product_cut_geometry=geometry_rows,
                native_Casimir_derivative_commutations=24)


def dot(left, right):
    return sum((v*w for v, w in zip(left, right)), F(0))


def matvec(matrix, vector):
    return [dot(row, vector) for row in matrix]


def residual_controls():
    rates = (F(2), F(3), F(5))
    free = a.diag(rates)
    perturbation = adapted.matrix(((2, 1, 0), (1, 2, 1), (0, 1, 2)))
    mixing = adapted.matrix(((0, 1, 0), (-1, 0, 1), (0, -1, 0)))
    assert a.psd(perturbation)
    assert a.multiply(free, perturbation) != a.multiply(perturbation, free)
    rows = []
    power_checks = 0
    for t in (F(1, 11), F(1, 101)):
        D = a.add(free, perturbation, t)
        assert a.psd(a.add(D, a.identity(3), -2))
        columns = a.transpose(a.add(a.identity(3), mixing, t))
        H = [sum((column[i] for column in columns), F(0)) for i in range(3)]
        q = [t*value+t*t*F(i-1, 7) for i, value in enumerate(H)]
        for s in (F(0), F(1, 7), F(1000)):
            shifted = a.add(D, a.identity(3), s)
            x = [t*sum((column[i]/(rate+s) for column, rate in zip(columns, rates)), F(0))
                 for i in range(3)]
            residual = [v-w for v, w in zip(q, matvec(shifted, x))]
            J = 2*dot(q, x)-dot(x, matvec(shifted, x))
            sigma = dot(q, matvec(a.inverse(shifted), q))
            correction = dot(residual, matvec(a.inverse(shifted), residual))
            assert sigma-J == correction
            assert 0 <= correction <= dot(residual, residual)/(2+s)
            rows.append(dict(parameter=t, damping=s, variational_value=J,
                             exact_response=sigma, residual_correction=correction))
        for column, rate in zip(columns, rates):
            e = [left-rate*right for left, right in zip(matvec(D, column), column)]
            for j in range(1, 6):
                left = [v-rate**(-j)*w for v, w in zip(matvec(a.power(a.inverse(D), j), column), column)]
                right = [sum((-rate**(-(j-k+1))*matvec(a.power(a.inverse(D), k), e)[i]
                              for k in range(1, j+1)), F(0)) for i in range(3)]
                assert left == right
                power_checks += 1
    return dict(fixture='Positive noncommuting rational matrices; not a Yang--Mills truncation',
                energy_residual_identities=rows, inverse_power_telescope_checks=power_checks)


def stencil_controls():
    atoms = ((F(2), F(1, 7)), (F(3), F(2, 5)), (F(5), F(1, 3)))
    moments = previous.atom_moments(atoms, 7)
    rows = []
    for order in range(7):
        for step in (F(1, 100), F(1, 3), F(2)):
            responses = [previous.atom_response(atoms, k*step)[0][0] for k in range(order+1)]
            bound = stencil_interval([(v, v) for v in responses], step, 2)
            direct = sum((weight / product(rate+k*step for k in range(order+1))
                          for rate, weight in atoms), F(0))
            assert bound['divided_difference_lower'] == direct == bound['divided_difference_upper']
            actual = moments[order][0][0]
            assert bound['lower'] <= actual <= bound['upper']
            padded = stencil_interval([(v*(1-F(1, 10**12)), v*(1+F(1, 10**12)))
                                       for v in responses], step, 2)
            assert padded['lower'] <= actual <= padded['upper']
            rows.append(dict(order=order, step=step, exact_moment=actual,
                             finite_difference=direct, support_factor=bound['support_correction']))
    delta, tiny = F(1, 100), F(1, 10000)
    alternative = ((F(1), 1-delta), (tiny, delta*tiny))
    for s in (F(0), F(1, 7), F(1), F(1000)):
        value = previous.atom_response(alternative, s)[0][0]
        assert (1-delta)/(1+s) <= value <= 1/(1+s)
    first = previous.atom_moments(alternative, 2)[1][0][0]
    assert first > 100
    return dict(positive_measure_stencil_checks=rows,
                response_bounds_cannot_be_differentiated=dict(relative_value_width=delta,
                    baseline_first_moment=1, alternative_first_moment=first, hidden_gap=tiny))


def product(values):
    result = F(1)
    for value in values:
        result *= value
    return result


def budget_controls():
    endpoint = parameters(F(1, 10**8), 3)
    assert endpoint['relative_lower_error'] < F(16, 1000)
    assert endpoint['relative_upper_error'] < F(61, 1000)
    examples = []
    displays = ((12, 1, F(996, 1000), F(1133, 1000)),
                (16, 2, F(996, 1000), F(1226, 1000)),
                (20, 3, F(997, 1000), F(1029, 1000)))
    for exponent, order, lower, upper in displays:
        theta = F(1, 10**exponent)
        row = moment_interval(theta, 3, order)
        assert lower*row['leading'] <= row['lower'] <= row['upper'] <= upper*row['leading']
        examples.append(dict(theta=theta, **row, displayed_ratio_lower=lower, displayed_ratio_upper=upper))
    theta, clock, s = F(1, 10**12), F(7, 3), F(5)
    scaled, base = response_interval(theta, 3, s, clock), response_interval(theta, 3, s/clock)
    assert all(scaled[key] == clock*base[key] for key in ('lower', 'upper', 'leading'))
    for order in range(4):
        scaled, base = moment_interval(theta, 3, order, clock), moment_interval(theta, 3, order)
        assert all(scaled[key] == clock**(1-order)*base[key] for key in ('lower', 'upper', 'leading'))
    broad = parameters(F(1, 16000), 3)
    assert broad['relative_lower_error'] > 1 and broad['relative_upper_error'] > 1000
    return dict(nonzero_memory_window=dict(max_abs_theta=F(1, 10**8), block_factor=3,
                every_admissible_finite_volume=True, every_nonnegative_s=True,
                simple_relative_lower=F(49, 50), simple_relative_upper=F(107, 100), constants=endpoint),
                interacting_moment_intervals=examples,
                broad_window_is_not_a_small_relative_error=dict(theta=F(1, 16000),
                    relative_lower_error=broad['relative_lower_error'], relative_upper_error=broad['relative_upper_error']),
                clock_covariance_checks=5)


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
        certificate_type='YM71_ACTUAL_INTERACTING_WILSON_MEMORY_WITH_VOLUME_UNIFORM_ERRORS',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        native=native_controls(), residuals=residual_controls(), stencils=stencil_controls(),
        budgets=budget_controls(),
        model=dict(actual_YM62_vacuum_and_YM64_spatial_cut_retained=True,
                   exact_interacting_hidden_trial_vectors=True,
                   volume_uniform_fixed_block_leading_remainder=True,
                   actual_scalar_interacting_response_enclosed=True,
                   actual_scalar_inverse_moments_enclosed=True,
                   nonzero_memory_on_explicit_coupling_interval=True,
                   positive_parameter_axis_uniform_relative_enclosure=True,
                   full_vacuum_numerically_computed=False,
                   exactly_four_rates_for_full_memory=False,
                   stable_interval_Gram_inversion_proved=False,
                   full_spatial_memory_locality=False,
                   uniform_as_block_factor_grows=False,
                   enlarged_bare_coupling_window=False,
                   continuum_R4_QFT=False, physical_mass_prediction=False),
        claim_status='FIXED_BLOCK_ACTUAL_MEMORY_AND_MOMENT_ENCLOSURES_UNIFORM_IN_FINITE_VOLUME',
        proof_scope='Written native score and hidden-lift proof with directed constants and separate algebraic controls',
        open='Coupled spatial readouts, growing blocks, generated interactions and the continuum trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM71 fresh certificate/pin mismatch')
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
    print('YM71 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], block_factor=3,
                          max_abs_theta='1/100000000', all_s_relative_bounds=['49/50', '107/100'],
                          uniform_in_finite_volume=True, continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
