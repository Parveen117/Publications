"""YM68: exterior-size-independent averaged conditioning and native scores.

Python 3.12 only. Exact arithmetic checks accompany the written proof;
they do not sample or discretize the interacting vacuum. --check/default
are read-only. --write changes this stage's result and hash only.
"""
import argparse
from fractions import Fraction as F
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym67_pointwise_coarse_vacuum as previous

native, p, y, spatial = previous.native, previous.p, previous.y, previous.spatial
exact, Iv = previous.exact, previous.Iv
RESULT = HERE / 'YM68_RESULT.json'
PIN = HERE / 'EXPECTED_YM68.sha256'
SOURCES = HERE / 'YM68_SOURCE_PINS.json'


def smoothing_budget(eta, target_links, covariance):
    """Averaged L1 bound; smoothing dimension is the target dimension only."""
    cell = previous.parameters(eta)
    native.width(target_links)
    covariance = spatial.nonnegative(covariance)
    if covariance > 1:
        raise ValueError('Use a covariance bound in [0,1]')
    if not cell['eta'] or not covariance:
        return dict(bound=F(0), radius=None, method='exact_zero')
    K = previous.power_upper(cell['kappa'].hi, target_links)
    curvature = target_links * cell['density_second_ratio'].hi * K
    cap = min(F(2), previous.up(K - 1 / K))
    candidates = [(cap, -1)]
    # Finite dyadic choices certify the displayed tables. The written
    # asymptotic theorem takes the infimum over ALL radii in (0,1].
    for order in range(81):
        r = F(1, 2**order)
        value = 40**target_links * covariance / r**(3*target_links) + curvature*r*r
        candidates.append((previous.up(value), order))
    value, order = min(candidates)
    return dict(bound=value, radius=None if order == -1 else F(1, 2**order),
                method='density_range' if order == -1 else 'symmetric_target_ball')


def exterior_budget(eta, distance, target_links, block_factor):
    """No exterior size or exterior density floor occurs in this interface."""
    cell = previous.parameters(eta)
    native.width(distance)
    native.width(target_links)
    native.width(block_factor)
    delta = previous.covariance_upper(eta, distance, block_factor*target_links)
    smooth = smoothing_budget(eta, target_links, delta)
    error = smooth['bound']
    log_cap = 8*target_links*spatial.pi_interval().hi*cell['wave_log_gradient'].hi
    a1, ell2 = cell['score'].hi, cell['log_second'].hi
    fisher_coefficient = 6*target_links*ell2 + 2*target_links*a1*a1
    fisher_cap = 4*target_links*a1*a1
    return dict(distance=distance, target_links=target_links, block_factor=block_factor,
                fine_target_links=block_factor*target_links, covariance=delta,
                smoothing=smooth, expected_L1=error, expected_TV=previous.up(error/2),
                mutual_information=previous.up(min(log_cap, log_cap*error)),
                expected_relative_Fisher=previous.up(min(fisher_cap, fisher_coefficient*error)),
                Fisher_coefficient=previous.up(fisher_coefficient),
                analytic_exponent=F(2, 3*target_links+2))


def exceptional_probability(mean_bound, threshold):
    mean_bound, threshold = map(spatial.nonnegative, (mean_bound, threshold))
    if not threshold:
        raise ValueError('A positive threshold is required')
    return min(F(1), previous.up(mean_bound/threshold))


def gamma_poly(f, g):
    return p.add(*(p.mul(y.native_generator(f, a), y.native_generator(g, a))
                   for a in range(1, 4)))


def laplace_poly(f):
    """Negative-spectrum Laplacian Delta, not the positive Casimir L."""
    return p.add(*(y.native_generator(y.native_generator(f, a), a)
                   for a in range(1, 4)))


def one_poly():
    return {(0, 0, 0, 0): F(1)}


def native_controls():
    X = p.COORD[0]
    assert laplace_poly(X) == p.scale(X, F(-3, 4))
    gamma = gamma_poly(X, X)
    assert y.phi(gamma) == F(3, 16)
    first = F(1, 3)
    q = p.add(one_poly(), p.scale(p.COORD[1], first))
    trials = [X, p.mul(X, X), p.mul(X, p.COORD[2])]
    weak_checks = 0
    for u in trials:
        for f in trials:
            r_minus_one = p.add(u, p.scale(one_poly(), -y.phi(p.mul(q, u))))
            left = y.phi(p.mul(q, gamma_poly(r_minus_one, f)))
            div = p.add(p.mul(q, laplace_poly(f)), gamma_poly(q, f))
            right = -y.phi(p.mul(r_minus_one, div))
            assert left == right
            weak_checks += 1
    # Symmetric geodesic quadrature cancels the linear term exactly.
    # t -> q exp(-t u/2) has scalar second derivative -q0/4.
    directions = [(F(1), F(0), F(0)), (F(3,5), F(4,5), F(0))]
    for direction in directions:
        second = previous.derivative(previous.derivative(X, 0, direction), 0, direction)
        assert second == p.scale(X, F(-1, 4))
    # An even quadratic attains the 1/2 in the symmetric Taylor bound.
    paired = [(F(-1), F(1,4)), (F(0), F(1,2)), (F(1), F(1,4))]
    assert sum(t*w for t,w in paired) == 0
    assert sum(t*t*w for t,w in paired) == F(1,2)
    return dict(weighted_integration_by_parts_checks=weak_checks,
                native_coordinate_rate=F(3,4), native_coordinate_energy=F(3,16),
                ordered_unit_second_derivative_checks=len(directions),
                symmetric_linear_cancellation=True, Taylor_remainder_factor=F(1,2),
                target_real_dimension_per_link=3)


def coordinate_even_moment(order):
    native.width(order, 0)
    return F(comb(2*order, order), (order+1)*4**order)


def fisher_series(amplitude, terms):
    """Independent native integral for p=1+t q0, q=1; rigorous tail."""
    t = spatial.nonnegative(amplitude)
    native.width(terms)
    if t >= 1:
        raise ValueError('Require a strictly positive affine density')
    lower = t*t/4*sum(t**(2*k)*(coordinate_even_moment(k)-coordinate_even_moment(k+1))
                      for k in range(terms))
    remainder = t**(2*terms+2)/(4*(1-t*t))
    return Iv(lower, lower+remainder)


def fisher_controls():
    rows = []
    # Pythagorean parameters make sqrt(1-t^2), and thus the exact
    # semicircle resolvent integral, rational.
    for t, root in ((F(3,5), F(4,5)), (F(20,101), F(99,101)),
                    (F(200,10001), F(9999,10001))):
        G = 2/(1+root)
        value = (1+(t*t-1)*G)/4
        enclosure = fisher_series(t, 80)
        assert enclosure.lo <= value <= enclosure.hi
        score = t/(2*(1-t))
        second = t/(4*(1-t)) + t*t/(4*(1-t)**2)
        coefficient = 6*second + 2*score*score
        # E|q0| = 4/(3*pi). This exact native marginal provides
        # an independent non-polynomial Fisher/variation comparison.
        L1_upper = t*4/(3*spatial.pi_interval().lo)
        assert value <= coefficient*L1_upper
        rows.append(dict(amplitude=t, exact_Fisher=value, series_enclosure=enclosure,
                         L1_upper=L1_upper, comparison_upper=coefficient*L1_upper))
    return dict(scope='Affine native densities, not the interacting YM vacuum', cases=rows)


def rare_pin_case(pins, epsilon=F(1,100)):
    native.width(pins)
    epsilon = spatial.nonnegative(epsilon)
    if not 0 < epsilon < 1:
        raise ValueError('Require epsilon in (0,1)')
    mean = F(1,2)**pins
    variance = F(5,16)**pins-F(1,4)**pins
    return dict(pins=pins, epsilon=epsilon, hidden_mean=mean, hidden_variance=variance,
                average_L1_upper=2*epsilon*mean,
                chi_square=epsilon**2*variance/4,
                Fisher_upper=3*epsilon**2*variance/(16*(1-epsilon)),
                exact_pin_density_amplitude=epsilon*(1-mean),
                exact_pin_score=epsilon*(1-mean)/2,
                density_floor=1-epsilon)


def limitation_controls():
    a = p.scale(p.add(one_poly(), p.COORD[0]), F(1,2))
    assert y.phi(a) == F(1,2)
    assert y.phi(p.mul(a,a)) == F(5,16)
    assert y.phi(p.mul(p.COORD[0],p.COORD[0])) == F(1,4)
    rows = [rare_pin_case(n) for n in (1,2,8,32,128)]
    epsilon = F(1,100)
    # On Q^3, 1+epsilon*x*y*z has independent (x,z) marginal,
    # but after y=1 the x,z covariance is epsilon/16.
    covariance = epsilon*F(1,4)**2
    return dict(scope='Positive smooth comparison densities, not YM counterexamples',
                rare_pins=rows, nearby_conditioning=dict(
                    unconditioned_target_remote_covariance=F(0),
                    covariance_given_middle_coordinate_one=covariance,
                    score_at_target_i_remote_one=epsilon/2,
                    score_at_target_i_remote_minus_one=-epsilon/2))


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
    eta = F(1,16000)
    budgets = [exterior_budget(eta, d, 4, 2) for d in (64,128,256,512,1024)]
    return spatial.encode(dict(
        certificate_type='YM68_AVERAGED_EXTERIOR_LOCALITY_AND_NATIVE_FISHER_CONTROL',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        native=native_controls(), Fisher=fisher_controls(), limitations=limitation_controls(),
        exterior_budgets=budgets,
        exceptional_L1_probability_at_512=exceptional_probability(budgets[3]['expected_L1'],F(1,1000)),
        model=dict(actual_YM62_vacuum_retained=True, unbounded_link_spaces_retained=True,
                   averaged_bounds_independent_of_exterior_cardinality=True,
                   countably_infinite_exterior_conditionals_C1_almost_surely=True,
                   averaged_relative_Fisher_decay=True, averaged_target_diagonal_drift_decay=True,
                   coarse_spatial_tail_trivial=True,
                   worst_case_exterior_configuration_locality=False,
                   uniform_mixing_after_nearby_conditioning=False,
                   full_Schur_memory_kernel_locality=False, time_zero_Markov_identified=False,
                   continuum_R4_QFT=False, physical_mass_prediction=False),
        claim_status='AVERAGED_INFINITE_EXTERIOR_CONTROL__WORST_CASE_AND_SCHUR_MEMORY_OPEN',
        proof_scope='Written proofs from YM65-YM67 with exact and directed arithmetic checks'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM68 fresh certificate/pin mismatch')
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
    print('YM68 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], averaged_infinite_exterior=True,
                          worst_case_exterior=False, full_memory_locality=False,
                          continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
