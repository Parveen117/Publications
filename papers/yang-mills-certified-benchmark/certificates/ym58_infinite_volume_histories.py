"""YM58 exact marginal/history controls. Python 3.12 only.

Reuses the native polynomial engine. Finite controls are distinguished
from the written infinite-volume proof. --check is read-only.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym57_neighbour_turn_gap as prior

a, p, q, y = prior.a, prior.p, prior.q, prior.y
exact = prior.exact
RESULT = HERE / 'YM58_RESULT.json'
PIN = HERE / 'EXPECTED_YM58.sha256'
SOURCES = HERE / 'YM58_SOURCE_PINS.json'


def sites(values):
    values = tuple(values)
    if any(isinstance(i, bool) or not isinstance(i, int) for i in values):
        raise ValueError('Integer lattice sites required')
    if values != tuple(sorted(set(values))):
        raise ValueError('Sites must be sorted and distinct')
    return values


def rate(k):
    k = exact(k)
    if k < 0:
        raise ValueError('Nonnegative edge rate required')
    return k


def polynomial(poly, support):
    support = sites(support)
    for m, c in poly.items():
        if len(m) != 4 * len(support) or any(
                isinstance(v, bool) or not isinstance(v, int) or v < 0 for v in m):
            raise ValueError('Invalid packed native exponent')
        exact(c)
    return p.clean(poly)


def embed(poly, source, target):
    source, target = sites(source), sites(target)
    poly = polynomial(poly, source)
    if not set(source) <= set(target):
        raise ValueError('Target must contain source sites')
    positions = [target.index(i) for i in source]
    out = {}
    for m, c in poly.items():
        term = [0] * (4 * len(target))
        for old, new in enumerate(positions):
            term[4*new:4*new+4] = m[4*old:4*old+4]
        out[tuple(term)] = c
    return p.clean(out)


def incident_edges(support):
    support = sites(support)
    return tuple(sorted({(i-1, i) for i in support} | {(i, i+1) for i in support}))


def active_weight(support, k):
    support, k = sites(support), rate(k)
    return len(support) + k * len(incident_edges(support))


def directions(support, k):
    """Weighted derivatives; a missing endpoint contributes no derivative."""
    support, k = sites(support), rate(k)
    positions = {site: i for i, site in enumerate(support)}
    rows = [(F(1), ((i, 1, F(1)),)) for i in range(len(support))]
    for left, right in incident_edges(support):
        word = []
        for site, sign in ((left, F(1)), (right, F(-1))):
            if site in positions:
                word.append((positions[site], 2, sign))
        rows.append((k, tuple(word)))
    return rows


def apply_direction(poly, word):
    return p.add(*(p.scale(prior.derivative(poly, i, axis), sign)
                   for i, axis, sign in word))


def square_generator(poly, support, k):
    poly = polynomial(poly, support)
    return p.scale(p.add(*(p.scale(apply_direction(apply_direction(poly, word), word), w)
                           for w, word in directions(support, k))), F(-1, 2))


def marginal(poly, support, k):
    """Expanded generator, independently checked against incident-edge squares."""
    support, k = sites(support), rate(k)
    poly = polynomial(poly, support)
    terms = []
    for i in range(len(support)):
        terms.append(p.scale(prior.derivative(prior.derivative(poly, i, 1), i, 1), F(-1, 2)))
        terms.append(p.scale(prior.derivative(prior.derivative(poly, i, 2), i, 2), -k))
    for i in range(len(support)-1):
        if support[i+1] == support[i]+1:
            terms.append(p.scale(prior.derivative(prior.derivative(poly, i, 2), i+1, 2), k))
    return p.add(*terms)


def power(poly, support, k, count, operator_fn=None):
    prior.width(count, 0)
    operator_fn = operator_fn or (lambda f: marginal(f, support, k))
    for _ in range(count):
        poly = operator_fn(poly)
    return poly


def history_coefficient(observables, orders, support, k, operator_fn=None):
    """Coefficient of the independent time-gap powers, not a heat approximation."""
    support = sites(support)
    if not observables or len(orders) != len(observables)-1:
        raise ValueError('One order per history time gap required')
    for order in orders:
        prior.width(order, 0)
    for f in observables:
        polynomial(f, support)
    current = observables[-1]
    for j in reversed(range(len(orders))):
        r = orders[j]
        current = p.scale(power(current, support, k, r, operator_fn), F((-1)**r, factorial(r)))
        current = p.mul(observables[j], current)
    return prior.phi(current, len(support))


def refinement_budget(support, k, volume_width, steps, degree, time):
    support, k, time = sites(support), rate(k), exact(time)
    prior.width(volume_width)
    prior.width(steps)
    prior.width(degree, 0)
    if time < 0:
        raise ValueError('Nonnegative heat time required')
    if support and volume_width < support[-1]-support[0]+3:
        raise ValueError('Volume too small to contain the neighbour layer')
    W = volume_width + k*(volume_width-1)
    A = active_weight(support, k)
    return time*time*degree**4 / steps * (W*A/F(384) + A*A/F(128))


def fast_volume_multiplier_lower(steps):
    prior.width(steps, 2)
    return F(1)-F(6*steps, 2*steps*steps-1)


def boundary_controls():
    rows = []
    w = prior.previous.witness(1)
    for k in (F(1, 10), F(1, 2), F(1), F(3)):
        assert marginal(w, (0,), k) == p.scale(w, k)
        assert not prior.operator(w, 1, k)
        two_site = prior.lift(w, 0, 2)
        assert prior.operator(two_site, 2, k) == p.scale(two_site, k/2)
        for f in p.COORD:
            assert marginal(f, (0,), k) == p.scale(f, (1+2*k)/8)
        expected_gap = min((1+2*k)/8, min(F(1, 2), k))
        inherited = prior.previous.energy.rates(tuple(sorted((F(0), F(1, 2), k))))[0]
        assert expected_gap == inherited
        rows.append(dict(k=k,bulk_w_rate=k,endpoint_w_rate=k/2,
                         singleton_gap=expected_gap,linear_rate=(1+2*k)/8))
    assert marginal({(): F(7)}, (), 1) == {}
    assert not marginal(w, (0,), 0)
    return dict(rows=rows,missing_boundary_noise_refused=True,
                zero_edge_rate_has_centered_invariant=True)


def locality_controls():
    cells = (((0,), (2,)), ((0, 1), (1, 1)), ((0, 2), (1, 1)),
             ((-3, 1), (2, 0)), ((0, 1, 2), (1, 1, 0)))
    generator_checks = nested_checks = power_checks = 0
    for support, degrees in cells:
        _, polys, _, _ = prior.block(degrees)
        halo = tuple(range(support[0]-1, support[-1]+2))
        larger = tuple(range(support[0]-2, support[-1]+4))
        intermediate = tuple(sorted(set(support) | {support[0]-1}))
        for k in (F(1, 7), F(1), F(3)):
            for f in polys:
                h = marginal(f, support, k)
                assert h == square_generator(f, support, k)
                for volume in (halo, larger):
                    full = prior.operator(embed(f, support, volume), len(volume), k)
                    assert full == embed(h, support, volume)
                    generator_checks += 1
                assert marginal(embed(f, support, intermediate), intermediate, k) == embed(h, support, intermediate)
                nested_checks += 1
            # All coefficients of the heat series agree by induction; these
            # are finite exact controls of that induction, not a time cutoff.
            f = p.add(polys[0], p.scale(polys[-1], F(2, 3)))
            for r in range(6):
                local = power(f, support, k, r)
                full = power(embed(f, support, larger), larger, k, r,
                             lambda g: prior.operator(g, len(larger), k))
                assert embed(local, support, larger) == full
                power_checks += 1
    return dict(open_volume_generator_identities=generator_checks,
                arbitrary_subset_consistency_identities=nested_checks,
                power_identities=power_checks,largest_power=5,
                support_sets=[list(s) for s, _ in cells],
                all_times_proof='YM58-T1--T2 finite invariant blocks and density')


def energy_controls():
    cells = (((0,), (2,)), ((0, 1), (1, 1)), ((0, 1), (2, 1)),
             ((0, 2), (1, 1)), ((-2, 1), (1, 2)))
    inequalities = factorizations = invariance = 0
    for support, degrees in cells:
        mons, polys, gram, centered = prior.block(degrees)
        n = len(support)
        for k in (F(1, 10), F(1), F(3)):
            columns = [marginal(f, support, k) for f in polys]
            H = [[column.get(m, F(0)) for column in columns] for m in mons]
            form = a.multiply(gram, H)
            assert form == a.transpose(form)
            assert a.psd(a.add(form, a.scale(centered, -prior.lower_bound(2, k))))
            inequalities += 1
            factored = a.scale(gram, 0)
            for weight, word in directions(support, k):
                images = [apply_direction(f, word) for f in polys]
                factor = [[weight*prior.phi(p.mul(v, w), n)/2 for w in images] for v in images]
                factored = a.add(factored, factor)
            assert factored == form
            factorizations += 1
            for col in columns:
                assert prior.phi(col, n) == 0
                invariance += 1
    return dict(rational_gap_inequalities=inequalities,
                independent_incident_square_factorizations=factorizations,
                invariant_reference_identities=invariance,
                max_block_dimension=40,unit_edge_gap_lower=F(1, 630),
                all_degree_infinite_volume_proof='YM58-T4 transfer and form closure')


def history_controls():
    coefficients = nonzero = 0
    for support in ((0, 1), (0, 2)):
        first = prior.lift(p.COORD[0], 0, 2)
        second = prior.lift(p.COORD[2], 1, 2)
        observables = (first, second, p.mul(first, second))
        for k in (F(1, 2), F(1)):
            for orders in itertools.product(range(4), repeat=2):
                local = history_coefficient(observables, orders, support, k)
                nonzero += bool(local)
                for padding in (1, 2):
                    volume = tuple(range(support[0]-padding, support[-1]+padding+1))
                    lifted = [embed(f, support, volume) for f in observables]
                    value = history_coefficient(lifted, orders, volume, k,
                                                lambda f: prior.operator(f, len(volume), k))
                    assert value == local
                    coefficients += 1
    assert nonzero > 0
    f = p.mul(prior.lift(p.COORD[0], 0, 2), prior.lift(p.COORD[0], 1, 2))
    g = p.mul(prior.lift(p.COORD[2], 0, 2), prior.lift(p.COORD[2], 1, 2))
    for k in (F(1, 10), F(1), F(3)):
        plus, minus = p.add(f, g), p.add(f, p.scale(g, -1))
        assert marginal(plus, (0, 1), k) == p.scale(plus, (1+3*k)/4)
        assert marginal(minus, (0, 1), k) == p.scale(minus, (1+k)/4)
        assert marginal(f, (0, 2), k) == p.scale(f, (1+2*k)/4)
    assert prior.phi(p.mul(f, f), 2) == prior.phi(p.mul(g, g), 2) == F(1, 16)
    assert prior.phi(p.mul(f, g), 2) == 0
    # k=1, t=4 log 2. Exact eigenmodes turn these heat values rational.
    adjacent = (F(1, 2)**4 + F(1, 2)**2)/32
    separated = F(1, 2)**3/16
    assert adjacent == F(5, 512) and separated == F(1, 128)
    assert adjacent > separated
    # A singleton's linear heat correlation at t=log(2)/rho has r=1/2.
    r = F(1, 2)
    grid = (-2, -1, 0, 1, 2)
    stationary_gram = [[r**abs(i-j)/4 for j in grid] for i in grid]
    reflected_gram = [[r**(i+j)/4 for j in range(4)] for i in range(4)]
    assert a.psd(stationary_gram) and a.psd(reflected_gram)
    assert reflected_gram == [[(r**i/2)*(r**j/2) for j in range(4)] for i in range(4)]
    assert not a.psd([[F(1, 4), F(1, 2)], [F(1, 2), F(1, 4)]])
    return dict(multitime_coefficient_volume_identities=coefficients,
                nonzero_local_coefficients=nonzero,
                coefficient_orders_per_gap=[0, 1, 2, 3],
                adjacent_heat_correlation=adjacent,separated_heat_correlation=separated,
                exact_heat_time='4 log 2, k=1',
                stationary_gram=stationary_gram,reflection_gram=reflected_gram,
                reversed_heat_sign_refused=True,
                positivity_and_markov_proof='YM58-T3 positive transitions and exact semigroup law',
                finite_Taylor_polynomials_claimed_positive=False)


def refinement_controls():
    # Extract the fourth moment of the actual global-label source, retaining
    # its N-dependent weights. On one linear coordinate D^2=-1/4.
    moments = []
    for N in (3, 7, 23):
        for k in (F(1, 10), F(1), F(3)):
            W, active = N+k*(N-1), 1+2*k
            labels = [(F(1), i == 1) for i in range(N)]
            labels += [(k, i in (0, 1)) for i in range(N-1)]
            second = sum(weight/W*W*F(-1, 4)/2 for weight, acts in labels if acts)
            fourth = sum(weight/W*W**2*F(1, 16)/24 for weight, acts in labels if acts)
            assert second == -active/8
            assert fourth == W*active/384
            moments.append(dict(width=N,k=k,second_coefficient=second,
                                fourth_coefficient=fourth))
    budgets = []
    for N in (4, 16, 64, 256):
        value = refinement_budget((0, 1), 1, N, N*N, 2, 1)
        budgets.append(dict(width=N,steps=N*N,degree=2,upper=value))
    assert all(budgets[i+1]['upper'] < budgets[i]['upper'] for i in range(len(budgets)-1))
    assert budgets[-1]['upper'] < F(1, 500)
    # e^(3/8) > 1+3/8+(3/8)^2/2 > 10/7 implies e^(-3/8)<7/10.
    assert 1+F(3, 8)+F(3, 8)**2/2 > F(10, 7)
    lower = fast_volume_multiplier_lower(100)
    assert lower == F(19399, 19999) and lower > F(96, 100) > F(7, 10)
    counterexamples = [dict(steps=n,width=n*n,multiplier_lower=fast_volume_multiplier_lower(n))
                       for n in (10, 100, 1000)]
    return dict(global_label_moment_checks=moments,joint_gate_budgets=budgets,
                sufficient_gate='W_N/n -> 0 on every positive fixed time gap',
                unrestricted_joint_limit=False,
                fast_volume_counterexamples=counterexamples,
                exact_heat_multiplier_upper=F(7, 10),
                volume_stabilization_at_finite_global_step=False)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = prior.digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = prior.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = prior.digest(SOURCES)
    return checked


def run():
    prior.previous.runtime_check()
    return prior.previous.encode(dict(
        certificate_type='YM58_CONSISTENT_CORRELATED_INFINITE_VOLUME_HISTORIES',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        boundary=boundary_controls(),locality=locality_controls(),energy=energy_controls(),
        histories=history_controls(),refinement=refinement_controls(),
        claim_status='SELECTED_CORRELATED_KINETIC_INFINITE_VOLUME_HISTORY_AND_L2_GAP__4D_CLAY_OPEN',
        proof_scope='Written all-content proof; exact finite controls are not formal verification',
        model=dict(infinite_volume_histories=True,row_markov_for_this_model=True,
                   stationary_L2_sector='constants',scalar_interaction_imported=False,
                   physical_clock_calibrated=False,spatial_continuum_constructed=False),
        open='Scalar interaction; physical gauge action, constraints and clock; spatial continuum; 4D Yang-Mills'))


def check(cert, result=RESULT, pin=PIN):
    sha = prior.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM58 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = prior.previous.native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM58 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],unit_edge_gap_lower='1/630',
                          marginal_identities=cert['locality']['open_volume_generator_identities'],
                          history_identities=cert['histories']['multitime_coefficient_volume_identities'])))


if __name__ == '__main__':
    main()
