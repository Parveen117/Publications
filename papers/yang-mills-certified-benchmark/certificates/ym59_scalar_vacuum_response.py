"""YM59: exact scalar vacuum response and local Poisson controls.

Python 3.12 only. Reuses the native derivative, moment and linear solver.
Finite checks support the written uniform finite-order proof, not a
nonzero-coupling gap claim. Default/--check are read-only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym58_infinite_volume_histories as prior

native = prior.prior
a, p, y = prior.a, prior.p, prior.y
exact = prior.exact
RESULT = HERE / 'YM59_RESULT.json'
PIN = HERE / 'EXPECTED_YM59.sha256'
SOURCES = HERE / 'YM59_SOURCE_PINS.json'


def positive_rate(k):
    k = prior.rate(k)
    if not k:
        raise ValueError('Positive kinetic neighbour rate required')
    return k


def constant(value, n):
    value = exact(value)
    native.width(n, 0)
    return {(0,)*(4*n): value} if value else {}


def pairing(f, g, n):
    return native.phi(p.mul(f, g), n)


def center(f, n):
    return p.add(f, constant(-native.phi(f, n), n))


def overlap(i, j, n):
    native.width(n)
    if i == j:
        raise ValueError('Distinct overlap endpoints required')
    return p.add(*(p.mul(native.lift(f, i, n), native.lift(f, j, n)) for f in p.COORD))


def degrees_for(n, boundary='open'):
    native.width(n, 2)
    if boundary == 'bulk':
        return (2,)*n
    if boundary == 'open':
        return (1,)+(2,)*(n-2)+(1,)
    raise ValueError('Unknown boundary convention')


def omega(k, left=2, right=2):
    k = positive_rate(k)
    if any(isinstance(d, bool) or not isinstance(d, int) or d not in (1, 2) for d in (left, right)):
        raise ValueError('A chain bond endpoint has degree one or two')
    return (2+k*(left+right+2))/8


def local_operator(f, support, k, degrees=None):
    support, k = prior.sites(support), positive_rate(k)
    degrees = (2,)*len(support) if degrees is None else tuple(degrees)
    if len(degrees) != len(support) or any(isinstance(d, bool) or not isinstance(d, int) or d not in (1, 2) for d in degrees):
        raise ValueError('One physical degree per supported site required')
    for i, site in enumerate(support):
        internal = int(site-1 in support)+int(site+1 in support)
        if degrees[i] < internal:
            raise ValueError('Physical degree cannot delete an internal bond')
    full = prior.marginal(f, support, k)
    corrections = [p.scale(native.derivative(native.derivative(f, i, 2), i, 2), k*(2-d)/2)
                   for i, d in enumerate(degrees) if d != 2]
    return p.add(full, *corrections)


def site_average(f, site, n):
    out = {}
    for m, c in f.items():
        coefficient = c*y.phi({m[4*site:4*site+4]: F(1)})
        term = m[:4*site]+(0,)*4+m[4*site+4:]
        out[term] = out.get(term, F(0))+coefficient
    return p.clean(out)


def gamma(f, g, support, k, degrees=None):
    H = lambda h: local_operator(h, support, k, degrees)
    return p.scale(p.add(p.mul(f, H(g)), p.mul(g, H(f)), p.scale(H(p.mul(f, g)), -1)), F(1, 2))


def inverse_source(f, support, k, degrees=None, max_dimension=32):
    """Exact Krylov inverse in the native pairing; refuses noncentered sources."""
    support, k = prior.sites(support), positive_rate(k)
    native.width(max_dimension)
    n = len(support)
    prior.polynomial(f, support)
    if native.phi(f, n):
        raise ValueError('Poisson source must be centered')
    H = lambda g: local_operator(g, support, k, degrees)
    basis, norms = [], []
    current = f
    for _ in range(max_dimension+1):
        for b, norm in zip(basis, norms):
            current = p.add(current, p.scale(b, -pairing(b, current, n)/norm))
        norm = pairing(current, current, n)
        if not norm:
            break
        if len(basis) == max_dimension:
            raise ValueError('Finite fixture exceeds declared Krylov budget')
        basis.append(current)
        norms.append(norm)
        current = H(current)
    if not basis:
        return dict(inverse={},dimension=0,quadratic=F(0),residual_norm2=F(0))
    images = [H(b) for b in basis]
    form = [[pairing(b, image, n) for image in images] for b in basis]
    assert form == a.transpose(form) and a.psd(form)
    rhs = [pairing(b, f, n) for b in basis]
    coefficients = a.y37.solve(form, rhs)
    inverse = p.add(*(p.scale(b, c) for b, c in zip(basis, coefficients)))
    residual = p.add(H(inverse), p.scale(f, -1))
    assert pairing(residual, residual, n) == 0
    assert native.phi(inverse, n) == 0
    return dict(inverse=inverse,dimension=len(basis),quadratic=pairing(f, inverse, n),
                residual_norm2=F(0))


@lru_cache(maxsize=None)
def fixture(distance, k=F(1), degrees=None):
    native.width(distance, 0)
    if distance not in (0, 1, 2, 3):
        raise ValueError('Fixture distance must be zero through three')
    n = min(distance+2, 4)
    support = tuple(range(n)) if distance < 3 else (0, 1, 3, 4)
    v = overlap(0, 1, n)
    if distance == 0:
        radius = p.mul(native.lift(y.RADIUS, 0, n), native.lift(y.RADIUS, 1, n))
        source = p.add(p.mul(v, v), p.scale(radius, F(-1, 4)))
    else:
        second = distance if distance < 3 else 2
        source = p.mul(v, overlap(second, second+1, n))
    assert pairing(source, source, n) == F(1, 16)
    result = inverse_source(source, support, k, degrees)
    return dict(result,source=source,support=support)


def energy_coefficients(n, k=F(1), boundary='open'):
    native.width(n, 2)
    k = positive_rate(k)
    degrees = degrees_for(n, boundary)
    rates = [omega(k, degrees[i], degrees[i+1]) for i in range(n-1)]
    c2 = -sum((F(1)/(4*w) for w in rates), F(0))
    c4 = F(0)
    for i, w in enumerate(rates):
        value = fixture(0, k, degrees[i:i+2])['quadratic']
        c4 += F(1)/(16*w**3)-value/w**2
    for distance in (1, 2):
        for i in range(n-1-distance):
            u, v = rates[i], rates[i+distance]
            value = fixture(distance, k, degrees[i:i+distance+2])['quadratic']
            c4 += (u+v)/(16*u*u*v*v)-(1/u+1/v)**2*value
    return c2, c4


def weak_energy_budget(k, theta):
    k, theta = positive_rate(k), exact(theta)
    gamma0 = native.lower_bound(2, k)
    margin = gamma0-2*abs(theta)
    if margin <= 0:
        raise ValueError('This energy budget needs abs(theta) < gamma_k/2')
    return dict(gamma=gamma0,margin=margin,
                energy_lower_per_bond=-theta*theta/(4*margin),
                excitation_upper_per_bond=theta*theta/(margin*margin))


def finite_gap_bound(n, k, theta):
    native.width(n, 2)
    return native.lower_bound(2, positive_rate(k))-2*abs(exact(theta))*(n-1)


def coercivity_controls():
    matrices = projections = 0
    for degrees in ((1, 0), (1, 1), (2, 0), (2, 1), (1, 0, 1)):
        n = len(degrees)
        mons, polys, gram, _ = native.block(degrees)
        number_form = a.scale(gram, 0)
        for site in range(n):
            images = [p.add(f, p.scale(site_average(f, site, n), -1)) for f in polys]
            number_form = a.add(number_form, [[pairing(f, g, n) for g in images] for f in images])
        for k in (F(1, 10), F(1), F(3)):
            columns = [native.operator(f, n, k) for f in polys]
            H = [[col.get(m, F(0)) for col in columns] for m in mons]
            assert a.psd(a.add(a.multiply(gram, H), a.scale(number_form, -native.lower_bound(n, k))))
            matrices += 1
            for f in polys[:4]:
                for site in range(n):
                    assert site_average(native.operator(f, n, k), site, n) == native.operator(site_average(f, site, n), n, k)
                    projections += 1
    return dict(excited_site_form_inequalities=matrices,projection_commutations=projections,
                full_proof='YM59-T1 inherits the summed local echo estimate before tensorization')


def bond_controls():
    eigenvalues = parity = 0
    rows = []
    for n in (2, 3, 4, 7):
        degrees = degrees_for(n)
        for k in (F(1, 10), F(1), F(3)):
            rates = []
            for i in range(n-1):
                v = overlap(i, i+1, n)
                w = omega(k, degrees[i], degrees[i+1])
                assert native.operator(v, n, k) == p.scale(v, w)
                assert native.phi(v, n) == 0 and pairing(v, v, n) == F(1, 4)
                assert not site_average(v, i, n) and not site_average(v, i+1, n)
                rates.append(w)
                eigenvalues += 1
            assert sum(rates) == (n-1)*omega(k)-k/4
            rows.append(dict(width=n,k=k,rates=rates,c2=-sum((1/(4*w) for w in rates), F(0))))
        def flip(f):
            return {m:c*(-1)**sum(sum(m[4*i:4*i+4]) for i in range(0, n, 2)) for m, c in f.items()}
        V = p.add(*(overlap(i, i+1, n) for i in range(n-1)))
        assert flip(V) == p.scale(V, -1)
        for f in (V, p.mul(V, V), native.lift(p.COORD[0], 0, n)):
            assert flip(native.operator(f, n, 1)) == native.operator(flip(f), n, 1)
            parity += 1
    return dict(local_bond_eigenvalue_checks=eigenvalues,rows=rows,
                antipodal_intertwinings=parity,bulk_quadratic_at_unit_rate=F(-1, 4),
                bulk_susceptibility_at_unit_rate=F(1, 2))


def poisson_controls():
    cells = []
    expected = ((3, F(1, 30)), (4, F(31, 705)), (2, F(2, 63)), (1, F(1, 32)))
    for distance, (dimension, value) in enumerate(expected):
        data = fixture(distance)
        assert data['dimension'] == dimension and data['quadratic'] == value
        cells.append(dict(distance=distance,dimension=dimension,quadratic=value,
                          inverse_terms=len(data['inverse']),residual_norm2=data['residual_norm2']))
    # Independent logarithmic recursion on one bulk bond. It checks the
    # product-rule sign and obtains c4 without the Rayleigh formula (15).
    support, n = (0, 1), 2
    f1 = overlap(0, 1, n)
    g2 = gamma(f1, f1, support, 1)
    c2 = -native.phi(g2, n)
    f2 = inverse_source(center(g2, n), support, 1)['inverse']
    g3 = p.scale(gamma(f1, f2, support, 1), 2)
    assert native.phi(g3, n) == 0
    c4 = -4*pairing(f1, gamma(f1, f2, support, 1), n)-native.phi(gamma(f2, f2, support, 1), n)
    assert c2 == F(-1, 4) and c4 == F(7, 240)
    # Verify Gamma independently by the original incident-edge squares.
    direct = p.scale(p.add(*(p.scale(p.mul(prior.apply_direction(f1, word),
                                          prior.apply_direction(f2, word)), weight)
                             for weight, word in prior.directions(support, 1))), F(1, 2))
    assert direct == gamma(f1, f2, support, 1)
    # Order-two logarithmic residual has exactly its theta^2 constant,
    # theta^3=-2 Gamma(f1,f2), theta^4=-Gamma(f2,f2) on Q.
    for theta in (F(-1, 3), F(1, 10), F(2, 5)):
        Fp = p.add(p.scale(f1, theta), p.scale(f2, theta*theta))
        residual = p.add(local_operator(Fp, support, 1), p.scale(gamma(Fp, Fp, support, 1), -1),
                         p.scale(f1, -theta))
        predicted = p.add(constant(c2*theta**2, n), p.scale(g3, -theta**3),
                          p.scale(gamma(f2, f2, support, 1), -theta**4))
        difference = p.add(residual, p.scale(predicted, -1))
        assert pairing(difference, difference, n) == 0
    return dict(exact_local_inverses=cells,logarithmic_c2=c2,logarithmic_c4_single_bond=c4,
                logarithmic_residual_identities=3,independent_square_product_rule=True,
                all_fixed_orders_proof='YM59-T3 local Poisson recursion and scalar residual bound',
                numerical_order_five_remainder_constant_evaluated=False,
                infinite_series_convergence_asserted=False)


def fourth_order_controls():
    a0, a1, a2 = [fixture(d)['quadratic'] for d in range(3)]
    connected = (F(1, 16)-a0, F(1, 8)-4*a1, F(1, 8)-4*a2)
    expected = (F(7, 240), F(-287, 5640), F(-1, 504))
    assert connected == expected
    bulk = sum(connected)
    assert bulk == F(-1123, 47376)
    marginal_rows, open_rows = [], []
    for n in (4, 5, 6, 9):
        c2, c4 = energy_coefficients(n, boundary='bulk')
        assert c2 == F(-(n-1), 4)
        assert c4 == (n-1)*bulk+F(6497, 118440)
        marginal_rows.append(dict(width=n,c2=c2,c4=c4))
    for n in (2, 3, 4, 5, 6, 7, 9):
        c2, c4 = energy_coefficients(n)
        open_rows.append(dict(width=n,c2=c2,c4=c4))
    by_n = {row['width']:row for row in open_rows}
    intercept = by_n[5]['c4']-4*bulk
    for n in (6, 7, 9):
        assert by_n[n]['c4'] == (n-1)*bulk+intercept
    # Arbitrarily separated bonds cancel even at unequal endpoint rates.
    cancellations = 0
    for u, v in itertools.product((F(3, 4), F(7, 8), F(1)), repeat=2):
        source_value = F(1)/(16*(u+v))
        assert (1/u+1/v)**2*source_value == (u+v)/(16*u*u*v*v)
        cancellations += 1
    # Independent assembled full marginal inverse on three scalar bonds.
    n = 4
    bonds = [overlap(i, i+1, n) for i in range(n-1)]
    G, RG = {}, {}
    for i in range(n-1):
        data = fixture(0)
        shifted = tuple(range(i, i+2))
        G = p.add(G, prior.embed(data['source'], shifted, tuple(range(n))))
        RG = p.add(RG, prior.embed(data['inverse'], shifted, tuple(range(n))))
    for i, j in itertools.combinations(range(n-1), 2):
        data = fixture(j-i)
        shifted = tuple(range(i, j+2))
        G = p.add(G, p.scale(prior.embed(data['source'], shifted, tuple(range(n))), 2))
        RG = p.add(RG, p.scale(prior.embed(data['inverse'], shifted, tuple(range(n))), 2))
    residual = p.add(local_operator(RG, tuple(range(n)), 1), p.scale(G, -1))
    assert pairing(residual, residual, n) == 0
    V = p.add(*bonds)
    literal_source = p.add(p.mul(V, V), constant(F(-(n-1), 4), n))
    source_defect = p.add(literal_source, p.scale(G, -1))
    assert pairing(source_defect, source_defect, n) == 0
    assert pairing(V, V, n) == F(n-1, 4)
    c4_direct = -pairing(G, RG, n)+F((n-1)**2, 16)
    assert c4_direct == energy_coefficients(n, boundary='bulk')[1]
    wrong_disconnected = connected[0]+connected[1]
    assert wrong_disconnected != bulk
    unnormalized = [energy_coefficients(M+1, boundary='bulk')[1]-F(M*M, 16) for M in (3, 4, 5)]
    wrong_second_difference = unnormalized[2]-2*unnormalized[1]+unnormalized[0]
    assert wrong_second_difference == F(-1, 8)
    return dict(connected_terms=connected,bulk_quartic=bulk,
                marginal_rows=marginal_rows,open_rows=open_rows,open_boundary_intercept=intercept,
                unequal_rate_disconnected_cancellations=cancellations,
                assembled_full_source_inverse=True,assembled_c4=c4_direct,
                deleting_kinetic_bridge_gives=wrong_disconnected,
                normalization_removed_second_difference=wrong_second_difference,
                missing_normalization_has_spurious_quadratic_volume_growth=True)


def energy_bound_controls():
    incidence = 0
    for n in range(2, 9):
        for mask in itertools.product((0, 1), repeat=n):
            count = sum(mask)
            pairs = sum(mask[i]*mask[i+1] for i in range(n-1))
            active = sum(max(mask[i], mask[i+1]) for i in range(n-1))
            assert pairs <= count and active <= 2*count
            incidence += 1
    rows = []
    for k in (F(1, 10), F(1), F(3)):
        g = native.lower_bound(2, k)
        for sign in (-1, 1):
            theta = sign*g/8
            budget = weak_energy_budget(k, theta)
            assert budget['margin'] == 3*g/4
            assert budget['excitation_upper_per_bond'] == F(1, 36)
            # 2x2 Young budget for vacuum-to-pair amplitude 1/2.
            delta = budget['margin']
            assert a.psd([[theta*theta/(4*delta), -theta/2], [-theta/2, delta]])
            rows.append(dict(k=k,theta=theta,**budget))
    # Exact dimer Rayleigh quotient and a rational strictly negative trial.
    v = overlap(0, 1, 2)
    dimer_trials = []
    for theta in (F(-1, 10), F(1, 100), F(1, 2)):
        c = theta  # omega(1)=1; enough for a rational negative trial.
        psi = p.add(constant(1, 2), p.scale(v, c))
        Hpsi = p.add(local_operator(psi, (0, 1), 1), p.scale(p.mul(v, psi), -theta))
        actual = pairing(psi, Hpsi, 2)/pairing(psi, psi, 2)
        expected = (c*c/4-theta*c/2)/(1+c*c/4)
        assert actual == expected < 0
        dimer_trials.append(dict(theta=theta,rayleigh=actual))
    g = native.lower_bound(2, 1)
    theta = g/8
    assert finite_gap_bound(2, 1, theta) > 0
    assert finite_gap_bound(5, 1, theta) == 0
    assert finite_gap_bound(6, 1, theta) < 0
    return dict(excitation_incidence_cases=incidence,weak_energy_cells=rows,
                exact_negative_dimer_trials=dimer_trials,shrinking_gap_window_refused=True,
                fixed_nonzero_coupling_uniform_gap_proved=False)


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
    return native.previous.encode(dict(
        certificate_type='YM59_SCALAR_INTERACTION_LOCAL_POISSON_TOWER_AND_BULK_ENERGY',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        coercivity=coercivity_controls(),bonds=bond_controls(),poisson=poisson_controls(),
        fourth_order=fourth_order_controls(),energy=energy_bound_controls(),
        claim_status='ACTUAL_SCALAR_INTERACTION_THERMODYNAMIC_ENERGY_AND_UNIFORM_FIXED_ORDER_RESPONSE__4D_CLAY_OPEN',
        proof_scope='Written all-width and fixed-order proof; exact finite controls are not formal verification',
        model=dict(original_scalar_interaction_retained=True,thermodynamic_energy_exists=True,
                   uniform_finite_order_energy_asymptotics=True,interacting_uniform_gap=False,
                   interacting_history_limit=False,infinite_series_convergence=False),
        open='Nonzero-coupling uniform gap and interacting histories; physical gauge selection and continuum; 4D Yang-Mills'))


def check(cert, result=RESULT, pin=PIN):
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM59 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM59 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],bulk_quadratic='-1/4',
                          bulk_quartic=cert['fourth_order']['bulk_quartic'],
                          local_inverse_dimensions=[r['dimension'] for r in cert['poisson']['exact_local_inverses']])))


if __name__ == '__main__':
    main()
