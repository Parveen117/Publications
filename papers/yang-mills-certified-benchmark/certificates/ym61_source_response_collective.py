"""YM61: exact source-cycle response and collective spectral controls.

Python 3.12 only. Reuses the native polynomial, moment and matrix engines.
Fixed-order evidence does not assert tower convergence or a lower mass gap.
Default/--check are read-only. --write changes only this result and hash.
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
import ym60_interacting_histories_response as previous

source = previous.previous
native, heat = previous.native, previous.heat
a, p, y = previous.a, previous.p, previous.y
exact, pairing = previous.exact, previous.pairing
RESULT = HERE / 'YM61_RESULT.json'
PIN = HERE / 'EXPECTED_YM61.sha256'
SOURCES = HERE / 'YM61_SOURCE_PINS.json'


def cycle_edges(distance):
    native.width(distance)
    return tuple((i, i+1) for i in range(distance))+((0, distance),)


def incidences(edges, mask, width):
    return tuple(sum(1 for e, pair in enumerate(edges) if mask >> e & 1 and site in pair)
                 for site in range(width))


@lru_cache(maxsize=None)
def cycle_response(distance, k=F(1)):
    """Coefficient of one labelled cycle, with exact completed-site elimination."""
    edges, k = cycle_edges(distance), source.positive_rate(k)
    n, full = distance+1, (1 << len(edges))-1
    values = {0: source.constant(1, n)}
    potentials = tuple(source.overlap(i, j, n) for i, j in edges)
    witnesses = []
    for size in range(1, len(edges)+1):
        for mask in range(1, full+1):
            if mask.bit_count() != size:
                continue
            rhs = p.add(*(p.mul(potentials[e], values[mask ^ (1 << e)])
                          for e in range(len(edges)) if mask >> e & 1))
            if mask == full:
                return dict(distance=distance,leading= native.phi(rhs, n),witnesses=witnesses)
            degrees = incidences(edges, mask, n)
            active = tuple(i for i, d in enumerate(degrees) if d == 1)
            assert active and all(d <= 2 for d in degrees)
            for site, degree in enumerate(degrees):
                if degree == 2:
                    rhs = source.site_average(rhs, site, n)
            assert native.phi(rhs, n) == 0
            for monomial in rhs:
                assert tuple(sum(monomial[4*i:4*i+4]) for i in range(n)) == tuple(int(i in active) for i in range(n))
            inverse = source.inverse_source(rhs, tuple(range(n)), k)
            values[mask] = inverse['inverse']
            # Independent reduced multilinear kinetic expression (15).
            reduced = p.scale(values[mask], F(len(active))*(1+2*k)/8)
            mixed = [p.scale(native.derivative(native.derivative(values[mask], i, 2), i+1, 2), k)
                     for i in active if i+1 in active]
            difference = p.add(reduced, *mixed, p.scale(rhs, -1))
            assert pairing(difference, difference, n) == 0
            witnesses.append(dict(mask=mask,active=active,completed=tuple(i for i, d in enumerate(degrees) if d == 2),
                                  dimension=inverse['dimension'],source_terms=len(rhs),
                                  centered=True,inverse_residual_norm2=inverse['residual_norm2']))
    raise AssertionError('Cycle recursion did not reach the full set')


def separated_rate(k):
    return (1+2*source.positive_rate(k))/4


def triangular_response(k):
    omega, nu = source.omega(k), separated_rate(k)
    return 1/(4*omega*nu)+1/(8*omega*omega)


def band(width, distance):
    native.width(width)
    native.width(distance)
    return [[F(abs(i-j) == distance) for j in range(width)] for i in range(width)]


def matrix_sum(*terms):
    if not terms:
        raise ValueError('At least one matrix required')
    out = a.scale(terms[0], 0)
    for term in terms:
        out = a.add(out, term)
    return out


def matrix_series_product(left, right, order):
    native.width(order, 0)
    return [matrix_sum(*(a.multiply(left[i], right[n-i])
                         for i in range(n+1) if i < len(left) and n-i < len(right)))
            for n in range(order+1)]


def response_matrices(correlations, k=F(1)):
    correlations = tuple(map(exact, correlations))
    if not correlations or correlations[0] != 1:
        raise ValueError('Normalized covariance requires c_0=1')
    k, L = source.positive_rate(k), len(correlations)
    G = [[correlations[abs(i-j)] for j in range(L)] for i in range(L)]
    if not a.psd(G):
        raise ValueError('Correlations must give a positive semidefinite covariance')
    M = a.scale(a.identity(L), (1+2*k)/8)
    if L > 1:
        M = a.add(M, a.scale(band(L, 1), -k*correlations[1]/8))
    assert a.psd(M)
    return G, M


def quadratic(vector, matrix):
    vector = tuple(map(exact, vector))
    if len(vector) != len(matrix):
        raise ValueError('Vector and matrix dimensions differ')
    return sum((x*matrix[i][j]*y for i, x in enumerate(vector)
                for j, y in enumerate(vector)), F(0))


def probe_rate(vector, correlations, k=F(1)):
    G, M = response_matrices(correlations, k)
    variance = quadratic(vector, G)
    if variance <= 0:
        raise ValueError('A non-null readout is required')
    return quadratic(vector, M)/variance


def scalar_quotient(numerator, denominator):
    if not denominator or not denominator[0] or len(numerator) != len(denominator):
        raise ValueError('Equal finite series with nonzero constant denominator required')
    numerator, denominator = tuple(map(exact, numerator)), tuple(map(exact, denominator))
    result = []
    for n in range(len(numerator)):
        value = numerator[n]-sum((denominator[j]*result[n-j] for j in range(1, n+1)), F(0))
        result.append(value/denominator[0])
    return tuple(result)


def rational_probe_series(vector):
    vector = tuple(map(exact, vector))
    if not vector or not any(vector):
        raise ValueError('Nonzero finite readout required')
    L = len(vector)
    norm = sum((x*x for x in vector), F(0))
    ratios = [2*sum((vector[i]*vector[i+r] for i in range(L-r)), F(0))/norm
              for r in range(1, 4)]
    c1, c2, c3 = F(1, 2), cycle_response(2)['leading'], cycle_response(3)['leading']
    b1 = F(previous.response_controls()['bulk_correlation_cubic'])
    denominator = (F(1), ratios[0]*c1, ratios[1]*c2, ratios[0]*b1+ratios[2]*c3)
    numerator = (F(3, 8), -ratios[0]*c1/8, F(0), -ratios[0]*b1/8)
    return scalar_quotient(numerator, denominator)


@lru_cache(maxsize=None)
def matrix_tower(width):
    native.width(width, 2)
    I, J1, J2, J3 = a.identity(width), band(width, 1), band(width, 2), band(width, 3)
    zero = a.scale(I, 0)
    c2, c3 = cycle_response(2)['leading'], cycle_response(3)['leading']
    b1 = F(previous.response_controls()['bulk_correlation_cubic'])
    G1, G2 = a.scale(J1, F(1, 2)), a.scale(J2, c2)
    G3 = matrix_sum(a.scale(J1, b1), a.scale(J3, c3))
    G = [I, G1, G2, G3]
    G1square = a.multiply(G1, G1)
    G1cube = a.multiply(G1square, G1)
    S = [I, a.scale(G1, F(-1, 2)),
         matrix_sum(a.scale(G2, F(-1, 2)), a.scale(G1square, F(3, 8))),
         matrix_sum(a.scale(G3, F(-1, 2)),
                    a.scale(matrix_sum(a.multiply(G1, G2), a.multiply(G2, G1)), F(3, 8)),
                    a.scale(G1cube, F(-5, 16)))]
    normalized = matrix_series_product(matrix_series_product(S, G, 3), S, 3)
    assert normalized == [I, zero, zero, zero]
    M = [a.scale(I, F(3, 8)), a.scale(J1, F(-1, 16)), zero, a.scale(J1, -b1/8)]
    B = matrix_series_product(matrix_series_product(S, M, 3), S, 3)
    J1square = a.multiply(J1, J1)
    expected = [a.scale(I, F(3, 8)), a.scale(J1, F(-1, 4)),
                matrix_sum(a.scale(J1square, F(1, 8)), a.scale(J2, F(-11, 64))),
                matrix_sum(a.scale(J1, F(-1123, 23688)), a.scale(J3, F(-49, 288)),
                           a.scale(matrix_sum(a.multiply(J1, J2), a.multiply(J2, J1)), F(77, 768)),
                           a.scale(a.multiply(J1square, J1), F(-1, 16)))]
    assert B == expected and all(x == a.transpose(x) for x in B)
    return dict(normalization=normalized,response=B)


@lru_cache(maxsize=None)
def source_controls():
    cycles = [cycle_response(r) for r in range(1, 5)]
    assert [c['leading'] for c in cycles] == [F(1, 2), F(11, 24), F(49, 108), F(4853, 10368)]
    assert sum(len(c['witnesses']) for c in cycles) == 52
    assert [max(w['dimension'] for w in c['witnesses']) for c in cycles] == [1, 1, 2, 2]
    triangular = []
    for k in (F(1, 3), F(1), F(2)):
        result = cycle_response(2, k)['leading']
        assert result == triangular_response(k)
        triangular.append(dict(k=k,value=result))
    # A full degree-two source inverse independently checks safe projection.
    projections = []
    for k in (F(1, 3), F(1), F(2)):
        raw = source.fixture(1, k)
        reduced = source.site_average(raw['inverse'], 1, 3)
        expected = p.scale(source.overlap(0, 2, 3), 1/(4*separated_rate(k)))
        difference = p.add(reduced, p.scale(expected, -1))
        assert pairing(difference, difference, 3) == 0
        projections.append(dict(k=k,unprojected_dimension=raw['dimension'],projection_inverse_residual_norm2=F(0)))
    premature = source.site_average(source.overlap(0, 1, 2), 0, 2)
    assert not premature
    rates = []
    for k in (F(1, 3), F(1), F(2)):
        for r in range(1, 5):
            n = r+1
            w = source.overlap(0, r, n)
            rate = source.omega(k) if r == 1 else separated_rate(k)
            assert source.local_operator(w, tuple(range(n)), k) == p.scale(w, rate)
            assert pairing(w, w, n) == F(1, 4)
            rates.append(dict(k=k,distance=r,rate=rate,susceptibility=1/(2*rate)))
    return dict(cycles=cycles,triangular_formula_cells=triangular,
                inverse_projection_checks=projections,source_susceptibilities=rates,
                premature_averaging_erases_source=True,
                wrong_adjacent_rate_gives=F(3, 8),wrong_factorized_correlation_gives=F(1, 4))


def parity_controls():
    counts = []
    for r in range(1, 8):
        edges, n = cycle_edges(r), r+1
        even = []
        for mask in range(1 << len(edges)):
            if all(d % 2 == 0 for d in incidences(edges, mask, n)):
                even.append(mask)
        assert even == [0, (1 << len(edges))-1]
        counts.append(dict(distance=r,even_subsets=len(even)))
    # Explicit chain parity at orders below the source separation.
    cases = 0
    import itertools
    for r in range(1, 5):
        for powers in itertools.product(range(3), repeat=r):
            degree = [0]*(r+1)
            degree[0] = degree[-1] = 1  # closing source
            for i, power in enumerate(powers):
                degree[i] += power
                degree[i+1] += power
            allowed = all(d % 2 == 0 for d in degree)
            assert allowed == all(power % 2 == 1 for power in powers)
            if allowed:
                assert sum(powers) >= r and (sum(powers)-r) % 2 == 0
            cases += 1
        assert previous.secant_error_order(2*r+4, r+2) == r+2
    return dict(cycle_parity=counts,path_incidence_cases=cases,
                arbitrary_separation_first_possible_order='r',uniform_in_separation_error=False)


def collective_controls():
    first = []
    for k in (F(1, 7), F(1, 3), F(1), F(3)):
        rho, susceptibility = (1+2*k)/8, 1/(2*source.omega(k))
        assert susceptibility*(rho+k/8) == F(1, 4)
        first.append(dict(k=k,free_rate=rho,source_response=susceptibility,first_matrix_coefficient=F(-1, 4)))
    for L in (2, 3, 4, 5, 6):
        matrix_tower(L)
        J, J2 = band(L, 1), band(L, 2)
        diagonal = [[F((1 if i in (0, L-1) else 2) if i == j else 0) for j in range(L)] for i in range(L)]
        assert a.multiply(J, J) == a.add(J2, diagonal)
    J1, J2 = band(4, 1), band(4, 2)
    AB, BA = a.multiply(J1, J2), a.multiply(J2, J1)
    assert AB != BA and AB != a.transpose(AB)
    assert a.add(AB, BA) == a.transpose(a.add(AB, BA))
    # L=3 leading sine vector u=(1,sqrt(2),1): polynomial identities use
    # J^2 u=2u and normalized endpoint mass 1/2, avoiding rounded roots.
    L3_second = F(-3, 64)*2+F(11, 64)*F(3, 2)
    assert L3_second == F(21, 128)
    vectors = ((1, 1), (1, 1, 1), (1, 2, 2, 1), (1, 2, 3, 3, 2, 1))
    probes = [dict(weights=b,series=rational_probe_series(b)) for b in vectors]
    assert probes[0]['series'] == tuple(map(F, previous.response_controls()['trial_rate_series']))
    assert probes[2]['series'] == (F(3, 8), F(-2, 5), F(73, 400), F(-1293529, 11844000))
    # A bona fide Toeplitz covariance fixture, not an assumed true vacuum.
    c = tuple(F(1, 3)**r for r in range(5))
    G, M = response_matrices(c)
    b = (1, 2, 1, 0)
    assert probe_rate(b, c[:4]) == probe_rate(b+(0,), c)
    for s in (-1, 1):
        staggered_c = tuple(v*s**r for r, v in enumerate(c))
        staggered_b = tuple(v*s**i for i, v in enumerate(b+(0,)))
        assert probe_rate(staggered_b, staggered_c) == probe_rate(b+(0,), c)
    return dict(first_response_cells=first,congruence_widths=[2, 3, 4, 5, 6],
                boundary_matrices_do_not_commute=True,L3_optimal_quadratic=L3_second,
                rational_probes=probes,zero_padding_preserves_rayleigh_quotient=True,
                coupling_sign_staggering_preserves_quotient=True)


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
        certificate_type='YM61_NATIVE_SOURCE_CYCLE_RESPONSE_AND_COLLECTIVE_SPECTRAL_CONSTRAINT',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        source=source_controls(),parity=parity_controls(),collective=collective_controls(),
        claim_status='UNIFORM_FIXED_SEPARATION_RESPONSE_AND_COLLECTIVE_UPPER_CONSTRAINT__LOWER_GAP_OPEN',
        proof_scope='Written source/limit proof; exact finite controls do not establish infinite-series convergence',
        model=dict(actual_scalar_operator_retained=True,independent_sources_removed_after_response=True,
                   fixed_distance_uniform_asymptotics=True,collective_first_response_independent_of_k=True,
                   cubic_matrix_response_noncommutative=True,uniform_in_distance_series=False,
                   exact_mass_eigenvalue=False,uniform_interacting_lower_gap=False,physical_4d_yang_mills=False),
        open='Uniform source-tower estimates or full-vacuum Poincare control; physical gauge and continuum'))


def check(cert, result=RESULT, pin=PIN):
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM61 fresh certificate/pin mismatch')
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
    print('YM61 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],leading_spatial_responses=[x['leading'] for x in cert['source']['cycles']],
                          first_collective_matrix_coefficient='-1/4',uniform_interacting_lower_gap=False)))


if __name__ == '__main__':
    main()
