"""YM60: actual-vacuum history and response-to-spectrum finite controls.

Python 3.12 only. Exact finite evidence supports the written analytic proof.
No numerical ground vector or uniform interacting lower gap is claimed.
Default/--check are read-only; --write changes only this result and hash.
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
import ym59_scalar_vacuum_response as previous

native, heat = previous.native, previous.prior
a, p, y = previous.a, previous.p, previous.y
exact, pairing, gamma = previous.exact, previous.pairing, previous.gamma
RESULT = HERE / 'YM60_RESULT.json'
PIN = HERE / 'EXPECTED_YM60.sha256'
SOURCES = HERE / 'YM60_SOURCE_PINS.json'
UNITS = tuple(tuple(F(i == j) for i in range(4)) for j in range(4))


def signed_transform(poly, matrix, n):
    """Common signed permutation of the four coordinates; no new operator engine."""
    native.width(n)
    heat.polynomial(poly, tuple(range(n)))
    if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
        raise ValueError('Four-coordinate signed permutation required')
    entries = []
    for row in matrix:
        active = [(j, exact(c)) for j, c in enumerate(row) if c]
        if len(active) != 1 or active[0][1] not in (-1, 1):
            raise ValueError('Signed permutation required')
        entries.append(active[0])
    if len({j for j, _ in entries}) != 4:
        raise ValueError('Permutation must be invertible')
    out = {}
    for m, c in poly.items():
        exponent, sign = [0]*(4*n), F(1)
        for site in range(n):
            for coordinate, (target, factor) in enumerate(entries):
                power = m[4*site+coordinate]
                exponent[4*site+target] += power
                sign *= factor**power
        key = tuple(exponent)
        out[key] = out.get(key, F(0))+c*sign
    return p.clean(out)


def quaternion_average(poly, n):
    """Average over common left/right Q8; global minus removes odd degrees."""
    terms = []
    for left in UNITS:
        for right in UNITS:
            matrix = a.multiply(y.qmatrix(left), y.qmatrix(right, True))
            terms.append(signed_transform(poly, matrix, n))
            terms.append(signed_transform(poly, a.scale(matrix, -1), n))
    return p.scale(p.add(*terms), F(1, 32))


def trial_rate(k, correlation):
    k, x = previous.positive_rate(k), exact(correlation)
    if not 0 <= x <= 1:
        raise ValueError('Aligned correlation must lie in [0,1]')
    return (1+2*k-k*x)/(8*(1+x))


def finite_trial_rate(n, k, aligned_correlation):
    native.width(n, 2)
    k, x = previous.positive_rate(k), exact(aligned_correlation)
    trial_rate(k, x)  # validates the physical correlation
    return (1+2*k-k/F(n-1)-k*x)/(8*(1+x))


def secant_error_order(energy_remainder_order, step_power):
    """For bounded P'', secant errors are O(h)+O(theta^r/h)."""
    native.width(energy_remainder_order, 2)
    native.width(step_power)
    if step_power >= energy_remainder_order:
        raise ValueError('Secants must leave a vanishing remainder')
    return min(step_power, energy_remainder_order-step_power)


def increment_budget(poly, support, k, time):
    support, time = heat.sites(support), exact(time)
    heat.polynomial(poly, support)
    if time < 0:
        raise ValueError('Nonnegative time required')
    g = gamma(poly, poly, support, k)
    return 2*time*sum(map(abs, g.values()), F(0))


def spectral_summary(atoms):
    atoms = tuple((exact(rate), exact(weight)) for rate, weight in atoms)
    if not atoms or any(rate < 0 or weight <= 0 for rate, weight in atoms):
        raise ValueError('Positive spectral weights and nonnegative rates required')
    mass = sum((weight for _, weight in atoms), F(0))
    mean = sum((rate*weight for rate, weight in atoms), F(0))/mass
    moment = [[sum((weight*rate**(i+j) for rate, weight in atoms), F(0))
               for j in range(3)] for i in range(3)]
    assert a.psd(moment)
    return dict(mass=mass,mean=mean,threshold=min(rate for rate, _ in atoms),
                positive_moment_matrix=True)


@lru_cache(maxsize=None)
def bracket_controls():
    checks = 0
    for n in (2, 3, 4):
        for site in range(n):
            edge = min(site, n-2)
            epsilon = 1 if site == edge else -1
            D = lambda f: native.derivative(f, site, 1)
            B = lambda f: native.bond(f, edge)
            bracket = lambda f: p.add(D(B(f)), p.scale(B(D(f)), -1))
            for coordinate in p.COORD:
                f = native.lift(coordinate, site, n)
                assert bracket(f) == p.scale(native.derivative(f, site, 3), epsilon)
                nested = p.add(D(bracket(f)), p.scale(bracket(D(f)), -1))
                assert nested == p.scale(native.derivative(f, site, 2), -epsilon)
                checks += 1
    return dict(site_coordinate_bracket_checks=checks,all_physical_endpoints_covered=True)


@lru_cache(maxsize=None)
def symmetry_controls():
    intertwinings = bilinears = means = 0
    for n in (2, 3):
        V = p.add(*(previous.overlap(i, i+1, n) for i in range(n-1)))
        specimens = [native.lift(p.COORD[0], 0, n), V, p.mul(V, V)]
        for unit in UNITS:
            for right in (False, True):
                matrix = y.qmatrix(unit, right)
                T = lambda f: signed_transform(f, matrix, n)
                assert T(V) == V
                for f in specimens:
                    assert T(native.operator(f, n, 1)) == native.operator(T(f), n, 1)
                    intertwinings += 1
        for site in range(n):
            for f in p.COORD:
                assert not quaternion_average(native.lift(f, site, n), n)
                means += 1
        for i, j in ((0, 0), (0, n-1)):
            overlap = native.lift(y.RADIUS, i, n) if i == j else previous.overlap(i, j, n)
            for alpha in range(4):
                for beta in range(4):
                    product = p.mul(native.lift(p.COORD[alpha], i, n),
                                    native.lift(p.COORD[beta], j, n))
                    expected = p.scale(overlap, F(1, 4)) if alpha == beta else {}
                    assert quaternion_average(product, n) == expected
                    bilinears += 1
    # Right invariance alone leaves a non-scalar real commutant.
    J = y.qmatrix(UNITS[1])
    assert J != a.scale(a.identity(4), J[0][0])
    for unit in UNITS:
        R = y.qmatrix(unit, True)
        assert a.multiply(J, R) == a.multiply(R, J)
    return dict(kinetic_symmetry_intertwinings=intertwinings,
                bilinear_average_identities=bilinears,zero_mean_identities=means,
                right_symmetry_alone_insufficient=True)


@lru_cache(maxsize=None)
def cancellation_controls():
    identities = commutators = 0
    missing_ground_residual = None
    for n in (2, 3):
        support, degrees = tuple(range(n)), previous.degrees_for(n)
        V = p.add(*(previous.overlap(i, i+1, n) for i in range(n-1)))
        h = p.add(previous.constant(1, n), p.scale(V, F(1, 7)))
        theta = F(1, 5)
        H = lambda f: native.operator(f, n, 1)
        A = lambda f: p.add(H(f), p.scale(p.mul(V, f), -theta))
        specimens = (native.lift(p.COORD[0], 0, n),
                     native.lift(p.COORD[0], n-1, n), V)
        for f in specimens:
            for g in specimens:
                actual = pairing(p.mul(h, f), A(p.mul(h, g)), n)-pairing(p.mul(h, p.mul(f, g)), A(h), n)
                expected = native.phi(p.mul(p.mul(h, h), gamma(f, g, support, 1, degrees)), n)
                assert actual == expected
                identities += 1
                inner = lambda r: p.add(A(p.mul(g, r)), p.scale(p.mul(g, A(r)), -1))
                double = p.add(p.mul(f, inner(h)), p.scale(inner(p.mul(f, h)), -1))
                assert double == p.scale(p.mul(h, gamma(f, g, support, 1, degrees)), 2)
                commutators += 1
        Etrial = pairing(h, A(h), n)/pairing(h, h, n)
        residual = p.add(A(h), p.scale(h, -Etrial))
        norm2 = pairing(residual, residual, n)
        assert norm2 > 0
        missing_ground_residual = norm2
    return dict(weighted_cancellation_identities=identities,double_commutator_identities=commutators,
                chosen_polynomial_is_not_ground_vector_residual_norm2=missing_ground_residual)


@lru_cache(maxsize=None)
def form_controls():
    checks = 0
    for n in (2, 3):
        support, degrees = tuple(range(n)), previous.degrees_for(n)
        V = p.add(*(previous.overlap(i, i+1, n) for i in range(n-1)))
        for epsilon in (F(-1, 7), F(1, 5)):
            h = p.add(previous.constant(1, n), p.scale(V, epsilon))
            density = p.mul(h, h)
            norm = native.phi(density, n)
            mean = lambda f: native.phi(p.mul(density, f), n)/norm
            for k in (F(1, 3), F(1)):
                for i in range(n):
                    for j in range(i, n):
                        cij = F(1) if i == j else mean(previous.overlap(i, j, n))
                        for alpha, beta in ((0, 0), (0, 2), (1, 1)):
                            f, g = native.lift(p.COORD[alpha], i, n), native.lift(p.COORD[beta], j, n)
                            assert mean(p.mul(f, g)) == (cij/4 if alpha == beta else 0)
                            expected = ((1+k*degrees[i])/32 if i == j else
                                        -k*cij/32 if j == i+1 else F(0)) if alpha == beta else F(0)
                            assert mean(gamma(f, g, support, k, degrees)) == expected
                            checks += 1
                s = 1 if epsilon > 0 else -1
                cs = [mean(previous.overlap(i, i+1, n)) for i in range(n-1)]
                numerator = denominator = F(0)
                for i in range(n-1):
                    f = p.add(native.lift(p.COORD[0], i, n), p.scale(native.lift(p.COORD[0], i+1, n), s))
                    numerator += mean(gamma(f, f, support, k, degrees))
                    denominator += mean(p.mul(f, f))
                assert numerator/denominator == finite_trial_rate(n, k, s*sum(cs)/F(n-1))
    incidence = []
    for n in (2, 3, 4, 9, 100):
        degrees = previous.degrees_for(n)
        assert sum(degrees[i]+degrees[i+1] for i in range(n-1)) == 4*n-6
        incidence.append(n)
    return dict(weighted_form_entries=checks,physical_incidence_widths=incidence,
                weighted_sources='Positive invariant polynomial densities used only as identity controls',
                scalar_vacuum_not_numerically_approximated=True)


@lru_cache(maxsize=None)
def response_controls():
    # Recompute the three exact local sources; do not copy a fitted decimal.
    pieces = (F(1, 16)-previous.fixture(0)['quadratic'],
              F(1, 8)-4*previous.fixture(1)['quadratic'],
              F(1, 8)-4*previous.fixture(2)['quadratic'])
    energy4 = sum(pieces)
    c3 = -4*energy4
    assert c3 == F(1123, 11844)
    # Formal quotient in t=abs(theta), through degree three.
    corr = (F(0), F(1, 2), F(0), c3)
    denominator = [8*(1+corr[0])]+[8*c for c in corr[1:]]
    numerator = [3-corr[0]]+[-c for c in corr[1:]]
    result = []
    for order in range(4):
        value = numerator[order]-sum((denominator[j]*result[order-j] for j in range(1, order+1)), F(0))
        result.append(value/denominator[0])
    assert tuple(result) == (F(3, 8), F(-1, 4), F(1, 8), F(-5207, 47376))
    assert secant_error_order(10, 5) == 5
    assert secant_error_order(6, 3) == 3
    rows = []
    for k in (F(1, 10), F(1), F(3)):
        for x in (F(0), F(1, 7), F(1, 2), F(1)):
            rate = trial_rate(k, x)
            assert (1+k)/16 <= rate <= (1+2*k)/8
            if x:
                assert rate < trial_rate(k, 0)
            rows.append(dict(k=k,aligned_correlation=x,upper_rate=rate))
    return dict(bulk_energy_quartic=energy4,bulk_correlation_cubic=c3,
                trial_rate_series=result,exact_rate_cells=rows,
                derivative_remainder_requires_higher_energy_order=True,
                uniform_lower_gap_inferred=False)


def history_controls():
    # A convex mixture of reversible Markov histories is reflection positive,
    # but the observed time-zero process need not retain the Markov property.
    r, s = F(1, 2), F(1, 4)
    C = lambda n: (r**n+s**n)/2
    reflected = [[C(i+j) for j in range(4)] for i in range(4)]
    assert a.psd(reflected)
    defect = C(2)-C(1)**2
    assert defect == F(1, 64)
    atoms = [spectral_summary(((eps, F(1, 2)), (2-eps, F(1, 2))))
             for eps in (F(1, 2), F(1, 10), F(1, 1000))]
    atoms.append(spectral_summary(((F(0), F(1, 2)), (F(2), F(1, 2)))))
    assert all(row['mean'] == 1 for row in atoms)
    assert atoms[-1]['threshold'] == 0
    f = p.add(native.lift(p.COORD[0], 0, 2), native.lift(p.COORD[0], 1, 2))
    budget = increment_budget(f, (0, 1), 1, F(1, 7))
    assert budget > 0 and increment_budget(f, (0, 1), 1, 0) == 0
    assert increment_budget(f, (10, 11), 1, F(1, 7)) == budget
    return dict(reflected_mixture_matrix_psd=True,markov_composition_defect=defect,
                equal_first_moment_counterexamples=atoms,local_increment_budget=budget,
                reconstructed_generator_nonnegative=True,time_zero_markov_claimed=False,
                unique_history_limit_claimed=False,sample_path_continuity_claimed=False)


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
        certificate_type='YM60_INTERACTING_HISTORY_LIMIT_AND_RESPONSE_TO_SPECTRAL_UPPER_BOUND',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        brackets=bracket_controls(),symmetry=symmetry_controls(),cancellation=cancellation_controls(),
        forms=form_controls(),response=response_controls(),histories=history_controls(),
        claim_status='ACTUAL_SCALAR_INTERACTING_SUBSEQUENTIAL_HISTORIES_AND_VARIATIONAL_UPPER_CONSTRAINT__LOWER_GAP_OPEN',
        proof_scope='Written analytic limiting proofs; exact finite controls are not formal verification',
        model=dict(original_scalar_interaction_retained=True,finite_volume_true_vacuum=True,
                   interacting_history_limit='Subsequential after spatial averaging',
                   reflection_positive=True,continuous_heat_time=True,
                   history_limit_unique=False,time_zero_markov=False,
                   uniform_interacting_lower_gap=False,physical_4d_yang_mills=False),
        open='Uniform weighted Poincare inequality for the actual vacuum; uniqueness/Markov limit; physical gauge and continuum'))


def check(cert, result=RESULT, pin=PIN):
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM60 fresh certificate/pin mismatch')
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
    print('YM60 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],trial_rate_series=cert['response']['trial_rate_series'],
                          uniform_interacting_lower_gap=False)))


if __name__ == '__main__':
    main()
