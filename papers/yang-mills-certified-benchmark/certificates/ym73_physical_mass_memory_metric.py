"""YM73: physical gap brackets, memory normalization and slow-mode retention.

Python 3.12 only. Rational matrix examples verify the written spectral bridge;
they do not evaluate a Yang--Mills vacuum or discharge continuum assumptions.
Default/--check are read-only; --write updates only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from math import factorial, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym72_growing_block_matrix_memory as previous

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative = previous.positive, previous.nonnegative
adapted = previous.previous.adapted
matrix, sub = adapted.matrix, adapted.sub
RESULT = HERE/'YM73_RESULT.json'
PIN = HERE/'EXPECTED_YM73.sha256'
SOURCES = HERE/'YM73_SOURCE_PINS.json'
GRID = previous.previous.GRID


def symmetric(rows):
    out = matrix(rows)
    if len(out) != len(out[0]) or out != a.transpose(out):
        raise ValueError('An exact symmetric square matrix is required')
    return out


def positive_definite(rows):
    out = symmetric(rows)
    if not a.psd(out):
        raise ValueError('Positive definite matrix required')
    a.inverse(out)  # Reject a singular positive semidefinite matrix.
    return out


def block_data(C, B, D):
    C, B, D = symmetric(C), matrix(B), positive_definite(D)
    n, m = len(C), len(D)
    if len(B) != n or len(B[0]) != m:
        raise ValueError('Incompatible retained and hidden dimensions')
    Bt, invD = a.transpose(B), a.inverse(D)
    X = a.multiply(invD, Bt)
    M0, M1 = a.multiply(B, X), a.multiply(a.transpose(X), X)
    S0, Z = a.add(C, M0, -1), a.add(a.identity(n), M1)
    if not a.psd(S0):
        raise ValueError('This physical gap criterion requires a nonnegative full operator')
    H = [C[i]+B[i] for i in range(n)]+[Bt[i]+D[i] for i in range(m)]
    W = a.identity(n)+a.scale(X, -1)
    return dict(C=C, B=B, D=D, H=H, X=X, M0=M0, M1=M1, S0=S0, Z=Z,
                harmonic_extension=W, retained_dimension=n, hidden_dimension=m)


def pencil_interval(S, Z, tolerance=F(1, 10**12)):
    """Exact bisection for min <f,Sf>/<f,Zf>; no eigenvector fitting."""
    S, Z, tolerance = symmetric(S), positive_definite(Z), positive(tolerance)
    if len(S) != len(Z) or not a.psd(S):
        raise ValueError('A positive semidefinite numerator of matching size is required')
    lower, upper = F(0), min(S[i][i]/Z[i][i] for i in range(len(S)))
    while upper-lower > tolerance:
        midpoint = (lower+upper)/2
        if a.psd(a.add(S, Z, -midpoint)):
            lower = midpoint
        else:
            upper = midpoint
    return dict(lower=lower, upper=upper)


def gap_bracket(effective_lower, effective_upper, hidden_physical_floor):
    lo, hi, hidden = nonnegative(effective_lower), nonnegative(effective_upper), positive(hidden_physical_floor)
    if hi < lo:
        raise ValueError('Reversed effective energy interval')
    return dict(lower=lo*hidden/(lo+hidden), upper=hi, hidden_physical_floor=hidden)


def gap_certificate(C, B, D, hidden_floor, clock=F(1), tolerance=F(1, 10**12)):
    cell = block_data(C, B, D)
    hidden_floor, clock = positive(hidden_floor), positive(clock)
    if not a.psd(a.add(cell['D'], a.identity(len(cell['D'])), -hidden_floor)):
        raise ValueError('The asserted hidden floor is false')
    effective = pencil_interval(a.scale(cell['S0'], clock), cell['Z'], tolerance)
    return dict(effective_energy=effective,
                full_gap=gap_bracket(effective['lower'], effective['upper'], clock*hidden_floor))


def interval_pencil(C_lower, C_upper, M0_lower, M0_upper, M1_lower, M1_upper,
                    hidden_floor, clock=F(1), probe=None, tolerance=F(1, 10**12)):
    cells = [symmetric(value) for value in (C_lower, C_upper, M0_lower, M0_upper, M1_lower, M1_upper)]
    n = len(cells[0])
    if any(len(cell) != n for cell in cells):
        raise ValueError('Interval matrices must have the same dimension')
    CL, CU, L0, U0, L1, U1 = cells
    for lo, hi in ((CL, CU), (L0, U0), (L1, U1)):
        if not a.psd(a.add(hi, lo, -1)):
            raise ValueError('Reversed Loewner interval')
    ZL, ZU = positive_definite(a.add(a.identity(n), L1)), positive_definite(a.add(a.identity(n), U1))
    SL, SU = a.add(CL, U0, -1), a.add(CU, L0, -1)
    if not a.psd(SL):
        raise ValueError('These intervals do not certify a nonnegative static Schur form')
    clock, hidden_floor = positive(clock), positive(hidden_floor)
    if probe is None:
        probe = [F(i == 0) for i in range(n)]
    probe = [exact(value) for value in probe]
    if len(probe) != n or not any(probe):
        raise ValueError('A nonzero retained probe of the right size is required')
    lo = pencil_interval(a.scale(SL, clock), ZU, tolerance)['lower']
    hi = clock*a.dot(probe, a.apply(SU, probe))/a.dot(probe, a.apply(ZL, probe))
    return dict(effective_lower=lo, effective_upper=hi,
                full_gap=gap_bracket(lo, hi, clock*hidden_floor))


def resolvent_control(C, B, D, hidden_floor, clock, spectral):
    clock, spectral, hidden_floor = positive(clock), positive(spectral), positive(hidden_floor)
    cell = block_data(C, B, D)
    n, m = cell['retained_dimension'], cell['hidden_dimension']
    if not a.psd(a.add(cell['D'], a.identity(m), -hidden_floor)):
        raise ValueError('The asserted hidden floor is false')
    Hs = a.add(a.scale(cell['H'], clock), a.identity(n+m), spectral)
    actual = sub(a.inverse(Hs), tuple(range(n)))
    static = a.add(a.scale(cell['S0'], clock), cell['Z'], spectral)
    approximation = a.inverse(static)
    hidden = a.add(cell['D'], a.identity(m), spectral/clock)
    sigma = a.multiply(a.multiply(cell['B'], a.inverse(hidden)), a.transpose(cell['B']))
    schur = a.add(a.add(a.scale(cell['C'], clock), a.identity(n), spectral), sigma, -clock)
    remainder = a.add(static, schur, -1)
    invD = a.inverse(cell['D'])
    exact_remainder = a.scale(a.multiply(a.multiply(a.multiply(a.multiply(
        cell['B'], invD), invD), a.inverse(hidden)), a.transpose(cell['B'])), spectral*spectral/clock)
    difference = a.add(actual, approximation, -1)
    ceiling = a.scale(a.inverse(cell['Z']), 1/(clock*hidden_floor))
    assert actual == a.inverse(schur)
    assert remainder == exact_remainder
    assert a.psd(remainder)
    assert a.psd(a.add(a.scale(cell['Z'], spectral*spectral/(clock*hidden_floor+spectral)), remainder, -1))
    assert a.psd(difference) and a.psd(a.add(ceiling, difference, -1))
    return dict(spectral=spectral, clock=clock, full_compressed_resolvent=actual,
                two_moment_resolvent=approximation, positive_difference=difference,
                unweighted_difference_upper=ceiling, normalized_norm_error_upper=1/(clock*hidden_floor))


def slow_probe_budget(theta, factor, readouts=1):
    cell = previous.parameters(theta, factor, readouts)
    dynamics = previous.dynamics_budget(theta, factor, readouts)
    residual = dynamics['normalized_generator_residual']
    lost_micro = F(200, 57*factor)
    overlap = previous.up(F(2, 3)*residual+previous.sqrt_upper(lost_micro))
    return dict(microscopic_low_weight_lower=max(F(0), 1-lost_micro),
                large_loop_low_weight_upper=previous.up(F(4, 9)*residual*residual),
                micro_large_overlap_upper=min(F(1), overlap),
                enlarged_Gram_lower=max(F(0), 1-overlap), enlarged_Gram_upper=1+overlap,
                independent_slow_direction_certified=overlap < 1, rho=cell['rho'])


def sqrt_lower(value):
    value = nonnegative(value)
    return F(isqrt(value.numerator*GRID*GRID//value.denominator), GRID)


def physical_escape_budget(spacing, factor, readouts=1, electric=F(1), magnetic=F(1)):
    spacing, electric, magnetic = positive(spacing), positive(electric), positive(magnetic)
    native.width(factor, 3)
    native.width(readouts)
    squared = 6400*electric*magnetic*factor*readouts/(spacing*spacing)
    return dict(physical_gap_lower_squared=squared, physical_gap_lower=sqrt_lower(squared),
                assumption='0 < theta*r*N <= 1/100000000 with the declared fixed positive calibration',
                calibrated_clock_lower_squared=10**8*electric*magnetic*factor*readouts/(spacing*spacing))


def escape_heat_upper(spacing, factor, readouts, time, norm_squared=F(1), order=8):
    time, norm_squared = positive(time), nonnegative(norm_squared)
    native.width(order)
    energy = physical_escape_budget(spacing, factor, readouts)['physical_gap_lower']
    if not energy:
        return norm_squared
    return norm_squared*min(F(1), F(factorial(order))/(time*energy)**order)


def fixture(clock, target=F(3, 2)):
    clock, target = positive(clock), positive(target)
    D = a.diag((2, 5))
    X = matrix(((F(1, 2), F(1, 3)), (F(1, 4), -F(1, 5))))
    B = a.multiply(a.transpose(X), D)
    Z = a.add(a.identity(2), a.multiply(a.transpose(X), X))
    M0 = a.multiply(B, X)
    # This algebraic family has exactly K=target*I for every clock.
    C = a.add(M0, Z, target/clock)
    return block_data(C, B, D)


@lru_cache(maxsize=None)
def native_controls():
    model = adapted.model
    plaquette = model.plaquette_controls()
    assert plaquette['free_electric_rate'] == 3
    geometry = []
    for side in (4, 5):
        cell = model.lattice(side)
        neighbors = {v: set() for v in cell['vertices']}
        undirected = set()
        for left, right in cell['ends']:
            assert left != right
            edge = frozenset((left, right))
            assert edge not in undirected
            undirected.add(edge)
            neighbors[left].add(right)
            neighbors[right].add(left)
        assert all(not neighbors[left].intersection(neighbors[right]) for left, right in cell['ends'])
        assert all(len({edge for edge, sign in face}) == 4 for face in cell['faces'])
        geometry.append(dict(periodic_side=side, no_self_loops_parallel_edges_or_triangles=True,
                             shortest_cycle=4, physical_fundamental_cycle_energy=F(3)))
    return dict(native_plaquette=plaquette, geometry=geometry,
                arbitrary_volume_free_gap_uses_written_gauge_representation_argument=True)


@lru_cache(maxsize=None)
def matrix_controls():
    rows, commutators, resolvents = [], [], []
    for clock in (F(4), F(16), F(64), F(256)):
        cell = fixture(clock)
        W, H, Z, S0 = (cell[key] for key in ('harmonic_extension', 'H', 'Z', 'S0'))
        assert a.multiply(a.transpose(W), W) == Z
        assert a.multiply(a.multiply(a.transpose(W), H), W) == S0
        assert a.scale(S0, clock) == a.scale(Z, F(3, 2))
        commutator = a.add(a.multiply(cell['C'], cell['M1']), a.multiply(cell['M1'], cell['C']), -1)
        assert any(value for row in commutator for value in row)
        commutators.append(commutator)
        bounds = gap_bracket(F(3, 2), F(3, 2), 2*clock)
        physical = a.scale(H, clock)
        assert a.psd(a.add(physical, a.identity(4), -bounds['lower']))
        actual = pencil_interval(physical, a.identity(4))
        assert actual['upper'] >= bounds['lower'] and actual['lower'] <= bounds['upper']
        rows.append(dict(clock=clock, exact_effective_energy=F(3, 2), physical_hidden_floor=2*clock,
                         full_gap_bracket=bounds, independently_bisected_full_gap=actual))
        for spectral in (F(1, 100), F(1), F(100), F(10**6)):
            resolvents.append(resolvent_control(cell['C'], cell['B'], cell['D'], 2, clock, spectral))
        scaled = block_data(a.scale(cell['C'], clock), a.scale(cell['B'], clock), a.scale(cell['D'], clock))
        assert scaled['M0'] == a.scale(cell['M0'], clock)
        assert scaled['M1'] == cell['M1'] and scaled['Z'] == Z
    return dict(fixture_scope='Positive rational operator families, not Yang--Mills vacua',
                gap_convergence=rows, nonzero_static_metric_commutators=commutators,
                exact_resolvent_checks=resolvents)


def interval_controls():
    rows = []
    for clock in (F(4), F(16), F(64)):
        cell = fixture(clock)
        radius = F(1, 100000)/clock
        I = a.identity(2)
        CL, CU = a.add(cell['C'], I, -radius), a.add(cell['C'], I, radius)
        L0, U0 = a.add(cell['M0'], I, -radius), a.add(cell['M0'], I, radius)
        L1, U1 = a.add(cell['M1'], I, -radius), a.add(cell['M1'], I, radius)
        cert = interval_pencil(CL, CU, L0, U0, L1, U1, 2, clock)
        assert cert['effective_lower'] <= F(3, 2) <= cert['effective_upper']
        assert cert['effective_upper']-cert['effective_lower'] < F(1, 1000)
        assert a.psd(a.add(a.scale(cell['H'], clock), a.identity(4), -cert['full_gap']['lower']))
        rows.append(dict(clock=clock, per_matrix_radius=radius, **cert))
    return rows


def sharpness_controls():
    h, target = F(2), F(3, 2)
    boundary = target*h/(target+h)
    rows = []
    for x in (F(1), F(4), F(16), F(64)):
        H = matrix(((h*x*x+target*(1+x*x), h*x), (h*x, h)))
        exact_trace = (h+target)*(1+x*x)
        exact_det = h*target*(1+x*x)
        assert H[0][0]+H[1][1] == exact_trace
        assert H[0][0]*H[1][1]-H[0][1]**2 == exact_det
        actual = pencil_interval(H, a.identity(2))
        assert actual['lower'] >= boundary
        rows.append(dict(mixing=x, lower_bound=boundary, actual_gap_interval=actual,
                         distance_upper=actual['upper']-boundary))
    assert all(rows[i]['distance_upper'] > rows[i+1]['distance_upper'] for i in range(len(rows)-1))
    return rows


def scope_counterexamples():
    # Missing first moment changes a limiting mass by a factor of five.
    wrong_normalization = []
    for clock in (F(4), F(16), F(64), F(256)):
        D, X = [[F(1)]], [[F(2)]]
        B, Z, M0 = [[F(2)]], [[F(5)]], [[F(4)]]
        C = a.add(M0, Z, F(3, 2)/clock)
        cell = block_data(C, B, D)
        actual = pencil_interval(a.scale(cell['H'], clock), a.identity(2))
        wrong_normalization.append(dict(clock=clock, correct_effective_energy=F(3, 2),
            unnormalized_static_energy=F(15, 2), actual_gap_interval=actual))
    # The retained scalar response and every memory moment are identical,
    # while an uncoupled hidden state can have any positive small energy.
    invisible = []
    for delta in (F(1, 2), F(1, 16), F(1, 256)):
        cell = block_data([[F(3)]], [[F(0)]], [[delta]])
        cert = gap_certificate(cell['C'], cell['B'], cell['D'], delta)
        assert cell['M0'] == [[0]] and cell['M1'] == [[0]]
        assert a.inverse(a.add(cell['H'], a.identity(2)))[0][0] == F(1, 4)
        assert cert['full_gap']['lower'] <= delta
        invisible.append(dict(hidden_gap=delta, retained_rate=F(3), all_memory_moments_zero=True,
                              retained_response_at_one=F(1, 4), honest_gap_bracket=cert['full_gap']))
    # A finite lattice gap does not ensure the observed finite-energy state survives.
    visibility = []
    for fast in (4, 16, 64):
        # At t=log(2), observed heat is exactly 2^(-fast), while full gap is one.
        visibility.append(dict(full_gap=F(1), selected_state_norm_squared=F(1),
                               selected_state_energy=F(fast), correlation_at_log_two=F(1, 2**fast)))
    return dict(fixtures_are_not_Yang_Mills_solutions=True, missing_first_moment=wrong_normalization,
                invisible_sector=invisible, finite_gap_without_local_visibility=visibility)


def slow_and_scale_controls():
    slow = []
    for factor in (16, 64, 256, 1024):
        count = factor
        theta = F(1, 10**8*factor*count)
        row = slow_probe_budget(theta, factor, count)
        assert row['enlarged_Gram_lower'] > F(53, 100)
        assert row['enlarged_Gram_upper'] < F(147, 100)
        assert row['independent_slow_direction_certified']
        slow.append(dict(block_factor=factor, readouts=count, theta=theta, **row))
    assert slow[0]['microscopic_low_weight_lower'] == F(89, 114)
    scales = []
    for inverse_spacing in (4, 16, 64, 256):
        spacing = F(1, inverse_spacing)
        row = physical_escape_budget(spacing, 4, 1)
        assert row['physical_gap_lower'] == 160*inverse_spacing
        upper = escape_heat_upper(spacing, 4, 1, 1)
        scales.append(dict(spacing=spacing, **row, time_one_correlation_upper=upper))
    assert all(scales[i]['time_one_correlation_upper'] > scales[i+1]['time_one_correlation_upper']
               for i in range(len(scales)-1))
    # General positive convention constants: kappa^2*theta*a^2=c_E*c_B.
    calibrations = []
    for theta, spacing, electric, magnetic in ((F(1, 10**8), F(1, 10), F(1), F(1)),
            (F(1, 10**12), F(1, 100), F(2, 3), F(7, 5))):
        squared_clock = electric*magnetic/(theta*spacing*spacing)
        assert squared_clock*theta*spacing*spacing == electric*magnetic
        calibrations.append(dict(theta=theta, spacing=spacing, electric=electric, magnetic=magnetic,
                                 physical_clock_squared=squared_clock))
    return dict(actual_model_slow_probe_enclosures=slow, current_trajectory_spectrum_escape=scales,
                calibration_identities=calibrations,
                all_slow_modes_retained_by_one_added_probe=False)


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
    return previous.spatial.encode(dict(
        certificate_type='YM73_PHYSICAL_MASS_MEMORY_METRIC_AND_CONTINUUM_CRITERION',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        native=native_controls(), matrices=matrix_controls(), interval_pencils=interval_controls(),
        sharp_lower_boundary=sharpness_controls(), counterexamples=scope_counterexamples(),
        actual_model_budgets=slow_and_scale_controls(),
        model=dict(actual_microscopic_slow_probe_enrichment=True,
                   full_physical_observation_cut_declared=True,
                   untruncated_hidden_space_in_written_theorem=True,
                   memory_normalization_Z_equals_I_plus_M1=True,
                   noncircular_conditional_physical_gap_bracket=True,
                   uniform_normalized_resolvent_error_one_over_hidden_physical_floor=True,
                   actual_YM72_calibrated_trajectory_has_escaping_spectrum=True,
                   continuum_gap_passage_on_a_dense_core_is_conditional=True,
                   independently_certified_hidden_floor_on_weak_coupling_trajectory=False,
                   actual_weak_coupling_static_pencil_computed=False,
                   full_low_energy_space_retained_by_one_extra_probe=False,
                   finite_fixture_is_Yang_Mills_vacuum=False,
                   lattice_gap_convergence_alone_constructs_continuum_fields=False,
                   continuum_R4_QFT=False, physical_mass_prediction=False),
        claim_status='PHYSICAL_SPECTRAL_BRIDGE_PROVED__YANG_MILLS_CONTINUUM_INPUTS_REMAIN_OPEN',
        proof_scope='Written operator and scaling proofs with exact rational matrix controls and directed native budgets',
        open='Weak-coupling whole-pencil bounds, independent hidden separation, surviving local states and 4D field construction'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM73 fresh certificate/pin mismatch')
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
    print('YM73 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], metric='I+M1',
                          physical_gap_bracket='m_eff*h/(m_eff+h) <= gap <= m_eff',
                          continuum_inputs_discharged=False, continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
