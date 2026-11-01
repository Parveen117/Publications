"""YM75: source-resolved Schur counts, quadratic solve errors and band decoding.

Python 3.12 only. Actual four-cube enclosures use the complete YM74 physical
cutoff; other finite matrices are labelled algebra controls, not gauge vacua.
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
import ym74_all_coupling_electric_enclosures as previous

native, p, a, exact = previous.native, previous.p, previous.a, previous.exact
positive, nonnegative = previous.positive, previous.nonnegative
matrix, symmetric, sub, encode = previous.matrix, previous.symmetric, previous.sub, previous.encode
RESULT = HERE/'YM75_RESULT.json'
PIN = HERE/'EXPECTED_YM75.sha256'
SOURCES = HERE/'YM75_SOURCE_PINS.json'
BASE_TREE = '24bbd9072bd1c7070d95a929f7ea6d06d3e088e3'
MODEL = dict(source_resolved_matrix_Schur_lower_bound=True,
             exact_full_retained_inertia_not_scalarized=True,
             quadratic_residual_Loewner_memory_enclosures=True,
             shifted_hidden_resolvents_supported=True,
             stable_band_decoder_under_verified_leakage_bound=True,
             actual_four_cube_theta_one_over_48_improvement=True,
             all_241_physical_retained_directions_counted=True,
             hidden_auxiliary_modes_may_be_discarded=False,
             source_injectivity_alone_implies_uniform_gap=False,
             finite_source_norm_determines_source_orientation=False,
             low_energy_coverage_proves_nonvacuum_state_exists=False,
             uniform_weak_coupling_positive_gap=False,
             continuum_R4_QFT=False, physical_mass_prediction=False,
             matrix_fixtures_are_Yang_Mills_vacua=False)


def inertia(rows):
    """Exact rational symmetric congruence, including zero-diagonal 2x2 pivots."""
    work = symmetric(rows)
    negative = zero = positive_count = 0
    while work:
        n = len(work)
        pivot = next((i for i in range(n) if work[i][i]), None)
        if pivot is not None:
            value = work[pivot][pivot]
            negative += int(value < 0)
            positive_count += int(value > 0)
            keep = [i for i in range(n) if i != pivot]
            work = [[work[i][j]-work[i][pivot]*work[pivot][j]/value for j in keep] for i in keep]
            continue
        pair = next(((i, j) for i in range(n) for j in range(i+1, n) if work[i][j]), None)
        if pair is None:
            zero += n
            break
        left, right = pair
        value = work[left][right]
        keep = [i for i in range(n) if i not in pair]
        work = [[work[i][j]-(work[i][left]*work[right][j]+work[i][right]*work[left][j])/value
                 for j in keep] for i in keep]
        negative += 1
        positive_count += 1
    return dict(negative=negative, zero=zero, positive=positive_count)


def scalar_inertia(values):
    values = tuple(map(exact, values))
    return dict(negative=sum(value < 0 for value in values),
                zero=sum(value == 0 for value in values),
                positive=sum(value > 0 for value in values))


def source_resolved_schur(C, G, hidden_floor, energy, index=0):
    """C-zI-G/(d-z), under the independently proved premises D>=d and BB*<=G."""
    C, G = symmetric(C), symmetric(G)
    d, z = exact(hidden_floor), exact(energy)
    native.width(index, 0)
    if len(C) != len(G) or index >= len(C):
        raise ValueError('Matching matrix dimensions and valid retained index required')
    if not a.psd(G):
        raise ValueError('A positive semidefinite source envelope is required')
    if z >= d:
        raise ValueError('The tested energy must lie strictly below the hidden floor')
    lower = a.add(a.add(C, a.identity(len(C)), -z), G, -1/(d-z))
    counts = inertia(lower)
    return dict(energy=z, hidden_floor=d, lower_Schur_matrix=lower, inertia=counts,
                eigenvalue_index=index, lower_eigenvalue_certified=counts['negative'] <= index,
                requires_independent_hidden_and_source_envelopes=True)


def finite_block_control(C, B, D, hidden_floor, energy, index=0, envelope=None):
    """Finite-fixture check that the claimed premises are true before counting."""
    C, B, D = symmetric(C), matrix(B), symmetric(D)
    n, m = len(C), len(D)
    if len(B) != n or len(B[0]) != m:
        raise ValueError('Incompatible block dimensions')
    d, z = exact(hidden_floor), exact(energy)
    if not a.psd(a.add(D, a.identity(m), -d)):
        raise ValueError('False hidden lower bound')
    BB = a.multiply(B, a.transpose(B))
    G = BB if envelope is None else symmetric(envelope)
    if len(G) != n or not a.psd(a.add(G, BB, -1)):
        raise ValueError('False source envelope')
    bound = source_resolved_schur(C, G, d, z, index)
    inverse = a.inverse(a.add(D, a.identity(m), -z))
    exact_Schur = a.add(a.add(C, a.identity(n), -z),
                        a.multiply(a.multiply(B, inverse), a.transpose(B)), -1)
    assert a.psd(a.add(exact_Schur, bound['lower_Schur_matrix'], -1))
    H = [C[i]+B[i] for i in range(n)]+[a.transpose(B)[i]+D[i] for i in range(m)]
    full_counts = inertia(a.add(H, a.identity(n+m), -z))
    exact_counts = inertia(exact_Schur)
    assert full_counts['negative'] == exact_counts['negative'] <= bound['inertia']['negative']
    assert full_counts['zero'] == exact_counts['zero']
    return dict(**bound, exact_full_inertia=full_counts, exact_Schur=exact_Schur,
                source_envelope=G, actual_source_Gram=BB,
                algebraic_fixture_not_gauge_vacuum=True)


def residual_enclosure(lower_form, residual_Gram, shifted_hidden_floor, clock=1):
    """Memory T<=M<=T+R*R/d, with T and R*R independently supplied exactly."""
    T, gram = symmetric(lower_form), symmetric(residual_Gram)
    d, clock = positive(shifted_hidden_floor), positive(clock)
    if len(T) != len(gram) or not a.psd(gram):
        raise ValueError('Matching dimensions and a positive residual Gram are required')
    error = a.scale(gram, 1/d)
    return dict(lower=T, upper=a.add(T, error), error_Gram_upper=error,
                physical_norm_error_upper=clock*sum(error[i][i] for i in range(len(error))),
                shifted_hidden_floor=d, clock=clock)


def finite_residual_memory(D, sources, approximate_solve, hidden_floor, shift=0, clock=1):
    """Exact finite controls for F* (D-z)^-1 F and a finite solve X."""
    D, sources, X = symmetric(D), matrix(sources), matrix(approximate_solve)
    d, z = positive(hidden_floor), exact(shift)
    if len(sources) != len(D) or len(X) != len(D) or len(X[0]) != len(sources[0]):
        raise ValueError('Source and solve matrices must have the hidden dimension')
    if not a.psd(a.add(D, a.identity(len(D)), -d)) or z >= d:
        raise ValueError('A true positive shifted hidden floor is required')
    A = a.add(D, a.identity(len(D)), -z)
    R = a.add(sources, a.multiply(A, X), -1)
    Ft, Xt = a.transpose(sources), a.transpose(X)
    T = a.add(a.add(a.multiply(Ft, X), a.multiply(Xt, sources)),
              a.multiply(a.multiply(Xt, A), X), -1)
    inverse = a.inverse(A)
    M = a.multiply(a.multiply(Ft, inverse), sources)
    residual_Gram = a.multiply(a.transpose(R), R)
    exact_error = a.multiply(a.multiply(a.transpose(R), inverse), R)
    assert M == a.add(T, exact_error)
    result = residual_enclosure(T, residual_Gram, d-z, clock)
    assert a.psd(a.add(M, result['lower'], -1))
    assert a.psd(a.add(result['upper'], M, -1))
    return dict(**result, exact_memory=M, residual=R, exact_error=exact_error,
                source=sources, approximate_solve=X, shift=z,
                algebraic_fixture_not_gauge_vacuum=True)


def decoder_budget(leakage_squared):
    zeta = nonnegative(leakage_squared)
    if zeta >= 1:
        raise ValueError('Stable decoding requires a strictly subunit leakage bound')
    return dict(leakage_squared=zeta, retained_Gram_floor=1-zeta,
                inverse_decoder_norm_squared_upper=1/(1-zeta),
                theorem_is_on_entire_declared_band=True,
                nonvacuum_state_existence_asserted=False)


def cube_parameters(theta):
    theta = exact(theta)
    return dict(theta=theta, d=F(9, 2)-192*abs(theta),
                vacuum_uniform_coupling_squared=48*theta*theta,
                q_uniform=F(383, 4)*theta*theta,
                q_transverse=F(191, 4)*theta*theta,
                q_winding=48*theta*theta)


def cube_schur_inertia(theta, energy):
    """Exact count of all 241 retained directions; the hidden sector is infinite."""
    cell, z = cube_parameters(theta), exact(energy)
    d, k2 = cell['d'], cell['vacuum_uniform_coupling_squared']
    if z >= d:
        raise ValueError('Four-cube threshold must lie strictly below d')
    transverse = 3-z-cell['q_transverse']/(d-z)
    winding = 3-z-cell['q_winding']/(d-z)
    if z:
        pivots = [-z, 3-z-cell['q_uniform']/(d-z)+k2/z]
        uniform = scalar_inertia(pivots)
    elif k2:
        pivots, uniform = None, dict(negative=1, zero=0, positive=1)
    else:
        pivots, uniform = [F(0), F(3)], dict(negative=0, zero=1, positive=1)
    counts = {key: uniform[key]+191*scalar_inertia([transverse])[key]+48*scalar_inertia([winding])[key]
              for key in ('negative', 'zero', 'positive')}
    assert sum(counts.values()) == 241
    return dict(**cell, energy=z, uniform_pivots=pivots, uniform_inertia=uniform,
                transverse_pivot=transverse, transverse_multiplicity=191,
                winding_pivot=winding, winding_multiplicity=48, inertia=counts,
                all_retained_directions=241)


def cube_threshold_bisection(theta, index, lower, upper, tolerance=F(1, 10**12)):
    native.width(index, 0)
    lower, upper, tolerance = exact(lower), exact(upper), positive(tolerance)
    if index >= 241 or lower >= upper:
        raise ValueError('Valid eigenvalue index and ordered thresholds required')
    if cube_schur_inertia(theta, lower)['inertia']['negative'] > index:
        raise ValueError('Initial lower threshold is not certified')
    if cube_schur_inertia(theta, upper)['inertia']['negative'] <= index:
        raise ValueError('Initial upper threshold does not bracket a comparison root')
    while upper-lower > tolerance:
        middle = (lower+upper)/2
        if cube_schur_inertia(theta, middle)['inertia']['negative'] <= index:
            lower = middle
        else:
            upper = middle
    return dict(certified_full_eigenvalue_lower=lower,
                comparison_threshold_upper=upper, width=upper-lower, index=index)


@lru_cache(maxsize=None)
def actual_cube_controls():
    geometry = previous.geometry_controls()
    assert geometry['physical_R3_dimension'] == 241
    theta, z0, z1 = F(1, 48), -F(712209, 10**8), F(48376192, 10**8)
    ground, first = cube_schur_inertia(theta, z0), cube_schur_inertia(theta, z1)
    assert ground['inertia'] == dict(negative=0, zero=0, positive=241)
    assert first['inertia'] == dict(negative=1, zero=0, positive=240)
    old = previous.four_cube_gap(theta)
    c0_upper = -F(692844, 10**8)
    assert c0_upper >= old['retained_ground_interval']['upper']
    lower, upper = z1-c0_upper, 3-z0
    assert lower == F(49069036, 10**8) and upper == F(300712209, 10**8)
    assert lower > old['gap']['lower'] and upper < old['gap']['upper']
    roots = [cube_threshold_bisection(theta, 0, -F(1, 100), -F(1, 1000)),
             cube_threshold_bisection(theta, 1, F(0), F(49, 100))]
    assert roots[0]['certified_full_eigenvalue_lower'] >= z0
    assert roots[1]['certified_full_eigenvalue_lower'] >= z1
    return dict(theta=theta, actual_periodic_cube_side=4,
                actual_complete_physical_retained_dimension=241,
                ground_threshold=ground, first_excited_threshold=first,
                retained_ground_energy_upper=c0_upper,
                certified_bare_gap=dict(lower=lower, upper=upper, width=upper-lower),
                previous_scalarized_gap=old['gap'], comparison_root_brackets=roots,
                full_hidden_space_untruncated=True, actual_interacting_vacuum_computed=False,
                volume_uniform_gap_asserted=False, continuum_mass_asserted=False)


@lru_cache(maxsize=None)
def matrix_controls():
    controls = []
    for row in previous.matrix_controls()['rational_fixtures_not_Yang_Mills_vacua']:
        C, B, D = row['retained_matrix'], row['coupling_matrix'], row['hidden_matrix']
        floor = previous.two_by_two_eigenvalues(D)[0]['lower']
        for z in (min(F(-10), floor-1), min(F(-1), floor-1), floor-F(1, 10)):
            for index in (0, 1):
                cell = finite_block_control(C, B, D, floor, z, index)
                exact_counts = scalar_inertia([F(value)-z for value in row['exact_full_spectrum']])
                assert cell['exact_full_inertia'] == exact_counts
                if cell['lower_eigenvalue_certified']:
                    assert row['exact_full_spectrum'][index] >= z
                controls.append(cell)
    # Zero source leaves genuine low hidden states; they cannot be discarded.
    invisible = finite_block_control([[3, 0], [0, 5]], [[0], [0]], [[1]], 1, F(9, 10), 1)
    assert invisible['lower_eigenvalue_certified']
    # Same ||B||^2, different source orientation; the threshold counts differ.
    C, D, z = [[0, 0], [0, 3]], [[F(1, 2)]], F(49, 100)
    orient = [finite_block_control(C, B, D, F(1, 2), z, 1)
              for B in ([[F(1, 5)], [0]], [[0], [F(1, 5)]])]
    assert orient[0]['inertia']['negative'] == 1 and orient[1]['inertia']['negative'] == 2
    assert orient[0]['lower_eigenvalue_certified'] and not orient[1]['lower_eigenvalue_certified']
    return dict(full_matrix_inertia_controls=controls, invisible_uncoupled_hidden_sector=invisible,
                equal_norm_different_source_orientation=orient,
                finite_fixtures_are_actual_Yang_Mills=False)


@lru_cache(maxsize=None)
def residual_controls():
    D = matrix(((2, 1), (1, 3)))
    sources = matrix(((1, F(1, 3)), (F(-1, 2), 2)))
    rows = []
    for shift in (-F(1), F(0), F(1, 2)):
        A = a.add(D, a.identity(2), -shift)
        exact_solve = a.multiply(a.inverse(A), sources)
        direction = matrix(((F(1, 7), F(-1, 3)), (F(1, 5), F(1, 11))))
        previous_error = None
        for scale in (F(1), F(1, 2), F(1, 4), F(0)):
            X = a.add(exact_solve, direction, scale)
            row = finite_residual_memory(D, sources, X, 1, shift=shift, clock=3)
            if previous_error is not None and scale:
                assert row['exact_error'] == a.scale(previous_error, F(1, 4))
            if not scale:
                assert row['lower'] == row['upper'] == row['exact_memory']
            previous_error = row['exact_error']
            rows.append(dict(**row, solve_error_scale=scale))
    # Uniform closeness of the symmetrized raw Gram is not a PSD lower certificate.
    bad_raw = finite_residual_memory([[2]], [[1]], [[1]], 2)
    assert bad_raw['exact_memory'] == [[F(1, 2)]] and bad_raw['lower'] == [[0]]
    assert [[1]] != bad_raw['lower']
    return dict(quadratic_residual_examples=rows,
                symmetrized_raw_Gram_not_a_lower_bound=bad_raw,
                raw_Gram_value_for_counterexample=[[F(1)]])


@lru_cache(maxsize=None)
def decoder_controls():
    # Two orthonormal band states, retaining only their first two coordinates.
    J = matrix(((F(3, 5), 0), (0, F(4, 5)), (F(4, 5), 0), (0, F(3, 5))))
    assert a.multiply(a.transpose(J), J) == a.identity(2)
    retained = sub(J, (0, 1), (0, 1))
    hidden = sub(J, (2, 3), (0, 1))
    retained_Gram = a.multiply(a.transpose(retained), retained)
    hidden_Gram = a.multiply(a.transpose(hidden), hidden)
    budget = decoder_budget(F(16, 25))
    decoder = a.inverse(retained)
    assert a.multiply(decoder, retained) == a.identity(2)
    decoder_Gram = a.multiply(a.transpose(decoder), decoder)
    assert max(decoder_Gram[i][i] for i in range(2)) == budget['inverse_decoder_norm_squared_upper']
    # Injectivity at every finite stage does not imply a uniform inverse margin.
    collapsing = []
    for n in (2, 4, 8, 16):
        zeta = 1-F(1, n*n)
        row = decoder_budget(zeta)
        assert row['inverse_decoder_norm_squared_upper'] == n*n
        collapsing.append(dict(stage=n, retained_scalar=F(1, n), **row))
    actual_capture = previous.combined_cutoff(256, 192, F(1, 100), 64,
                                               reference_energy=2, physical_band=3)
    actual = decoder_budget(actual_capture['actual_budget']['entire_low_energy_subspace_leakage_squared_upper'])
    return dict(band_isometry=J, retained_Gram=retained_Gram, hidden_Gram=hidden_Gram,
                exact_band_decoder=decoder, sharp_budget=budget,
                injective_without_uniform_stability=collapsing,
                actual_YM74_cutoff_decoder_bound=actual,
                state_existence_or_continuum_limit_not_implied=True)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if len(manifest['upstream_sha256']) != 260 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM75 source manifest shape mismatch')
    if manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM75 base tree mismatch')
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
    return encode(dict(certificate_type='YM75_RH_SOURCE_RESOLVED_TRANSFER',
        verdict='PASS', runtime_policy='Python 3.12 only', input_sha256=source_checks(),
        actual_cube=actual_cube_controls(), spectral_controls=matrix_controls(),
        residual_memory=residual_controls(), band_decoding=decoder_controls(), model=MODEL,
        claim_status='ACTUAL_FINITE_GAP_IMPROVED_AND_ERROR_TRANSFER_PROVED__CONTINUUM_OPEN',
        proof_scope='Written Schur and decoder proofs with exact native geometry and rational matrix controls',
        open='Uniform positive physical retained edge and nontrivial local four-dimensional field construction'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM75_RH_SOURCE_RESOLVED_TRANSFER':
        raise ValueError('YM75 scope or model claim mismatch')
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM75 fresh certificate/pin mismatch')
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
    print('YM75 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0], source_resolved_inertia=True,
                          actual_four_cube_gap_interval=['0.49069036', '3.00712209'],
                          continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
