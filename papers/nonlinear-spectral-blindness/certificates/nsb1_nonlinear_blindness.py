"""NSB1 nonlinear spectral-blindness evidence, Python 3.12 only.

Exact finite jet, response-geometry and compact-gauge controls accompany
written proofs. --write changes only NSB1_RESULT.json and EXPECTED_NSB1.sha256.
Default and --check reproduce evidence without rewriting it.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import nsb1_jets as jets
import nsb1_geometry as geometry

RESULT = HERE/'NSB1_RESULT.json'
PIN = HERE/'EXPECTED_NSB1.sha256'
SOURCES = HERE/'NSB1_SOURCE_MANIFEST.json'
BASE_TREE = 'cf7c60e802fed9f87823d473ea303007783167f5'
MODEL = dict(
    nonlinear_common_fixed_point_jet_characteristic_theorem=True,
    all_unmarked_jet_word_spectra_recover_nonlinear_coefficients=False,
    Jordan_structure_claimed_blind=False,
    coefficient_markers_recover_declared_finite_jets=True,
    finite_probe_minimum_requires_the_declared_chart_and_C1_decoder=True,
    arbitrary_smooth_germs_recovered_from_all_jets=False,
    entire_Hessian_eigenvalue_field_determines_flatness=False,
    entire_Maxwell_speed_field_determines_calibrated_curvature=False,
    flat_curved_pair_has_gauge_inequivalent_selected_connections=True,
    exponential_pair_has_gauge_inequivalent_selected_connections=False,
    all_physical_symmetries_imply_all_curvatures_zero=False,
    dual_flat_affine_connections_force_their_mean_flat=False,
    general_nonlinear_Onsager_reciprocity_derived=False,
    native_gauge_parent_matches_free_potential_jets_to_arbitrarily_high_finite_order=True,
    original_homogeneous_Yang_Mills_Hamiltonian_modified=False,
    native_counterfamily_shares_the_full_quantum_spectrum=False,
    negative_compression_bound_from_SPECTRAL2_used_uncorrected=False,
    infinite_Koopman_spectrum_identified=False,
    experimental_material_effect_validated=False,
    positive_four_dimensional_continuum_Yang_Mills_gap_proved=False,
    finite_controls_formally_prove_universal_analytic_theorems=False,
)


def runtime_check():
    jets.runtime_check()
    if sys.flags.optimize:
        raise RuntimeError('Run without optimization: evidence assertions must execute')


def nonnegative(value):
    value = jets.exact(value)
    if value < 0:
        raise ValueError('A nonnegative exact parameter is required')
    return value


def exp_negative_interval(value, terms=None):
    """Reciprocal positive Taylor series with a geometric remainder."""
    value = nonnegative(value)
    if terms is None:
        terms = max(24, 2*value.__ceil__())
    jets.width(terms, 1)
    if value >= terms+2:
        raise ValueError('The geometric exponential tail ratio must be below one')
    total, term = F(1), F(1)
    for k in range(1, terms+1):
        term *= value/k
        total += term
    first_omitted = term*value/(terms+1)
    tail_upper = first_omitted/(1-value/(terms+2))
    return dict(lower=1/(total+tail_upper), upper=1/total, terms=terms)


def gauge_parent(N, t):
    """One physical four-link plaquette; no overlapping-vacuum ansatz."""
    jets.width(N)
    t = nonnegative(t)
    w = {(1,): F(1)}
    one = {(0,): F(1)}
    w2 = jets.poly_mul(w, w)
    s = jets.poly_add(one, jets.poly_scale(w2, -1))
    f = jets.poly_add(one, jets.poly_scale(jets.poly_pow(s, N+1), -1))
    df = jets.poly_derivative(f, 0)
    ddf = jets.poly_derivative(df, 0)
    electric = jets.poly_add(jets.poly_scale(jets.poly_mul(w, df), 3),
                             jets.poly_scale(jets.poly_mul(s, ddf), -1))
    gamma = jets.poly_mul(s, jets.poly_mul(df, df))
    expected_electric = jets.poly_scale(
        jets.poly_mul(jets.poly_pow(s, N),
                      jets.poly_add(jets.poly_scale(w2, 2*N+4),
                                    jets.poly_scale(one, -1))), 2*(N+1))
    expected_gamma = jets.poly_scale(
        jets.poly_mul(w2, jets.poly_pow(s, 2*N+1)), 4*(N+1)**2)
    assert electric == expected_electric and gamma == expected_gamma
    potential = jets.poly_add(jets.poly_scale(electric, -t/2),
                              jets.poly_scale(gamma, t*t/4))

    # Independently express this even polynomial in the local radial variable
    # s=1-W^2. Since s=O(|A|^2), its first power gives the local jet order.
    radial = {}
    S = {(1,): F(1)}
    for (power,), coefficient in potential.items():
        assert power % 2 == 0
        radial = jets.poly_add(radial, jets.poly_scale(
            jets.poly_pow(jets.poly_add(one, jets.poly_scale(S, -1)), power//2),
            coefficient))
    first_radial_order = min((power[0] for power in radial), default=None)
    if t:
        assert first_radial_order == N
        assert radial[(N,)] == -t*(N+1)*(2*N+3)
    else:
        assert radial == {}
    cN = F(3, 4)**(N+1)-F(1, 4)**(N+1)
    assert cN > 0
    exponential = exp_negative_interval(t*cN)
    cap = F(2, 81)
    assert 4/cap == 162
    return dict(N=N, t=t, f=jets.polynomial_records(f),
                electric_f=jets.polynomial_records(electric),
                Gamma_f=jets.polynomial_records(gamma),
                multiplication_potential=jets.polynomial_records(potential),
                potential_in_s=jets.polynomial_records(radial),
                first_radial_order=first_radial_order,
                first_possible_local_coordinate_order=2*N,
                matched_free_potential_jet_order=2*N-1 if N else None,
                potential_matches_all_orders_when_t_zero=(t == 0),
                exact_ground_state='h proportional to exp(t f_N/2)',
                exact_ground_energy=F(0),
                numerator_density_exponent=1-F(3, 4)**(N+1),
                denominator_density_exponent=1-F(1, 4)**(N+1),
                Haar_cap_lower=cap, decay_exponent=cN,
                exponential_interval=exponential,
                physical_gauge_invariant_gap_upper=162*exponential['upper'],
                extra_W_even_restriction_would_exclude_the_trial=True,
                exact_vacuum_requires_one_or_edge_disjoint_active_plaquettes=True,
                original_homogeneous_Yang_Mills_model=False)


def corrected_compression_bound(restricted_supremum, matched_tail):
    return max(F(0), jets.exact(restricted_supremum))+nonnegative(matched_tail)


def source_erratum_controls():
    mu, tail, actual_top = F(-100), F(2), F(2)
    incorrect = mu+tail
    corrected = corrected_compression_bound(mu, tail)
    assert incorrect <= 1 < actual_top and corrected == actual_top
    return dict(source='SPECTRAL_2.pdf, Section 16, printed pages 28-29',
                compact_diagonal_counterexample=[-100, 2, 0],
                projection_indices=[0], restricted_supremum=mu,
                matched_tail_norm=tail, actual_spectral_supremum=actual_top,
                printed_upper=-98, corrected_upper=corrected,
                printed_condition_would_falsely_accept=True,
                corrected_condition_accepts=False,
                I_minus_K_negative_eigenvalue=F(-1),
                correction='max(0,mu_n)+norm(K-P_n K P_n)',
                correct_Schur_complement_identity_rejected=False)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['base_tree'] != BASE_TREE or manifest['stage'] != 'NSB1':
        raise ValueError('NSB1 source base or stage mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed upstream source: '+path)
        checked[path] = actual
    if not manifest['upstream_sha256'] or len(set(manifest['local_inputs'])) != len(manifest['local_inputs']):
        raise ValueError('Invalid NSB1 source manifest')
    for path in manifest['local_inputs']:
        checked[path] = hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = hashlib.sha256(SOURCES.read_bytes()).hexdigest()
    return checked


@lru_cache(maxsize=1)
def all_controls():
    gauge = [gauge_parent(N, t) for N in (0, 1, 2, 4) for t in (F(0), F(16), F(64))]
    assert gauge_parent(2, 16)['decay_exponent'] == F(13, 32)
    assert gauge_parent(2, 16)['physical_gauge_invariant_gap_upper'] < 1
    assert gauge_parent(0, 0)['decay_exponent'] == F(1, 2)
    return dict(jet_and_mechanical_controls=jets.all_controls(),
                nonlinear_response_geometry=geometry.all_controls(),
                compact_gauge_finite_jet_counterfamily=gauge,
                source_finite_section_correction=source_erratum_controls(),
                jet_scope=jets.MODEL, geometry_scope=geometry.MODEL)


def canonical_sha(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def run():
    runtime_check()
    return jets.encode(dict(certificate_type='NSB1_NONLINEAR_SPECTRAL_BLINDNESS',
        verdict='PASS', runtime_policy='Python 3.12 only', model=MODEL,
        input_sha256=source_checks(), controls=all_controls(),
        claim_status='NONLINEAR_JET_BLINDNESS_AND_CALIBRATED_PHYSICS_RECOVERY__CONTINUUM_YM_OPEN',
        proof_scope='Written universal theorems with exact finite arithmetic and independent controls',
        open='Actual homogeneous weak-coupling vacuum estimates and four-dimensional continuum Yang-Mills gap'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'NSB1_NONLINEAR_SPECTRAL_BLINDNESS':
        raise ValueError('NSB1 model or scope overclaim')
    if cert.get('controls', {}).get('jet_scope') != jets.MODEL or cert.get('controls', {}).get('geometry_scope') != geometry.MODEL:
        raise ValueError('NSB1 component scope mismatch')
    sha = canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('NSB1 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--write', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('NSB1 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],
        same_spectral_field_flat_curvature='0', same_spectral_field_curved_anchor='1/16',
        gauge_N2_upper='162 exp(-13t/32)', continuum_Yang_Mills_mass_gap=False)))


if __name__ == '__main__':
    main()
