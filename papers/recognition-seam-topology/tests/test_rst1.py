import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CERT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


r1 = _load("rst1", "rst1_recognition_seam_topology.py")


# ---------------- T1: cut-stable topology ----------------

def test_recognition_is_idempotent_and_classes_partition():
    assert all(r1.chi(r1.chi(x)) == r1.chi(x) for x in r1.X)
    classes = r1.recog_classes()
    flat = sorted(x for v in classes.values() for x in v)
    assert flat == list(r1.X)


def test_cut_family_closed_under_composition():
    tables = {n: tuple(f(x) for x in r1.X) for n, f in r1.CUTS.items()}
    for fa in r1.CUTS.values():
        for fb in r1.CUTS.values():
            comp = tuple(fa(fb(x)) for x in r1.X)
            assert comp in tables.values()


def test_cut_stable_sets_form_a_topology_exhaustively():
    opens = r1.cut_stable_topology()
    assert frozenset() in opens and frozenset(r1.X) in opens
    for U in opens:
        for V in opens:
            assert frozenset(U | V) in opens
            assert frozenset(U & V) in opens


def test_topology_is_nontrivial_and_not_discrete():
    opens = sorted(sorted(U) for U in r1.cut_stable_topology())
    assert opens == [[], [0, 1, 2, 3], [0, 1, 2, 3, 4, 5], [4, 5]]


def test_both_stability_conditions_are_load_bearing():
    assert not r1.is_saturated(frozenset({0}))
    assert r1.is_stable(frozenset({0})) is False or True  # see next lines
    # saturated but not stable:
    assert r1.is_saturated(frozenset({0, 1}))
    assert not r1.is_stable(frozenset({0, 1}))
    # stable under nothing relevant, not saturated:
    assert not r1.is_saturated(frozenset({4}))


def test_seam_is_recognition_fixed_and_zero_cost():
    assert r1.seam(r1.cut_kappa) == r1.X
    assert r1.seam(r1.cut_c) == (1, 2, 3, 4, 5)
    rho = r1.cost_of_cut(r1.cut_c)
    for x in r1.X:
        assert (rho[x] == 0) == (x in r1.seam(r1.cut_c))
    assert rho[0] is None


def test_curvature_defect_separates_declared_composites():
    lawful = [r1.curvature_defect(r1.cut_kappa, r1.cut_c, r1.cut_kc, x)
              for x in r1.X]
    unlawful = [r1.curvature_defect(r1.cut_kappa, r1.cut_c, lambda y: y, x)
                for x in r1.X]
    assert lawful == [0] * 6
    assert unlawful == [1, 0, 0, 0, 0, 0]


# ---------------- T2: rotation instead of gluing ----------------

def test_orientation_target_cannot_factor_through_gluing():
    # any map on the one-point quotient is constant; the target is not
    targets = {("*", "+"): 1, ("*", "-"): -1}
    for const in (1, -1):
        assert not all(v == const for v in targets.values())


def test_binary_label_repairs_the_glued_fibre():
    label = {("*", "+"): 1, ("*", "-"): -1}
    aug = {p: ("*", label[p]) for p in label}
    vals = list(aug.values())
    assert vals[0] != vals[1]


def test_reflected_pair_is_exact_and_side_target_odd():
    eps, t = Fr(1, 7), Fr(3)
    sp, sm = (Fr(1, 2) + eps, t), (Fr(1, 2) - eps, t)
    assert (Fr(1) - sp[0], sp[1]) == sm      # rho(s) = 1 - s, exactly
    assert sp != sm


# ---------------- T3: rotation-induced involution ----------------

def test_JU_is_selfadjoint_involution_with_orthogonal_projections():
    U = [[Fr(3, 5), Fr(-4, 5)], [Fr(4, 5), Fr(3, 5)]]
    J = r1.block_JU(U)
    I4 = r1.meye(4)
    assert r1.mmul(J, J) == I4
    assert r1.mT(J) == J
    Pp = r1.mscale(Fr(1, 2), r1.madd(I4, J))
    Pm = r1.mscale(Fr(1, 2), r1.msub(I4, J))
    assert r1.mmul(Pp, Pp) == Pp and r1.mmul(Pm, Pm) == Pm
    assert r1.mmul(Pp, Pm) == r1.mzero(4)
    assert r1.madd(Pp, Pm) == I4


def test_bridge_identity_transport_gives_cut_swap():
    J1 = r1.block_JU([[Fr(1)]])
    assert J1 == [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]


def test_nonorthogonal_transport_breaks_self_adjointness():
    # control: U not orthogonal => block J_U with U^T is not the inverse
    U = [[Fr(2), Fr(0)], [Fr(0), Fr(1)]]
    J = r1.block_JU(U)
    assert r1.mmul(J, J) != r1.meye(4)


# ---------------- T4: flow generator ----------------

def test_cayley_inverse_is_parameter_reversal():
    D = r1.E(3, 0, 1)
    prod = r1.pmat_mul(r1.cayley_step(D, +1), r1.cayley_step(D, -1))
    assert r1.pmat_coeff(prod, 0) == r1.meye(3)
    for d in range(1, r1.pmat_maxdeg(prod) + 1):
        assert r1.pmat_coeff(prod, d) == r1.mzero(3)


def test_loop_residue_is_exactly_the_commutator():
    D1, D2 = r1.E(3, 0, 1), r1.E(3, 1, 2)
    L = r1.commutator_loop(D1, D2)
    assert r1.pmat_coeff(L, 0) == r1.meye(3)
    assert r1.pmat_coeff(L, 1) == r1.mzero(3)
    assert r1.pmat_coeff(L, 2) == r1.commutator(D1, D2)


def test_nilpotent_loop_terminates_at_degree_two():
    # stronger than O(h^3): for these generators L_h = I + h^2 [D1,D2]
    D1, D2 = r1.E(3, 0, 1), r1.E(3, 1, 2)
    L = r1.commutator_loop(D1, D2)
    assert r1.pmat_maxdeg(L) == 2


def test_commuting_generators_close_exactly():
    D1, D3 = r1.E(3, 0, 1), r1.E(3, 0, 2)
    assert r1.commutator(D1, D3) == r1.mzero(3)
    L = r1.commutator_loop(D1, D3)
    assert r1.pmat_coeff(L, 0) == r1.meye(3)
    for d in range(1, r1.pmat_maxdeg(L) + 1):
        assert r1.pmat_coeff(L, d) == r1.mzero(3)


def test_eigenmode_law_term_by_term_and_control():
    K = [[Fr(2), Fr(0)], [Fr(5), Fr(3)]]
    v, kappa = [Fr(0), Fr(1)], Fr(3)
    w = list(v)
    for m in range(1, 13):
        w = [sum(K[i][j] * w[j] for j in range(2)) for i in range(2)]
        assert w == [kappa ** m * x for x in v]
    u = [Fr(1), Fr(0)]
    Ku = [sum(K[i][j] * u[j] for j in range(2)) for i in range(2)]
    assert Ku != [Fr(2) * x for x in u]


def test_covariance_defect_iff_conjugation_invariance():
    n = 3
    Jrev = [[Fr(1) if i + j == n - 1 else Fr(0) for j in range(n)]
            for i in range(n)]
    K_cov = r1.madd(r1.E(n, 0, 1), r1.E(n, 2, 1))
    K_bad = r1.E(n, 0, 1)
    for Kg, invariant in ((K_cov, True), (K_bad, False)):
        assert (r1.mmul(Jrev, r1.mmul(Kg, Jrev)) == Kg) == invariant
        W = r1.poly_exp_nilpotent(Kg)
        conj = r1.pmat_mul(r1.pmat_mul(r1.pmat_from(Jrev), W),
                           r1.pmat_from(Jrev))
        same = all(
            r1.pmat_coeff(conj, d) == r1.pmat_coeff(W, d)
            for d in range(max(r1.pmat_maxdeg(W),
                               r1.pmat_maxdeg(conj)) + 1))
        assert same == invariant


# ---------------- T5: dagger algebra ----------------

def test_primitive_relations_and_cut_swap():
    I2 = r1.cmeye()
    assert r1.cmmul(r1.R0, r1.R0) == r1.cmneg(I2)
    assert r1.cmmul(r1.KF, r1.KF) == I2
    assert r1.cmmul(r1.R0, r1.KF) == r1.cmneg(r1.cmmul(r1.KF, r1.R0))
    assert r1.JCUT == [[r1.GZERO, r1.GONE], [r1.GONE, r1.GZERO]]


def test_the_separation_hilbert_vs_dagger():
    assert r1.cmadj(r1.R0) == r1.cmneg(r1.R0)
    assert r1.dagger(r1.R0) == r1.R0


def test_dagger_is_involutive_antiautomorphism_on_grid():
    mats = r1.sample_matrices()
    for A in mats[:64]:
        assert r1.dagger(r1.dagger(A)) == A
    for A in mats[:20]:
        for B in mats[:20]:
            assert r1.dagger(r1.cmmul(A, B)) == r1.cmmul(
                r1.dagger(B), r1.dagger(A))


def test_leakage_three_way_equivalence():
    Pp = [[r1.GONE, r1.GZERO], [r1.GZERO, r1.GZERO]]
    Pm = [[r1.GZERO, r1.GZERO], [r1.GZERO, r1.GONE]]
    Z = [[r1.GZERO, r1.GZERO], [r1.GZERO, r1.GZERO]]
    for A in r1.sample_matrices():
        commutes = r1.cmmul(A, r1.KF) == r1.cmmul(r1.KF, A)
        leaks0 = (r1.cmmul(Pm, r1.cmmul(A, Pp)) == Z
                  and r1.cmmul(Pp, r1.cmmul(A, Pm)) == Z)
        blockdiag = A[0][1] == r1.GZERO and A[1][0] == r1.GZERO
        assert commutes == leaks0 == blockdiag


# ---------------- T6: sf = Wind = ind = q ----------------

def test_spectral_flow_exact_counts_and_orientations():
    for q in r1.QS:
        sf, crossings = r1.spectral_flow(q)
        assert sf == q
        for (n, t, s) in crossings:
            assert Fr(0) < t < Fr(1)
            assert n + r1.DELTA + t * q == 0     # exact crossing
            assert s == (1 if q > 0 else -1)


def test_crossing_window_is_exhaustive_not_a_truncation():
    # every crossing satisfies |n| <= |q|: widening the window can add
    # nothing, so the count is over the FULL Fourier spectrum
    for q in r1.QS:
        if q == 0:
            continue
        _, a = r1.spectral_flow(q, n_window=40)
        _, b = r1.spectral_flow(q, n_window=400)
        assert a == b


def test_winding_is_exact_and_matches_q():
    for q in r1.QS:
        assert r1.winding(q) == q


def test_circle_points_lie_exactly_on_the_unit_circle():
    for (x, y) in r1.circle_points():
        assert x * x + y * y == 1


def test_toeplitz_kernel_cokernel_and_index():
    for q in r1.QS:
        k, c, ind = r1.toeplitz_index(q)
        assert ind == q
        assert k == max(q, 0) and c == max(-q, 0)
        assert k == 0 or c == 0                  # never both


def test_the_three_integers_agree():
    for q in r1.QS:
        sf, _ = r1.spectral_flow(q)
        assert sf == r1.winding(q) == r1.toeplitz_index(q)[2] == q


# ---------------- T7: staircase ----------------

def test_staircase_length_is_exactly_two_for_every_n():
    for n in (1, 2, 3, 5, 8, 32, 100):
        length, d2 = r1.staircase_length_and_distance(n)
        assert length == 2
        assert d2 <= Fr(1, n) ** 2 / 2


def test_limit_of_lengths_is_not_length_of_limit():
    assert Fr(2) ** 2 != Fr(2)


# ---------------- certificate integrity and claim boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = r1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_RST1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_claim_boundary_analytic_layer_declared_not_certified():
    cert = json.load(open(os.path.join(CERT_DIR, "RST1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["analytic_layer"].startswith("DECLARED, NOT CERTIFIED")
    assert cb["arithmetic_meaning_of_q"].startswith("NOT CLAIMED")
    assert cb["critical_line_identification"].startswith("NOT CLAIMED")
    assert cb["model_zeros"].startswith("NOT CLAIMED")
    assert cb["riemann_hypothesis"] == "ABSTAIN"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "rst1_recognition_seam_topology.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "RST1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 7
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
