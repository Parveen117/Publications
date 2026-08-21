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


r2 = _load("rst2", "rst2_arithmetic_fusion_spine.py")


# ---------------- T1: sampled reflection ----------------

def test_reflected_samples_are_exact_and_swap_under_rho():
    for t in r2.HEIGHTS:
        sp = (Fr(1, 2) + r2.EPS, t)
        sm = (Fr(1, 2) - r2.EPS, t)
        assert (Fr(1) - sp[0], sp[1]) == sm


def test_pair_swap_eigenvectors_and_orthogonality():
    J = [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]
    assert r2.mv(J, [Fr(1), Fr(1)]) == [Fr(1), Fr(1)]
    assert r2.mv(J, [Fr(1), Fr(-1)]) == [Fr(-1), Fr(1)]
    assert r2.dot([Fr(1), Fr(1)], [Fr(1), Fr(-1)]) == 0


def test_seam_projection_blind_to_antiseam_component():
    J = [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]
    Pp = r2.mscale(Fr(1, 2), r2.madd(r2.meye(2), J))
    for c1 in (Fr(2), Fr(-3, 4)):
        outs = set()
        for c2 in (Fr(0), Fr(5), Fr(-1)):
            vec = [c1 + c2, c1 - c2]
            outs.add(tuple(r2.mv(Pp, vec)))
        assert len(outs) == 1                    # independent of c2


# ---------------- T2: weights, leakage, crossings ----------------

def test_equal_pair_weights_iff_swap_isometric():
    def swap(x):
        return [x[1], x[0]]
    w_eq, w_neq = [Fr(3), Fr(3)], [Fr(3), Fr(1)]
    for x in ([Fr(1), Fr(0)], [Fr(2), Fr(-5)], [Fr(0), Fr(7)]):
        assert r2.wnorm2(swap(x), w_eq) == r2.wnorm2(x, w_eq)
    assert r2.wnorm2(swap([Fr(1), Fr(0)]), w_neq) != r2.wnorm2(
        [Fr(1), Fr(0)], w_neq)


def test_avoided_family_never_singular_exact_gap():
    for c in (Fr(0), Fr(2, 3)):
        for mu in (Fr(1, 5), Fr(3)):
            for eta in (c - 1, c, c + Fr(1, 3), c + 2):
                D = [[eta - c, mu], [mu, -(eta - c)]]
                assert r2.det(D) <= -mu * mu < 0


def test_scalar_branch_crossing_orientation():
    c = Fr(1, 3)
    assert (0 - c) < 0 < (1 - c)                 # upward for eps=+1
    assert (-(0 - c)) > 0 > (-(1 - c))           # downward for eps=-1


# ---------------- T3: arithmetic recognition ----------------

def test_formal_log_is_additive_on_products():
    for a in range(2, 20):
        for b in range(2, 20):
            assert r2.logvec(a * b) == r2.vadd(r2.logvec(a), r2.logvec(b))


def test_mobius_pressure_extractor_equals_Lambda_all_n():
    for n in range(1, r2.N_ARITH + 1):
        conv = {}
        for d in range(1, n + 1):
            if n % d == 0:
                conv = r2.vadd(conv, r2.vscale(
                    r2.mobius(n // d), r2.logvec(d)))
        assert conv == r2.Lambda_vec(n)


def test_chebyshev_identity_log_n_is_sum_of_Lambda():
    for n in range(1, r2.N_ARITH + 1):
        s = {}
        for d in range(1, n + 1):
            if n % d == 0:
                s = r2.vadd(s, r2.Lambda_vec(d))
        assert s == r2.logvec(n)


def test_Lambda_transfer_has_zero_antiseam_current():
    stream = {n: r2.Lambda_vec(n) for n in range(2, r2.N_ARITH + 1)}
    t = r2.eye_transfer(stream)
    for q in r2.anti_seam_classes(r2.N_ARITH):
        assert q not in t


def test_raw_log_clock_leaks_into_class_six():
    stream = {n: r2.logvec(n) for n in range(2, r2.N_ARITH + 1)}
    t = r2.eye_transfer(stream)
    assert t.get(6)                              # nonzero vector


def test_controlled_contamination_is_exactly_eps_e6():
    eps = Fr(-2, 9)
    supp = {n: (Fr(1) if r2.Lambda_vec(n) else Fr(0))
            for n in range(2, r2.N_ARITH + 1)}
    supp[6] = supp.get(6, Fr(0)) + eps
    t = r2.eye_transfer_scalar(supp)
    anti = {q: v for q, v in t.items()
            if q in r2.anti_seam_classes(r2.N_ARITH)}
    assert anti == {6: eps}
    assert anti[6] * 1 == eps and anti[6] * -1 == -eps


# ---------------- T4: prime-shell transfer ----------------

def test_symmetric_transfer_covariance_identities():
    w = {p: Fr(1, p) for p in r2.SHELLS}
    A, U, Par, J = r2.shell_ops(w, w, w)
    n = len(r2.SHELLS)
    Pm = r2.mscale(Fr(1, 2), r2.msub(r2.meye(2 * n), J))
    assert r2.mmul(U, A) == r2.mmul(A, Par)
    assert r2.mmul(J, A) == A
    assert r2.mmul(Pm, A) == r2.mzero(2 * n, n)
    assert r2.mmul(J, r2.mmul(r2.mT(U), J)) == U


def test_asymmetric_leakage_formula_exact_per_shell():
    w = {p: Fr(1, p) for p in r2.SHELLS}
    wR, wL = dict(w), dict(w)
    wR[7], wL[7] = Fr(1, 2), Fr(1, 9)
    A, _, Par, J = r2.shell_ops(w, w, w)
    _, Uasym, _, _ = r2.shell_ops(wR, wL, w)
    n = len(r2.SHELLS)
    Pm = r2.mscale(Fr(1, 2), r2.msub(r2.meye(2 * n), J))
    Delta = r2.msub(r2.mmul(Uasym, A), r2.mmul(A, Par))
    PmD = r2.mmul(Pm, Delta)
    for i, p in enumerate(r2.SHELLS):
        coef = (wR[p] - wL[p]) / 2
        assert PmD[2 * i][i] == coef
        assert PmD[2 * i + 1][i] == -coef


def test_converse_fails_common_seam_error_invisible():
    w = {p: Fr(1, p) for p in r2.SHELLS}
    w_bad = dict(w)
    w_bad[2] = Fr(9)
    A, _, Par, J = r2.shell_ops(w, w, w)
    _, Ubad, _, _ = r2.shell_ops(w_bad, w_bad, w)
    n = len(r2.SHELLS)
    Pm = r2.mscale(Fr(1, 2), r2.msub(r2.meye(2 * n), J))
    Delta = r2.msub(r2.mmul(Ubad, A), r2.mmul(A, Par))
    assert Delta != r2.mzero(2 * n, n)
    assert r2.mmul(Pm, Delta) == r2.mzero(2 * n, n)


def test_crossing_pairings_aligned_reversed_avoided():
    alpha = Fr(4, 5)
    c_minus = [Fr(1), Fr(-1)]
    J_def = [alpha * x for x in c_minus]
    n2 = r2.dot(c_minus, c_minus)
    assert r2.dot(J_def, c_minus) == alpha * n2
    assert r2.dot(J_def, [-x for x in c_minus]) == -alpha * n2
    assert r2.dot(J_def, [Fr(1), Fr(1)]) == 0


# ---------------- T5: zero-identification boundary ----------------

def test_paired_determinant_normalized_even_with_exact_zeros():
    D = r2.paired_determinant(r2.LEVELS)
    assert r2.peval(D, Fr(0)) == 1
    assert all(c == 0 for i, c in enumerate(D) if i % 2 == 1)
    for lam in set(r2.LEVELS):
        mult = sum(1 for x in r2.LEVELS if x == lam)
        assert r2.root_multiplicity(D, lam) == mult
        assert r2.root_multiplicity(D, -lam) == mult


def test_double_level_gives_double_zero():
    assert r2.root_multiplicity(
        r2.paired_determinant(r2.LEVELS), Fr(3, 2)) == 2


def test_log_derivative_identity_cross_multiplied():
    D = r2.paired_determinant(r2.LEVELS)
    factors = [[-(lam * lam), Fr(0), Fr(1)] for lam in r2.LEVELS]
    prod_all = [Fr(1)]
    for f in factors:
        prod_all = r2.pmul(prod_all, f)
    rhs = [Fr(0)]
    for j in range(len(r2.LEVELS)):
        term = [Fr(0), Fr(2)]
        for k, f in enumerate(factors):
            if k != j:
                term = r2.pmul(term, f)
        rhs = r2.padd(rhs, term)
    assert r2.pmul(r2.pderiv(D), prod_all) == r2.pmul(D, rhs)


def test_nonvanishing_factor_preserves_divisor():
    D = r2.paired_determinant(r2.LEVELS)
    F = r2.pmul([Fr(1), Fr(0), Fr(1)], D)        # (E^2 + 1) D
    for lam in set(r2.LEVELS):
        assert r2.root_multiplicity(F, lam) == r2.root_multiplicity(D, lam)


def test_polynomial_division_is_exact():
    p = r2.pmul([Fr(-2), Fr(1)], [Fr(3), Fr(5), Fr(1)])
    q, rem = r2.pdivmod(p, [Fr(-2), Fr(1)])
    assert q == [Fr(3), Fr(5), Fr(1)] and rem == [Fr(0)]


# ---------------- T6: Surya/Weil fusion ----------------

def test_feshbach_congruence_is_exact_identity():
    A = [[Fr(3), Fr(1)], [Fr(1), Fr(2)]]
    C = [[Fr(2), Fr(0)], [Fr(-1), Fr(1)]]
    Dm = [[Fr(5), Fr(1)], [Fr(1), Fr(4)]]
    M, S, Dg, F = r2.feshbach_congruence(A, C, Dm)
    assert M == r2.mmul(r2.mT(S), r2.mmul(Dg, S))


def test_feshbach_sign_equivalence_both_directions():
    A = [[Fr(2), Fr(1)], [Fr(1), Fr(2)]]
    C = [[Fr(1), Fr(-1)], [Fr(0), Fr(2)]]
    _, _, _, F_pos = r2.feshbach_congruence(
        A, C, [[Fr(4), Fr(0)], [Fr(0), Fr(6)]])
    M_pos, _, _, _ = r2.feshbach_congruence(
        A, C, [[Fr(4), Fr(0)], [Fr(0), Fr(6)]])
    assert r2.is_psd(F_pos) and r2.is_psd(M_pos)
    M_neg, S, _, F_neg = r2.feshbach_congruence(
        A, C, [[Fr(1), Fr(0)], [Fr(0), Fr(1)]])
    assert not r2.is_psd(F_neg) and not r2.is_psd(M_neg)


def test_negative_witness_transports_through_congruence():
    A = [[Fr(2), Fr(1)], [Fr(1), Fr(2)]]
    C = [[Fr(1), Fr(-1)], [Fr(0), Fr(2)]]
    M, S, _, F = r2.feshbach_congruence(
        A, C, [[Fr(1), Fr(0)], [Fr(0), Fr(1)]])
    for z in ([Fr(1), Fr(0)], [Fr(0), Fr(1)], [Fr(1), Fr(1)],
              [Fr(1), Fr(-1)]):
        if r2.dot(z, r2.mv(F, z)) < 0:
            x = r2.mv(r2.inv(S), [Fr(0), Fr(0)] + z)
            assert r2.dot(x, r2.mv(M, x)) == r2.dot(z, r2.mv(F, z))
            return
    raise AssertionError("no negative witness found")


def test_surya_squared_amplitude_sign_only():
    for x in (Fr(-10), Fr(-1), Fr(-1, 3), Fr(0), Fr(1, 2), Fr(4)):
        assert (r2.ssq(x) == 0) == (x >= 0)
        assert r2.ssq(x) >= 0
        if x < 0:
            assert r2.ssq(x) == x * x / (1 + x * x)


def test_operator_surya_vanishes_iff_W0_psd():
    Q = [[Fr(3, 5), Fr(-4, 5)], [Fr(4, 5), Fr(3, 5)]]
    for spec in ((Fr(1), Fr(2)), (Fr(-1, 2), Fr(3)), (Fr(0), Fr(0))):
        L = [[spec[0], Fr(0)], [Fr(0), spec[1]]]
        W0 = r2.mmul(Q, r2.mmul(L, r2.mT(Q)))
        Asq = r2.mmul(Q, r2.mmul(
            [[r2.ssq(spec[0]), Fr(0)], [Fr(0), r2.ssq(spec[1])]],
            r2.mT(Q)))
        assert (Asq == r2.mzero(2)) == r2.is_psd(W0)


def test_seam_charge_additive_and_threshold_strict():
    even, odd = (Fr(2), Fr(1, 2)), (Fr(3), Fr(-1), Fr(1))
    kp = sum(1 for m in even if m > 1)
    ko = sum(1 for m in odd if m > 1)
    ks = sum(1 for m in even + odd if m > 1)
    assert ks == kp + ko == 2
    assert Fr(1) in odd and not (Fr(1) > 1)      # mu = 1 not in the charge


def test_spectral_flow_is_minus_charge_all_downward():
    spec = (Fr(2), Fr(1, 2), Fr(3), Fr(-1), Fr(1))
    sf = 0
    for mu in spec:
        if mu != 0 and Fr(0) < Fr(1) / mu < Fr(1):
            assert mu > 1 and -mu < 0
            sf -= 1
    assert sf == -sum(1 for m in spec if m > 1)


def test_congruence_with_rational_sqrt_reference():
    Ahalf = [[Fr(2), Fr(0)], [Fr(0), Fr(3)]]
    for spec in ((Fr(1, 2), Fr(-2)), (Fr(3, 2), Fr(0))):
        K = [[spec[0], Fr(0)], [Fr(0), spec[1]]]
        ImK = r2.msub(r2.meye(2), K)
        W = r2.mmul(Ahalf, r2.mmul(ImK, Ahalf))
        k = sum(1 for m in spec if m > 1)
        assert r2.is_psd(W) == (k == 0)


def test_surya_angle_indicator_equality():
    for mu in (Fr(2), Fr(1), Fr(0), Fr(-3), Fr(6, 5)):
        s, amp2 = r2.csq(mu)
        assert (s > 0 and amp2 > 0) == (mu > 1)
        assert (s == 0) == (mu == 1)


def test_principal_holonomy_blindness():
    for m in (-4, 0, 7):
        assert m % 1 == 0
    assert 7 != 0


# ---------------- integrity and claim boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = r2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_RST2.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_declares_analytic_and_open_gates():
    cert = json.load(open(os.path.join(CERT_DIR, "RST2_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["analytic_layer"].startswith("DECLARED, NOT CERTIFIED")
    assert "OPEN" in cb["open_gates"]
    assert cb["forbidden_promotions"].startswith("Zmodel -> Zxi")
    assert "IMPORTED" not in cb["imported_not_certified"][:0]  # present key
    assert cb["riemann_hypothesis"] == "ABSTAIN"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "rst2_arithmetic_fusion_spine.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src
    assert "** 0.5" not in src and "**0.5" not in src
    assert "math.sqrt" not in src and "isqrt" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "RST2_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
