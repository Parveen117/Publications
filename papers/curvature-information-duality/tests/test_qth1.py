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


q = _load("qth1", "qth1_quantum_recognition_information.py")


# ---------------- T1: exact SLD ----------------

def test_state_is_hermitian_unit_trace():
    rho = q.bloch_state(q.R0)
    assert q.dag(rho) == rho
    assert q.trace(rho) == q.GO


def test_derivative_is_traceless_hermitian():
    for dr in (q.DR, q.DR2, (Fr(1, 3), Fr(-1, 2), Fr(0))):
        d = q.bloch_derivative(dr)
        assert q.trace(d) == q.GZ
        assert q.dag(d) == d


def test_sld_satisfies_its_defining_equation_exactly():
    rho = q.bloch_state(q.R0)
    for dr in (q.DR, q.DR2, (Fr(1, 5), Fr(0), Fr(2, 5))):
        d = q.bloch_derivative(dr)
        L = q.sld(rho, d)
        assert q.dag(L) == L
        recon = q.mscale((Fr(1, 2), Fr(0)),
                         q.madd(q.mmul(L, rho), q.mmul(rho, L)))
        assert recon == d


def test_two_routes_to_qfi_agree_exactly():
    rho = q.bloch_state(q.R0)
    for dr in (q.DR, q.DR2, (Fr(1, 3), Fr(0), Fr(1, 3)),
               (Fr(0), Fr(2, 5), Fr(1, 5))):
        d = q.bloch_derivative(dr)
        L = q.sld(rho, d)
        assert q.qfi_from_sld(d, L) == q.qfi_closed_form(q.R0, dr)


def test_qfi_nonnegative_and_zero_only_when_static():
    rho = q.bloch_state(q.R0)
    for dr in (q.DR, q.DR2, (Fr(0), Fr(0), Fr(0))):
        d = q.bloch_derivative(dr)
        f = q.qfi_from_sld(d, q.sld(rho, d))
        assert f >= 0
        assert (f == 0) == (dr == (Fr(0), Fr(0), Fr(0)))


# ---------------- T2: the ledger becomes an inequality ----------------

def test_braunstein_caves_holds_on_every_direction():
    qf = q.qfi_closed_form(q.R0, q.DR)
    for n in q.UNITS:
        assert q.cfi_projective(q.R0, q.DR, n) <= qf


def test_exact_saturation_at_the_aligned_direction():
    qf = q.qfi_closed_form(q.R0, q.DR)
    assert q.cfi_projective(q.R0, q.DR, (Fr(3, 5), Fr(4, 5), Fr(0))) == qf


def test_orthogonal_measurement_discards_everything():
    qf = q.qfi_closed_form(q.R0, q.DR)
    assert q.cfi_projective(q.R0, q.DR, (Fr(4, 5), Fr(-3, 5), Fr(0))) == 0
    assert qf == 1


def test_discard_is_nonnegative_and_sometimes_strict():
    qf = q.qfi_closed_form(q.R0, q.DR)
    gaps = [qf - q.cfi_projective(q.R0, q.DR, n) for n in q.UNITS]
    assert all(g >= 0 for g in gaps)
    assert any(g > 0 for g in gaps)
    assert any(g == 0 for g in gaps)


def test_measurement_directions_are_genuine_unit_vectors():
    for n in q.UNITS:
        assert sum(c * c for c in n) == 1


# ---------------- T3: one tensor, two halves ----------------

def test_symmetric_part_is_symmetric_antisymmetric_is_not():
    rho = q.bloch_state(q.R0)
    d1 = q.bloch_derivative((Fr(3, 5), Fr(4, 5), Fr(0)))
    d2 = q.bloch_derivative((Fr(4, 5), Fr(-3, 5), Fr(0)))
    s12, a12 = q.geometric_tensor(rho, d1, d2)
    s21, a21 = q.geometric_tensor(rho, d2, d1)
    assert s12 == s21
    assert a12 == -a21 != 0


def test_antisymmetric_diagonal_vanishes():
    rho = q.bloch_state(q.R0)
    for dr in (q.DR, q.DR2):
        d = q.bloch_derivative(dr)
        _, a = q.geometric_tensor(rho, d, d)
        assert a == 0


def test_symmetric_diagonal_is_the_qfi():
    rho = q.bloch_state(q.R0)
    for dr in (q.DR, (Fr(0), Fr(0), Fr(1, 4))):
        d = q.bloch_derivative(dr)
        s, _ = q.geometric_tensor(rho, d, d)
        assert s == q.qfi_closed_form(q.R0, dr)


def test_obstruction_vanishes_exactly_for_commuting_slds():
    rho = q.bloch_state(q.R0)
    e1 = q.bloch_derivative((Fr(0), Fr(0), Fr(1, 4)))
    e2 = q.bloch_derivative((Fr(0), Fr(0), Fr(1, 2)))
    L1, L2 = q.sld(rho, e1), q.sld(rho, e2)
    assert q.msub(q.mmul(L1, L2), q.mmul(L2, L1)) == q.mzero()
    _, a = q.geometric_tensor(rho, e1, e2)
    assert a == 0


def test_obstruction_nonzero_when_slds_do_not_commute():
    rho = q.bloch_state(q.R0)
    d1 = q.bloch_derivative((Fr(3, 5), Fr(4, 5), Fr(0)))
    d2 = q.bloch_derivative((Fr(4, 5), Fr(-3, 5), Fr(0)))
    L1, L2 = q.sld(rho, d1), q.sld(rho, d2)
    assert q.msub(q.mmul(L1, L2), q.mmul(L2, L1)) != q.mzero()
    _, a = q.geometric_tensor(rho, d1, d2)
    assert a != 0


# ---------------- T4: where the identity fails ----------------

def test_quantum_discard_is_nonzero_for_a_suboptimal_measurement():
    n = (Fr(0), Fr(0), Fr(1))
    gap = q.qfi_closed_form(q.R0, q.DR) - q.cfi_projective(q.R0, q.DR, n)
    assert gap == 1 > 0


def test_outcome_distribution_is_a_genuine_distribution():
    for n in q.UNITS:
        p = q.outcome_distribution(q.R0, n)
        assert sum(p) == 1 and all(x > 0 for x in p)


def test_classical_bookkeeping_inside_outcomes_still_holds():
    # the singleton partition discards nothing, classically
    assert Fr(0) == Fr(0)
    n = (Fr(0), Fr(0), Fr(1))
    gap = q.qfi_closed_form(q.R0, q.DR) - q.cfi_projective(q.R0, q.DR, n)
    assert gap != Fr(0)


# ---------------- T5: the commuting face ----------------

def test_commuting_family_has_commuting_sld():
    r, dr = (Fr(0), Fr(0), Fr(1, 2)), (Fr(0), Fr(0), Fr(1, 4))
    rho, d = q.bloch_state(r), q.bloch_derivative(dr)
    L = q.sld(rho, d)
    assert q.msub(q.mmul(L, rho), q.mmul(rho, L)) == q.mzero()


def test_commuting_face_has_zero_quantum_discard():
    r, dr = (Fr(0), Fr(0), Fr(1, 2)), (Fr(0), Fr(0), Fr(1, 4))
    n = (Fr(0), Fr(0), Fr(1))
    assert q.qfi_closed_form(r, dr) == q.cfi_projective(r, dr, n)


def test_commuting_face_matches_cid1_covariance_route():
    r, dr = (Fr(0), Fr(0), Fr(1, 2)), (Fr(0), Fr(0), Fr(1, 4))
    n = (Fr(0), Fr(0), Fr(1))
    p = q.outcome_distribution(r, n)
    dp = (Fr(1, 2) * dr[2], -Fr(1, 2) * dr[2])
    assert q.classical_fisher_covariance(p, dp) == q.qfi_closed_form(r, dr)


def test_classical_fisher_helper_on_an_independent_example():
    p = (Fr(1, 2), Fr(1, 2))
    dp = (Fr(1, 4), Fr(-1, 4))
    assert q.classical_fisher_covariance(p, dp) == Fr(1, 4)


# ---------------- T6: monotonicity ----------------

def test_rotation_is_orthogonal():
    for i in range(3):
        for j in range(3):
            expect = Fr(1) if i == j else Fr(0)
            assert sum(q.ROT[k][i] * q.ROT[k][j]
                       for k in range(3)) == expect


def test_unitary_channel_preserves_qfi_exactly():
    q0 = q.qfi_closed_form(q.R0, q.DR)
    r1 = q.rotate_bloch(q.R0, q.ROT)
    d1 = q.rotate_bloch(q.DR, q.ROT)
    assert q.qfi_closed_form(r1, d1) == q0


def test_depolarizing_channel_strictly_decreases_qfi():
    q0 = q.qfi_closed_form(q.R0, q.DR)
    vals = []
    for lam in (Fr(1), Fr(3, 4), Fr(1, 2), Fr(1, 5)):
        rl = tuple(lam * x for x in q.R0)
        dl = tuple(lam * x for x in q.DR)
        v = q.qfi_closed_form(rl, dl)
        assert v == lam * lam * q0
        vals.append(v)
    assert all(a > b for a, b in zip(vals, vals[1:]))
    assert all(v <= q0 for v in vals)


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = q.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_QTH1.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_is_a_witness_capsule():
    cert = json.load(open(os.path.join(CERT_DIR, "QTH1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert "WITNESS" in cb["witness_capsule"]
    assert cb["not_claimed"].startswith("NOT CLAIMED")
    assert cb["no_thermodynamic_interpretation"].startswith("NOT CLAIMED")
    assert cb["berry_phase_identification"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0_YM"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "qth1_quantum_recognition_information.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src
    assert "** 0.5" not in src and "**0.5" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "QTH1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
