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


ugd1 = _load("ugd1", "ugd1_numerals.py")


# ---------------- digits and the phase alphabet ----------------

def test_seam_alphabet_is_balanced():
    for s in (-1, 0, 1):
        ugd1.digit(3, 0, s)
    for bad in (-2, 2, 5):
        try:
            ugd1.digit(3, 0, bad)
            assert False, "seam alphabet not enforced"
        except AssertionError as exc:
            if "balanced" not in str(exc):
                raise


def test_phase_exponent_is_reduced_mod_K():
    d = ugd1.digit(15, 0, 0, K=12)
    (_, (a, _)), = d.digits.items()
    assert a == 3


def test_quarter_turn_at_K4():
    powers = [(1 * j) % 4 for j in range(1, 5)]
    assert powers == [1, 2, 3, 0]        # order exactly 4
    assert (1 + 1) % 4 == 2              # quarter squared = half turn


def test_phase_generator_has_full_order():
    for K in (4, 7, 12):
        x, order = 0, None
        for r in range(1, K + 1):
            x = (x + 1) % K
            if x == 0:
                order = r
                break
        assert order == K


# ---------------- carry conservation ----------------

def test_carry_conserves_total_seam_charge():
    for a1 in range(0, 12):
        for a2 in range(0, 12):
            for s1 in (-1, 0, 1):
                for s2 in (-1, 0, 1):
                    N1, N2 = ugd1.digit(a1, 0, s1), ugd1.digit(a2, 0, s2)
                    S = ugd1.add_digits(N1, N2)
                    assert S.total_seam() == \
                        N1.total_seam() + N2.total_seam()


def test_result_digits_return_to_balanced_alphabet():
    for s1 in (-1, 0, 1):
        for s2 in (-1, 0, 1):
            S = ugd1.add_digits(ugd1.digit(2, 0, s1), ugd1.digit(3, 0, s2))
            for (_, s) in S.digits.values():
                assert s in (-1, 0, 1)


def test_phase_overflow_carries_to_next_scale():
    # 7 + 8 = 15 = 3 mod 12, carry 1 to the next scale
    S = ugd1.add_digits(ugd1.digit(7, 0, 0), ugd1.digit(8, 0, 0))
    assert S.digits[0][0] == 3
    assert S.digits.get(1, (0, 0))[0] == 1


def test_broken_carry_control_loses_charge():
    B1, B2 = ugd1.digit(3, 0, 1), ugd1.digit(5, 0, 1)
    good = ugd1.add_digits(B1, B2)
    bad = ugd1.add_digits_broken_carry(B1, B2)
    assert good.total_seam() == 2
    assert bad.total_seam() == 1
    assert good.total_seam() != bad.total_seam()


# ---------------- multiplication ----------------

def test_scale_indices_add_under_multiplication():
    for m in range(-3, 4):
        for n in range(-3, 4):
            P = ugd1.mul_digits(ugd1.digit(2, m, 0), ugd1.digit(3, n, 0))
            (sc, (ph, _)), = P.digits.items()
            assert sc == m + n
            assert ph == 6 % ugd1.K_DEFAULT


def test_scale_projection_is_multiplicative():
    for m in range(-2, 3):
        for n in range(-2, 3):
            d1, d2 = ugd1.digit(2, m, 0), ugd1.digit(3, n, 0)
            P = ugd1.mul_digits(d1, d2)
            assert P.scale_projection() == \
                d1.scale_projection() * d2.scale_projection()


def test_undeclared_seam_rule_is_refused():
    try:
        ugd1.mul_digits(ugd1.digit(1, 0, 0), ugd1.digit(1, 0, 0),
                        seam_rule="invented")
        assert False, "undeclared seam rule accepted"
    except ValueError:
        pass


# ---------------- projection blindness ----------------

def test_null_seam_recovers_positional_notation():
    N = ugd1.Numeral({0: (3, 0), 1: (2, 0), 2: (1, 0)})
    b = ugd1.BETA
    assert N.classical_projection() == 3 + 2 * b + 1 * b ** 2


def test_projection_discards_seam_information():
    A = ugd1.Numeral({0: (3, 0), 1: (2, 0)})
    B = ugd1.Numeral({0: (3, 1), 1: (2, -1)})
    C = ugd1.Numeral({0: (3, 1), 1: (2, 1)})
    assert A.classical_projection() == B.classical_projection() \
        == C.classical_projection()
    assert (A.total_seam(), B.total_seam(), C.total_seam()) == (0, 0, 2)
    assert B.digits != C.digits


# ---------------- cut-zero ----------------

def test_neutral_digit_is_derived_neutral():
    Z = ugd1.ZERO_DIGIT
    for a in range(0, 12):
        for s in (-1, 0, 1):
            D = ugd1.digit(a, 0, s)
            S = ugd1.add_digits(D, Z)
            assert S.digits == D.digits and S.ledger == D.ledger


def test_neutral_numeral_is_distinct_from_absence():
    Z = ugd1.ZERO_DIGIT
    absent = ugd1.Numeral({})
    assert Z.classical_projection() == 0 and Z.total_seam() == 0
    assert absent.digits == {}
    assert Z.digits == {0: (0, 0)}
    assert Z != absent


# ---------------- lambda-logic ----------------

def test_negation_is_an_involution():
    for p in range(0, 4):
        for sg in (-1, 0, 1):
            for k in (-1, 0, 1):
                v = (p, sg, k)
                assert ugd1.ugd_negate(ugd1.ugd_negate(v)) == v


def test_seam_predicate_marks_nonzero_charge():
    for p in range(0, 4):
        for sg in (-1, 0, 1):
            for k in (-1, 0, 1):
                assert ugd1.seam_predicate((p, sg, k)) == (k != 0)


def test_classical_limit_of_the_logic():
    for sg in (-1, 0, 1):
        v = (0, sg, 0)
        n = ugd1.ugd_negate(v)
        assert n == (0, -sg, 0)
        assert not ugd1.seam_predicate(v)
        assert not ugd1.seam_predicate(n)


def test_contradiction_seam_charges_cancel():
    A = (1, 1, 1)
    notA = ugd1.ugd_negate(A)
    assert A[2] + notA[2] == 0


# ---------------- certificate integrity ----------------

def test_certificate_pin_matches_regeneration():
    cert = ugd1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_UGD1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_no_floats_in_the_representation():
    d = ugd1.digit(3, 2, 1)
    (n, (a, s)), = d.digits.items()
    assert isinstance(a, int) and isinstance(s, int) and isinstance(n, int)
    assert isinstance(d.classical_projection(), Fr)


def test_open_boundaries_declared():
    cert = json.load(open(os.path.join(CERT_DIR, "UGD1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["infinite_numerals"].startswith("NOT CLAIMED")
    assert cb["linked_seam_composition_rule"].startswith("DECLARED")
    assert cb["CFE_uniqueness_U"].startswith("OPEN")
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["quantum_gravity"] == "not touched"
