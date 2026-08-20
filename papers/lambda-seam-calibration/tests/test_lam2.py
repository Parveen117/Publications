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


lam2 = _load("lam2", "lam2_derived_phase.py")
lam1 = lam2.lam1


def test_cyclotomic_ring_basics():
    # zeta_4^2 = -1 as an exact ring statement
    v = lam2.vec_new(4)
    lam2.vec_add_monomial(v, 1)          # zeta_4
    sq = lam2.vec_mul(v, v)
    assert lam2.elem_as_int(sq) == -1


def test_gauss_magnitude_prime_cases():
    for q, j in [(3, 1), (5, 1), (5, 2), (7, 3), (13, 5)]:
        ch = lam2.CyclicCharacter(q, j)
        assert ch.is_primitive()
        tau = ch.gauss_sum_vec()
        assert lam2.elem_as_int(
            lam2.vec_mul(tau, lam2.vec_conj(tau))) == q


def test_gauss_quarter_turn_q4():
    ch = lam2.CyclicCharacter(4, 1)
    tau = ch.gauss_sum_vec()
    assert lam2.elem_as_int(lam2.vec_mul(tau, tau)) == -4


def test_imprimitive_mod9_magnitude_zero():
    ch = lam2.CyclicCharacter(9, 3)
    assert not ch.is_primitive()
    tau = ch.gauss_sum_vec()
    assert lam2.elem_as_int(lam2.vec_mul(tau, lam2.vec_conj(tau))) == 0


def test_derivation_step_q5_all_residues():
    ch = lam2.CyclicCharacter(5, 1)
    tau_bar = ch.gauss_sum_vec(conjugate=True)
    for n in range(5):
        rhs = lam2.vec_new(ch.L)
        for a in lam2.unit_group(5):
            lam2.vec_add_monomial(
                rhs, ch.conj_value_exponent(a) + a * n * ch.step_add)
        e_n = ch.value_exponent(n)
        lhs = (lam2.vec_new(ch.L) if e_n is None
               else lam2.vec_shift(tau_bar, e_n))
        assert lam2.elem_is_zero(lam2.vec_sub(rhs, lhs))


def test_derivation_step_fails_imprimitive():
    ch = lam2.CyclicCharacter(9, 3)
    tau_bar = ch.gauss_sum_vec(conjugate=True)
    n = 3
    rhs = lam2.vec_new(ch.L)
    for a in lam2.unit_group(9):
        lam2.vec_add_monomial(
            rhs, ch.conj_value_exponent(a) + a * n * ch.step_add)
    e_n = ch.value_exponent(n)
    lhs = (lam2.vec_new(ch.L) if e_n is None
           else lam2.vec_shift(tau_bar, e_n))
    assert not lam2.elem_is_zero(lam2.vec_sub(rhs, lhs))


def test_twisted_theta_seam_overlap_and_budget():
    th_half = lam2.twisted_theta_iv(Fr(1, 2))
    th_two = lam2.twisted_theta_iv(Fr(2))
    rhs = lam1.iv_mul_pos(lam1.sqrt_brackets(8), th_two)
    assert lam1.iv_overlap(th_half, rhs)
    assert lam1.iv_width(th_half) < Fr(1, 10 ** 20)
    assert lam1.iv_width(rhs) < Fr(1, 10 ** 20)


def test_control_multiplier_tamper_separates():
    th_half = lam2.twisted_theta_iv(Fr(1, 2))
    th_two = lam2.twisted_theta_iv(Fr(2))
    rhs_bad = lam1.iv_mul_pos(lam1.sqrt_brackets(2), th_two)
    assert not lam1.iv_overlap(th_half, rhs_bad)


def test_no_conductor_witnesses():
    for q in (2, 7, 12, 30, 60):
        found = False
        for n in range(2, 4000):
            if lam1.omega_total(n) != lam1.omega_total(n + q):
                found = True
                break
        assert found


def test_certificate_pin_matches_regeneration():
    cert = lam2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_LAM2.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_claim_boundary_open():
    with open(os.path.join(CERT_DIR, "LAM2_RESULT.json")) as f:
        cert = json.load(f)
    cb = cert["claim_boundary"]
    assert cb["N1_native_continuation"] == "OPEN"
    assert cb["N2_native_functional_equation"].startswith("OPEN")
    assert cb["N3_identification"] == "OPEN"
    assert cb["K0_L0_RH"] == "OPEN"
    assert cb["YM_continuum_gates"] == "OPEN"
    assert cb["zeros_of_any_L_function"] == "no claims"
