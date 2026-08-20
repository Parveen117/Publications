import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load():
    spec = importlib.util.spec_from_file_location(
        "lam1", os.path.join(CERT_DIR, "lam1_seam_interface.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules["lam1"] = mod
    spec.loader.exec_module(mod)
    return mod


lam1 = _load()


def test_cyclotomic_phi12():
    # Phi_12(x) = x^4 - x^2 + 1
    assert lam1.cyclotomic(12) == [1, 0, -1, 0, 1]


def test_ramanujan_known_values():
    # c_1(n) = 1; prime p: c_p(n) = p-1 if p|n else -1
    assert lam1.ramanujan_direct(1, 7) == (1, True)
    assert lam1.ramanujan_direct(5, 10) == (4, True)
    assert lam1.ramanujan_direct(5, 3) == (-1, True)
    assert lam1.ramanujan_moebius(5, 10) == 4
    assert lam1.ramanujan_moebius(5, 3) == -1


def test_two_route_agreement_sample():
    for c in (6, 8, 9, 12, 30, 49):
        for n in (1, 2, 6, 12, 30):
            val, is_int = lam1.ramanujan_direct(c, n)
            assert is_int
            assert val == lam1.ramanujan_moebius(c, n)


def test_control_coprimality_tamper_separates():
    val, is_int = lam1.ramanujan_direct(4, 2, require_coprime=False)
    assert is_int and val == 0
    assert lam1.ramanujan_moebius(4, 2) == -2
    assert val != lam1.ramanujan_moebius(4, 2)


def test_omega_additivity_and_euler_marker():
    for m, n in [(2, 3), (4, 9), (12, 25), (30, 7)]:
        assert (lam1.omega_total(m * n)
                == lam1.omega_total(m) + lam1.omega_total(n))
    table = lam1.euler_marker_witness(60)
    assert set(table.keys()) == set(range(1, 61))
    assert table[12] == (3, 1)   # 12 = 2^2 * 3, Omega = 3, multiplicity 1
    assert table[59] == (1, 1)


def test_control_fake_winding_detected():
    fake = lambda n: lam1.omega_total(n) + (1 if n > 1 else 0)
    assert fake(6) != fake(2) + fake(3)


def test_theta_seam_overlap_and_budget():
    th_half = lam1.theta_iv(Fr(1, 2))
    th_two = lam1.theta_iv(Fr(2))
    s2 = lam1.sqrt_brackets(2)
    rhs = lam1.iv_mul_pos(s2, th_two)
    assert lam1.iv_overlap(th_half, rhs)
    budget = Fr(1, 10 ** 20)
    assert lam1.iv_width(th_half) < budget
    assert lam1.iv_width(rhs) < budget


def test_control_sqrt3_tamper_separates():
    th_half = lam1.theta_iv(Fr(1, 2))
    th_two = lam1.theta_iv(Fr(2))
    s3 = lam1.sqrt_brackets(3)
    rhs3 = lam1.iv_mul_pos(s3, th_two)
    assert not lam1.iv_overlap(th_half, rhs3)


def test_s_move_bijection_and_zero_net_weight():
    region = lam1.primitive_region(25)
    image = {lam1.s_move(p) for p in region}
    assert image == region
    w_c, w_d = lam1.weight_sums(region)
    assert w_c == w_d


def test_control_asymmetric_region_separates():
    from math import gcd
    C = 25
    asym = {(c, d) for c in range(1, C + 1) for d in range(1, 2 * C + 1)
            if gcd(c, d) == 1}
    aw_c = sum(lam1.omega_total(c) for (c, d) in asym)
    aw_d = sum(lam1.omega_total(d) for (c, d) in asym)
    assert aw_c != aw_d


def test_directed_rounding_soundness():
    x = Fr(1, 10 ** 400)   # e^{-900}-scale magnitudes must survive
    lo, hi = lam1.round_out(x, x)
    assert 0 < lo <= x <= hi


def test_certificate_pin_matches_regeneration():
    cert = lam1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_LAM1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_claim_boundary_present_and_open():
    with open(os.path.join(CERT_DIR, "LAM1_RESULT.json")) as f:
        cert = json.load(f)
    cb = cert["claim_boundary"]
    for key in ("N1_native_continuation", "N2_native_functional_equation",
                "N3_identification", "K0_L0_RH", "YM_continuum_gates"):
        assert cb[key] == "OPEN"
    assert cb["zeros_of_any_L_function"] == "no claims"
