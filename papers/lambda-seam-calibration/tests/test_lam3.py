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


lam3 = _load("lam3", "lam3_crossing_derivation.py")
lam1 = lam3.lam1


def _ref_bracket(digits30, slack=Fr(1, 10 ** 28)):
    v = Fr(digits30)
    return (v - slack, v + slack)


# External 30-digit reference anchors (classical witnesses, CIRC-1:
# verification channel, never rederived by the capsule under test)
REF_GAMMA = "0.577215664901532860606512090082"
REF_LN2 = "0.693147180559945309417232121458"
REF_LNPI = "1.144729885849400174143427351354"
REF_ZETA3 = "1.202056903159594285399738161511"
REF_ZETA5 = "1.036927755143369926331365486457"


def test_gamma_against_reference():
    g = lam3.gamma_euler_iv(128)
    assert lam1.iv_overlap(g, _ref_bracket(REF_GAMMA))
    assert lam1.iv_width(g) < Fr(1, 10 ** 24)


def test_ln2_lnpi_against_reference():
    assert lam1.iv_overlap(lam3.ln2_iv(), _ref_bracket(REF_LN2))
    assert lam1.iv_overlap(lam3.ln_pi_iv(), _ref_bracket(REF_LNPI))


def test_zeta_values_against_reference():
    assert lam1.iv_overlap(lam3.zeta_iv(3, 64), _ref_bracket(REF_ZETA3))
    assert lam1.iv_overlap(lam3.zeta_iv(5, 64), _ref_bracket(REF_ZETA5))


def test_zeta_two_depths_overlap():
    assert lam1.iv_overlap(lam3.zeta_iv(3, 64), lam3.zeta_iv(3, 128))


def test_body_only_zeta_separates():
    body = sum(Fr(1, n ** 3) for n in range(1, 64))
    z = lam3.zeta_iv(3, 64)
    assert body < z[0] and z[0] - body > Fr(1, 10 ** 5)


def test_bernoulli_tamper_separates():
    tampered = dict(lam3.BERNOULLI)
    tampered[4] = Fr(1, 30)
    bad = lam3.gamma_euler_iv(64, bern=tampered)
    good = lam3.gamma_euler_iv(128)
    assert not lam1.iv_overlap(bad, good)


def test_E1_at_pi_scale_positive_and_sane():
    # E_1(3) external float reference (scipy): 1.3048381094197039e-02
    e = lam3.E1_at(Fr(3))
    ref = _ref_bracket("0.013048381094197039", slack=Fr(1, 10 ** 14))
    assert lam1.iv_overlap(e, ref)
    assert lam1.iv_width(e) < Fr(1, 10 ** 24)


def test_exact_cancellation_large_c():
    # I_{-2}(c) stays positive with tight width even at c ~ 16 pi,
    # where floats would suffer ~22-digit catastrophic cancellation
    pi_iv = lam1.pi_brackets()
    c = (pi_iv[0] * 16, pi_iv[1] * 16)
    im2 = lam3.I_m2_iv(c)
    assert im2[0] > 0
    assert lam1.iv_width(im2) < Fr(1, 10 ** 26)


def test_split_flip_identity_certified():
    cert = json.load(open(os.path.join(CERT_DIR, "LAM3_RESULT.json")))
    t2 = cert["T2_split_flip_identity_s3"]
    assert t2["verdict"] == "PASS" and t2["overlap"] is True
    # leading digits match the scratch-independent value 0.17561767...
    assert t2["lhs_enclosure"][0].startswith("0.17561767")


def test_completed_identities_certified():
    cert = json.load(open(os.path.join(CERT_DIR, "LAM3_RESULT.json")))
    assert cert["T3_completed_identity_s3"]["verdict"] == "PASS"
    assert cert["T4_completed_identity_s5"]["verdict"] == "PASS"
    assert cert["T3_completed_identity_s3"]["lhs_enclosure"][0] \
        .startswith("0.19131329")
    assert cert["T4_completed_identity_s5"]["lhs_enclosure"][0] \
        .startswith("0.07879706")


def test_controls_recorded_and_separated():
    cert = json.load(open(os.path.join(CERT_DIR, "LAM3_RESULT.json")))
    for key in ("controls_C1_polar_drop_separates",
                "controls_C2_seam_weight_tamper_separates",
                "controls_C5_exponent_tamper_separates",
                "controls_C3_body_only_zeta_separates",
                "controls_C4_bernoulli_tamper_separates"):
        assert cert[key]["separated"] is True


def test_certificate_pin_matches_regeneration():
    cert = lam3.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_LAM3.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_claim_boundary_open_and_honest():
    cert = json.load(open(os.path.join(CERT_DIR, "LAM3_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["N1_native_continuation"] == "OPEN"
    assert cb["N2_native_functional_equation"] == "OPEN"
    assert cb["N3_identification"] == "OPEN"
    assert cb["K0_L0_RH"] == "OPEN"
    assert cb["YM_continuum_gates"] == "OPEN"
    assert cb["zeros_of_any_L_function"] == "no claims"
    assert "NOT a native theorem" in cb["continuation_into_strip"]
