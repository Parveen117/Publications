import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym8_all_coupling_capstone as ym8  # noqa: E402


def test_verdict_and_pin():
    cert = ym8.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM8.sha256")).read().strip()
    assert ym8.canonical_sha(cert) == pin


def test_kernel_floor_positive_at_strong_coupling():
    fl = ym8.kernel_floor(F(2))
    assert fl.lo > 0
    # floor decreases with coupling
    assert ym8.kernel_floor(F(0)).lo > ym8.kernel_floor(F(1)).hi


def test_separations_strictly_positive():
    cert = ym8.run()
    for kap, sep in cert["lambda1_lambda2_separations"].items():
        assert F(sep) > 0


def test_honest_remainder_names_gates():
    cert = ym8.run()
    assert cert["honest_remainder"]["clay_predicate"] == "OPEN"
    assert any("UV" in g for g in cert["honest_remainder"]["open_gates"])
    assert "not uniform" in cert["honest_remainder"]["statement"].lower()


def test_anchor_pinned_not_rederived():
    cert = ym8.run()
    assert "Jentzsch" in cert["anchor_pinned"]
    assert "not rederived" in cert["anchor_pinned"]
