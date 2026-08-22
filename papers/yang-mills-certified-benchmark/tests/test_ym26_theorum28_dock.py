import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym26_theorum28_dock as ym26  # noqa: E402


def test_verdict_and_pin():
    cert = ym26.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM26.sha256")).read().strip()
    assert ym26.canonical_sha(cert) == pin


def test_blind_error_crosses_one_at_mstar_plus_one():
    kap = F(1, 2)
    ms = ym26.predicted_m_star(kap)
    assert ym26.blind_error(kap, ms) < 1 <= ym26.blind_error(kap, ms + 1)


def test_ladder_tail_small_at_strong_coupling():
    assert ym26.ladder_tail_per_face(F(1, 8)).hi < F(1, 100)
