import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym31_uniform_floor as ym31  # noqa: E402


def test_verdict_and_pin():
    cert = ym31.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM31.sha256")).read().strip()
    assert ym31.canonical_sha(cert) == pin


def test_rung_truncation_is_monotone_lower_bound():
    zero, half, one = F(0), F(1, 2), F(1)
    u = {zero: F(1), half: F(1, 30), one: F(1, 1500)}
    k = {zero: F(1), half: F(43, 100), one: F(1, 10)}
    k2 = dict(k)
    k2[one] = F(0)
    assert ym31.ladder_partition(u, u, k2, 4) <= ym31.ladder_partition(u, u, k, 4)


def test_floor_rate_above_old_floor_and_e4d():
    cert = ym31.run()
    for kap, row in cert["grid"].items():
        assert F(row["per_face_floor_rate_all_m_ge_m0"]) > F(row["old_floor_rate_log_f0"])
        assert F(row["margin_per_face"]) > 0
