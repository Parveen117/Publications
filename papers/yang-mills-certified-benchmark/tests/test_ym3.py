import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym3_crossing_direction as ym3  # noqa: E402


def test_verdict_and_pin():
    cert = ym3.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM3.sha256")).read().strip()
    assert ym3.canonical_sha(cert) == pin


def test_theta_integral_exact_half():
    assert ym3.theta_integral() == F(1, 2)


def test_crossing_direction_rank_one_symmetric():
    cert = ym3.run()
    assert cert["exact_values"]["crossing_vector"] == [1, 1]
    assert cert["exact_values"]["receding_vector"] == [1, -1]
    assert cert["exact_values"]["block_eigen_derivatives"] == ["1/4", "-1/4"]


def test_slope_enclosure_and_slack():
    cert = ym3.run()
    lo = F(cert["enclosures"]["r_slope_lo"])
    hi = F(cert["enclosures"]["r_slope_hi"])
    # r'(0) = lambda_half / 4 ~ 0.108281856...
    assert F("0.1082") < lo <= hi < F("0.1083")
    assert cert["exact_values"]["slack_factor_vs_sandwich"] == "24"


def test_multiplicity_two_separation():
    cert = ym3.run()
    assert cert["enclosures"]["separation_next_levels"] is True


def test_run_all_pins_consistent():
    out = ym3.run_all_and_pin()
    assert set(out) == {"YM1", "YM2", "YM3"}
    assert all(v == "PASS" for v, _ in out.values())
