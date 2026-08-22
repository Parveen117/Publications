import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym22_cutoff_tiling_rate as ym22  # noqa: E402


def test_verdict_and_pin():
    cert = ym22.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM22.sha256")).read().strip()
    assert ym22.canonical_sha(cert) == pin


def test_limit_is_uniform_lower_bound():
    limit = F(23, 32)
    for a in (F(1), F(1, 1000), F(1, 10 ** 6)):
        g = ym22.gamma(a, F(1, 16))
        assert g.lo >= limit and g.hi <= F(3, 4)


def test_beats_ym9_factor_twelve():
    # YM-9: 3/4 - 6 theta ; tiling: 3/4 - theta/2
    theta = F(1, 16)
    assert ym22.gamma(F(1, 100), theta).lo > F(3, 4) - 6 * theta


def test_trajectory_ceiling():
    assert ym22.gamma(F(1, 10 ** 4), F(2)).hi < 0
    assert ym22.gamma(F(1, 10 ** 4), F(1)).lo > 0
