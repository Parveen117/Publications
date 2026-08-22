import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym29_vacuum_floor as ym29  # noqa: E402


def test_verdict_and_pin():
    cert = ym29.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM29.sha256")).read().strip()
    assert ym29.canonical_sha(cert) == pin


def test_floor_above_one_and_zero_recovers_old():
    kap = F(1, 4)
    lam = ym29.lam_half()
    f0 = ym29.f0_of(kap)
    r = ym29._r(ym29.fj(1, kap) / f0)
    r1 = ym29._r(ym29.fj(2, kap) / f0)
    for m in (2, 7, 30):
        c = ym29.optimal_c(m, lam, r, r1)
        assert ym29.rayleigh_Q(m, c, lam, r, r1).lo > 1
        assert ym29.rayleigh_Q(m, F(0), lam, r, r1).lo == 1
