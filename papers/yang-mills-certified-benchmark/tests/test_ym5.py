import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym5_two_sided_gap as ym5  # noqa: E402


def test_verdict_and_pin():
    cert = ym5.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM5.sha256")).read().strip()
    assert ym5.canonical_sha(cert) == pin


def test_two_sided_gap_beyond_sandwich_threshold():
    # kappa = 1/4 > kappa_0 = 0.1394...: YM-2 could not certify; YM-5 does
    r = ym5.certified_two_sided(F(1, 4))
    assert r["certifies"]
    assert r["ratio_hi"] < F(3, 5)
    assert r["gap_lower"] > F(1, 2)


def test_fail_closed_at_large_kappa():
    r = ym5.certified_two_sided(F(1))
    assert not r["certifies"]


def test_offblock_gram_upper_bound_positive_and_small():
    lam, _, _ = ym5.lam_levels()
    psq, gmin, gmax = ym5.offblock_bound(F(1, 8), (lam * lam).hi)
    assert gmax > 0
    assert psq < F(1, 10)          # off-block coupling genuinely small
    assert gmin > F(-1, 25)


def test_levels_certified():
    lam, lam1, lam32 = ym5.lam_levels()
    nxt = lam * lam
    assert lam1.hi < nxt.lo and lam32.hi < nxt.lo and nxt.hi < lam.lo


def test_outward_rounding_sound():
    x = ym5.Iv(F(1, 3), F(1, 3))
    r = ym5._r(x)
    assert r.lo <= F(1, 3) <= r.hi and r.hi - r.lo <= F(2, 10 ** 30)
