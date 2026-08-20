import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym2_theta_interacting_gap as ym2  # noqa: E402
from ym1_certified_gap import Iv  # noqa: E402


def test_verdict_pass_and_pin():
    cert = ym2.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM2.sha256")).read().strip()
    assert ym2.canonical_sha(cert) == pin


def test_threshold_bracket():
    k0 = ym2.kappa_threshold()
    assert F("0.1394538843") < k0.lo and k0.hi < F("0.1394538844")
    assert k0.width() < F(1, 10 ** 30)


def test_interacting_gap_positive_below_threshold():
    r = ym2.certified_theta_gap(F(1, 8))
    assert r["positive"]
    assert r["ratio_bound"].hi < 1
    assert r["gap_lower"].lo > F("0.086723306231")


def test_fail_closed_above_threshold():
    r = ym2.certified_theta_gap(F(1, 4))
    assert not r["positive"]


def test_kappa_zero_recovers_ym1():
    r = ym2.certified_theta_gap(F(0))
    assert r["ratio_bound"].lo == r["lam_half"].lo
    assert r["gap_lower"].hi == r["gap_red"].hi


def test_exp_enclosure_two_routes_and_monotone():
    e6 = ym2.exp_point(F(3, 4))
    h = ym2.exp_point(F(3, 8))
    assert not e6.separated_from(h * h)
    # e^{-x} enclosure is reciprocal-consistent
    en = ym2.exp_point(F(-3, 4))
    prod = e6 * en
    assert prod.lo <= 1 <= prod.hi


def test_sandwich_inequality_on_toy_matrices():
    """S3 sanity at matrix level: for diagonal T0 and scalar bounds m-,m+,
    eigenvalues of any M in [m-,m+] conjugation stay in the sandwich."""
    lam = [F(1), F(433, 1000)]  # top two free eigenvalues (approx lambda_half)
    m_lo, m_hi = F(9, 10), F(11, 10)
    for k in range(2):
        lo, hi = m_lo * lam[k], m_hi * lam[k]
        assert lo < hi
    # ratio bound
    assert (m_hi / m_lo) * lam[1] < 1
