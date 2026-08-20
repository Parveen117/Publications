import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym14_overlap_dock as ym14  # noqa: E402


def test_verdict_and_pin():
    cert = ym14.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM14.sha256")).read().strip()
    assert ym14.canonical_sha(cert) == pin


def test_bridged_pair_certified_counts():
    assert ym14.dock_cell(F(1, 8), F(3, 5))["count_exact"] == 1
    assert ym14.dock_cell(F(1, 2), F(9, 10))["count_exact"] == 1
    assert ym14.dock_cell(F(1), F(3, 2))["count_exact"] == 0


def test_first_order_anatomy_same_theta_integral():
    W, dp, dm = ym14.first_order_block()
    assert W[0][1] == F(1, 4)               # theta integral / 2
    assert W[2][2] == 0 and W[3][3] == 0    # B-plane inert
    assert dp.lo > 0 > dm.hi


def test_swap_and_free_limit():
    M, _ = ym14.m_bridge(F(1, 4))
    for i in range(5):
        for j in range(5):
            a, b = M[i][j], M[ym14.swap12(i)][ym14.swap12(j)]
            assert a.lo == b.lo and a.hi == b.hi
    M0, e0 = ym14.m_bridge(F(0))
    assert M0[0][0].lo == 1 and M0[1][2].hi == 0 and e0 == 0


def test_refusal_when_mu_below_complement():
    r = ym14.dock_cell(F(1, 2), F(1, 4))
    assert r["status"].startswith("REFUSED")
