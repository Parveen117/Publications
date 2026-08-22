import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym16_chain_dock as ym16  # noqa: E402


def test_verdict_and_pin():
    cert = ym16.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM16.sha256")).read().strip()
    assert ym16.canonical_sha(cert) == pin


def test_finite_volume_certified_and_refused():
    assert ym16.dock_cell(F(1, 2), 4)["count_exact"] == 1
    assert ym16.dock_cell(F(1, 2), 5)["status"].startswith("REFUSED")
    assert ym16.dock_cell(F(1, 8), 13)["count_exact"] == 1


def test_mstar_formula_matches_dock():
    for kap in (F(1, 8), F(1, 4), F(1, 2)):
        ms = ym16.predicted_m_star(kap)
        assert ym16.dock_cell(kap, ms)["count_exact"] == 1
        assert ym16.dock_cell(kap, ms + 1)["status"].startswith("REFUSED")


def test_free_chain_all_m():
    for m in (2, 3, 7):
        assert ym16.dock_cell(F(0), m)["count_exact"] == 1


def test_sup_price_exceeds_one():
    for kap in (F(1, 8), F(1)):
        lo = ym16.exp_point(kap).lo
        assert lo / ym16.f0_of(kap).hi > 1
