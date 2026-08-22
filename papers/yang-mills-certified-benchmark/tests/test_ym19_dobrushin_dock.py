import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym19_dobrushin_dock as ym19  # noqa: E402


def test_verdict_and_pin():
    cert = ym19.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM19.sha256")).read().strip()
    assert ym19.canonical_sha(cert) == pin


def test_certified_and_refused_cells():
    assert ym19.cell(F(6), F(1, 16))["status"] == "CERTIFIED_UNIFORM_IN_m"
    assert ym19.cell(F(6), F(1, 8))["status"].startswith("REFUSED")
    assert ym19.cell(F(10), F(1, 4))["status"].startswith("REFUSED")


def test_coefficient_bound_elementary():
    # delta = 1 -> coefficient 0 ; delta -> large -> coefficient -> 1
    assert ym19.coeff(ym19.Iv(F(1))).hi == 0
    assert ym19.coeff(ym19.Iv(F(100))).lo > F(99, 100)


def test_space_term_coupling_ceiling():
    # 2(1 - e^{-4k}) < 1 iff k < log2/4 ~ 0.1733
    assert (ym19.Iv(F(2)) * ym19.coeff(ym19.delta_s(F(17, 100)))).hi < 1
    assert (ym19.Iv(F(2)) * ym19.coeff(ym19.delta_s(F(18, 100)))).lo > 1


def test_exp_neg_large_argument():
    e = ym19.exp_neg(F(45))
    assert e.lo > 0 and e.hi < F(1, 10 ** 19)
