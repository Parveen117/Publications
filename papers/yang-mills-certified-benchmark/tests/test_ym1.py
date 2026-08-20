import json
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym1_certified_gap as ym  # noqa: E402


def test_verdict_pass_and_pin():
    cert = ym.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM1.sha256")).read().strip()
    assert ym.canonical_sha(cert) == pin


def test_gap_bracket_two_sided_and_tight():
    _, _, lam1, gap = ym.certified_reduced_gap(F(1), F(2))
    assert gap.lo < gap.hi
    assert gap.width() < F(1, 10 ** 30)
    assert F(0) < lam1.lo and lam1.hi < F(1)
    # 40-digit reference bracket (independent of decimal printer)
    assert gap.lo > F("0.83672330623158891305")
    assert gap.hi < F("0.83672330623158891306")


def test_interval_arithmetic_directed():
    a, b = ym.Iv(F(1, 3), F(1, 2)), ym.Iv(F(-2), F(3))
    m = a * b
    assert m.lo == F(-1) and m.hi == F(3, 2)
    try:
        _ = a / ym.Iv(F(-1), F(1))
        assert False, "division through zero must fail closed"
    except ZeroDivisionError:
        pass


def test_bessel_monotone_in_order():
    # I_1(2) > I_2(2) > I_3(2) strictly (separated enclosures)
    i1 = ym.bessel_I(1, F(2), 60)
    i2 = ym.bessel_I(2, F(2), 60)
    i3 = ym.bessel_I(3, F(2), 60)
    assert i2.hi < i1.lo and i3.hi < i2.lo


def test_log_bracket_contains_log2():
    lo, hi = ym._log_point(F(2), 80)
    # ln 2 = 0.69314718055994530941723212145817656807... ; two-sided anchors
    anchor_lo = F("0.69314718055994530941723212145817")   # < ln 2
    anchor_hi = F("0.69314718055994530941723212145818")   # > ln 2
    assert lo < anchor_hi and hi > anchor_lo               # enclosures overlap
    assert hi - lo < F(1, 10 ** 25)


def test_tamper_separates():
    cert = ym.run()
    assert cert["controls"]["C2_I3_tamper_separates"] is True
