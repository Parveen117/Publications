import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym17_interleaving_seam as ym17  # noqa: E402


def test_verdict_and_pin():
    cert = ym17.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM17.sha256")).read().strip()
    assert ym17.canonical_sha(cert) == pin


def test_pair_ratio_certified_and_edge_not():
    assert ym17.pair_ratio_bound(F(1, 4))[0] == ym17.NU
    assert ym17.pair_ratio_bound(F(2))[0] is None


def test_bracket_narrows():
    kap = F(1, 4)
    top = ym17.pair_top_upper(2 * kap)
    new_hi = (ym17.log_iv(ym17.Iv(top), ym17.LOG_TERMS) * ym17.Iv(F(1, 2))).hi
    f0 = ym17.f0_of(kap)
    assert ym17.log_iv(f0, ym17.LOG_TERMS).lo < new_hi < kap


def test_residual_price_above_one():
    for kap in (F(1, 8), F(1, 2)):
        p2 = ym17.iv_sqrt(ym17.Iv(ym17.pair_top_upper(2 * kap))).hi / ym17.f0_of(kap).lo
        assert p2 > 1
