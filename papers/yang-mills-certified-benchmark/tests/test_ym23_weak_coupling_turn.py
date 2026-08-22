import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym23_weak_coupling_turn as ym23  # noqa: E402


def test_verdict_and_pin():
    cert = ym23.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM23.sha256")).read().strip()
    assert ym23.canonical_sha(cert) == pin


def test_tiling_route_sign_change():
    a = F(1)
    assert ym23.leading_rate_step(a, F(1)).lo > 0
    assert ym23.leading_rate_step(a, F(2)).hi < 0


def test_product_law_commutator_is_the_only_nonadditive_term():
    p1, p2 = (F(1), F(0), F(0)), (F(0), F(1), F(0))
    h = ym23.quat_from_odd(p1) * ym23.quat_from_odd(p2)
    assert (h.b, h.c, h.d) == (F(1), F(1), F(1))   # p1 + p2 + p1 x p2
    assert h.a == 1


def test_rayleigh_exact():
    m = 7
    v = [F(i) - F(m + 1, 2) for i in range(1, m + 1)]
    num = sum((v[i] - v[i + 1]) ** 2 for i in range(m - 1))
    den = sum(x * x for x in v)
    assert num / den == F(12, m * (m + 1))
