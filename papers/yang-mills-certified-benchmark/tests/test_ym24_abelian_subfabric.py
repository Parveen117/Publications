import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym24_abelian_subfabric as ym24  # noqa: E402


def test_verdict_and_pin():
    cert = ym24.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM24.sha256")).read().strip()
    assert ym24.canonical_sha(cert) == pin


def test_native_addition_law_on_circle():
    u = (F(3, 7), F(2, 7), F(6, 7))
    c1, s1 = ym24.circle_point(F(1, 2))
    c2, s2 = ym24.circle_point(F(-1, 3))
    prod = ym24.abelian_face(c1, s1, u) * ym24.abelian_face(c2, s2, u)
    assert prod == ym24.abelian_face(c1 * c2 - s1 * s2, s1 * c2 + c1 * s2, u)


def test_lagrange_and_second_order():
    p, q = (F(1), F(2), F(0)), (F(0), F(1), F(3))
    x = ym24.cross(p, q)
    assert ym24.dot(x, x) == ym24.dot(p, p) * ym24.dot(q, q) - ym24.dot(p, q) ** 2
    odd2, even2 = ym24.second_order([p, q])
    assert odd2 == ym24.vadd(ym24.vadd(p, q), x) and even2 == 1 - ym24.dot(p, q)


def test_non_abelian_energy_vanishes_on_parallel():
    base = (F(1), F(-2), F(3))
    assert ym24.e_na([ym24.vscale(base, F(k)) for k in range(1, 5)]) == 0
    assert ym24.e_na([(F(1), F(0), F(0)), (F(0), F(1), F(0))]) == 1
