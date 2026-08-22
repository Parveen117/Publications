import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym21_tiling_law as ym21  # noqa: E402


def test_verdict_and_pin():
    cert = ym21.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM21.sha256")).read().strip()
    assert ym21.canonical_sha(cert) == pin


def test_negative_witness_boundary_zero_faces_positive():
    g = ym21.rational_unit(F(1), F(2), F(-1, 3))
    assert (g * g.inv()).residue() == 0
    assert g.residue() + g.inv().residue() > 0


def test_subtelescoping_exact():
    rungs = [ym21.rational_unit(F(k), F(1, k + 1), F(2, k + 2)) for k in range(1, 7)]
    faces = ym21.ladder_faces(rungs)
    prod = ym21.QONE
    for i in range(1, 4):      # faces l_2..l_4 -> A_2^{-1} A_5
        prod = prod * faces[i]
    assert prod == rungs[1].inv() * rungs[4]


def test_leading_order_rate_matches_ym15_bracket():
    r = ym21.r_of(F(1, 4))
    lam = ym21.lam_half()
    one = ym21.Iv(F(1))
    up = lam * (one + r) / (one - r)
    assert up.hi < 1
