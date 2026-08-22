import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym25_direction_change as ym25  # noqa: E402


def test_verdict_and_pin():
    cert = ym25.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM25.sha256")).read().strip()
    assert ym25.canonical_sha(cert) == pin


def test_recoupling_squares_exact():
    R = ym25.recoupling_squares()
    assert R == [[F(1, 4), F(3, 4)], [F(3, 4), F(1, 4)]]


def test_scheme_vectors_are_J_half_and_orthogonal():
    assert ym25.vdot(ym25.scheme_12(0), ym25.scheme_12(1)) == 0
    assert ym25.vdot(ym25.scheme_23(0), ym25.scheme_23(1)) == 0
    # singlet(12) x up(3): norm 2 ; triplet-coupled: norm 6
    assert ym25.vdot(ym25.scheme_12(0), ym25.scheme_12(0)) == 2
    assert ym25.vdot(ym25.scheme_12(1), ym25.scheme_12(1)) == 6


def test_tamper_breaks_unitarity():
    Rt = ym25.recoupling_squares(tamper=True)
    assert not all(sum(row) == 1 for row in Rt)
