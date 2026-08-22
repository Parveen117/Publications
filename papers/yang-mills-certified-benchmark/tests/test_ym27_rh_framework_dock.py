import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym27_rh_framework_dock as ym27  # noqa: E402


def test_verdict_and_pin():
    cert = ym27.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM27.sha256")).read().strip()
    assert ym27.canonical_sha(cert) == pin


def test_e5a_wrap_law_and_winding():
    N = ym27.chain_product([F(1, 2), F(3, 4), F(3, 4)])   # sum = 2 -> two turns
    assert N == (F(0), F(0), F(2))
    A = (F(1, 3), F(0), F(0))
    assert ym27.ugd_mul(A, ym27.ugd_inv(A)) == (F(0), F(0), F(0))


def test_aliasing_threshold():
    assert ym27.aliases(6, 8) == []          # j=3, n=8 > 6: exact
    assert ym27.aliases(8, 8) == [-8, 8]     # n <= 2j: aliases
