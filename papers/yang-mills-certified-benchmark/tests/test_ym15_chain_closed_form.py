import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym15_chain_closed_form as ym15  # noqa: E402


def test_verdict_and_pin():
    cert = ym15.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM15.sha256")).read().strip()
    assert ym15.canonical_sha(cert) == pin


def test_closed_form_symbolic_all_m():
    for m in range(2, 9):
        assert ym15.t1_symbolic_check(m)


def test_tamper_breaks_identity():
    assert not ym15.t1_symbolic_check(2, tamper=True)
    assert not ym15.t1_symbolic_check(5, tamper=True)


def test_vacuum_is_f0_power_and_mixed_zero():
    assert ym15.chain_entry(6, {}) == ym15.p_pow_f(0, 5)
    assert ym15.chain_entry(6, {3: 1}) == {}


def test_uniform_ceiling_and_refusal():
    _, _, floor, ceil = ym15.uniform_bounds(F(1, 2))
    assert ceil.hi < 1 and floor.lo > 0
    _, _, _, ceil2 = ym15.uniform_bounds(F(2))
    assert ceil2.hi >= 1


def test_lammax_monotone_and_below_ceiling():
    r = ym15.r_of(F(1))
    prev = None
    ceil = ((ym15.Iv(F(1)) + r) / (ym15.Iv(F(1)) - r)).hi
    for m in range(2, 7):
        b = ym15.lammax_bracket(m, r)
        assert b.hi < ceil
        if prev is not None:
            assert b.lo > prev.hi
        prev = b
