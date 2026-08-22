import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym28_e4d_on_chain as ym28  # noqa: E402


def test_verdict_and_pin():
    cert = ym28.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM28.sha256")).read().strip()
    assert ym28.canonical_sha(cert) == pin


def test_single_insertion_covariance_zero_formally():
    for m in (2, 5, 8):
        for i in range(1, m + 1):
            assert ym28.chain_entry(m, {i: 2}) == {}


def test_two_insertion_is_r1_power():
    m, i, j = 6, 2, 5
    pij = ym28.chain_entry(m, {i: 2, j: 2})
    mono = tuple([m - 1 - (j - i), 0, j - i] + [0] * (ym28.P - 2))
    assert pij == {mono: F(1)}
