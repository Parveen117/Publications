import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym32_dressed_excitation as ym32  # noqa: E402


def test_verdict_and_pin():
    cert = ym32.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM32.sha256")).read().strip()
    assert ym32.canonical_sha(cert) == pin


def test_one_rail_insertion_vanishes_and_reversal():
    z, h, o = F(0), F(1, 2), F(1)
    u = {z: F(1), h: F(1, 3), o: F(1, 10)}
    k = {z: F(1), h: F(2, 5), o: F(1, 10)}
    assert ym32.ladder_partition_ins(u, u, k, 3, {2: h}, {}) == 0
    assert ym32.ladder_partition_ins(u, u, k, 3, {1: h}, {1: h}) == \
        ym32.ladder_partition_ins(u, u, k, 3, {3: h}, {3: h})
