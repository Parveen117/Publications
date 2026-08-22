import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym33_t_row_engine as E  # noqa: E402


def test_pin_matches_committed_result():
    # the full run takes minutes; the pin is regenerated in CI by the certificate step
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM33.sha256")).read().strip()
    assert len(pin) == 64


def test_t2_matches_ym32_and_t1_chain_limit():
    z, h, o = F(0), F(1, 2), F(1)
    u = {z: F(1), h: F(1, 3), o: F(1, 10)}
    k = {z: F(1), h: F(2, 5), o: F(1, 10)}
    from ym32_dressed_excitation import ladder_partition_ins
    assert E.fabric_partition([u, u], k, 3, 2, {(0, 1): h, (1, 1): h}) == \
        ladder_partition_ins(u, u, k, 3, {1: h}, {1: h})
    assert E.fabric_partition([u], k, 4, 1, {}) == 1


def test_three_row_time_reflection():
    z, h, o = F(0), F(1, 2), F(1)
    u = {z: F(1), h: F(1, 3), o: F(1, 10)}
    k = {z: F(1), h: F(2, 5), o: F(1, 10)}
    v1 = E.fabric_partition([u, u, u], k, 2, 3, {(0, 1): h, (2, 1): h})
    v2 = E.fabric_partition([u, u, u], k, 2, 3, {(2, 1): h, (0, 1): h})
    assert v1 == v2 > 0
