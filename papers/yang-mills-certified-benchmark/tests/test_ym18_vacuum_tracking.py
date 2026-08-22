import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym18_vacuum_tracking as ym18  # noqa: E402


def test_verdict_and_pin():
    cert = ym18.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM18.sha256")).read().strip()
    assert ym18.canonical_sha(cert) == pin


def test_pair_element_first_order_is_f_half():
    B, P = ym18.compressed_symbolic(2)
    k = [i for i, b in enumerate(B) if b[0] == "P1"][0]
    assert P[0][k] == {tuple([0, 1] + [0] * (ym18.P_SYM - 1)): F(1)}


def test_engine_reduces_to_ym15_on_site_entries():
    import ym15_chain_closed_form as ym15
    B, P = ym18.compressed_symbolic(5)
    idx = {b[0]: i for i, b in enumerate(B)}
    assert P[idx["S1"]][idx["S4"]] == ym15.closed_form_entry(5, 1, 4)
    assert P[0][0] == ym15.closed_form_entry(5, 0, 0)


def test_enlarged_vacuum_exceeds_carrier_vacuum():
    M = ym18.compressed_T(3, F(1, 4))
    f0 = ym18.f0_of(F(1, 4))
    mu = ym18.iv_pow(f0, 2).hi * (1 + F(1, 10 ** 5))
    assert ym18.count_above(M, mu) >= 1
