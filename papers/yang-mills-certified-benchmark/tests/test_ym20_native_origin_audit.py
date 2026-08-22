import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym20_native_origin_audit as ym20  # noqa: E402


def test_verdict_and_pin():
    cert = ym20.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM20.sha256")).read().strip()
    assert ym20.canonical_sha(cert) == pin


def test_ledger_closed_and_tamper_rejected():
    assert ym20.ledger_ok(ym20.LEDGER)
    bad = dict(ym20.LEDGER)
    del bad["YM-7"]
    assert not ym20.ledger_ok(bad)


def test_cayley_and_native_log():
    y = F(1, 3)
    assert ym20.cayley((1 + y) / (1 - y)) == y
    a = ym20.native_odd_log(ym20.Iv(y))
    cl = ym20.log_iv(ym20.Iv(F(2)), ym20.LOG_TERMS)   # (1+1/3)/(1-1/3)=2
    assert not (a.hi < cl.lo or cl.hi < a.lo)


def test_demoted_anchors_are_exactly_two():
    demoted = [k for k, v in ym20.LEDGER.items()
               if v[2] == "DEMOTED_CLASSICAL_IMPORT"]
    assert sorted(demoted) == ["YM-19", "YM-8"]
