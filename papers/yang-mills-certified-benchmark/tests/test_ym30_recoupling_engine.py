import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym30_recoupling_engine as E  # noqa: E402


def test_verdict_and_pin():
    cert = E.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM30.sha256")).read().strip()
    assert E.canonical_sha(cert) == pin


def test_clebsch_gordan_exact():
    h = F(1, 2)
    assert E.cg(h, h, h, -h, F(0), F(0)).t == {2: F(1, 2)}          # 1/sqrt2
    assert E.cg(h, h, h, h, F(1), F(1)).t == {1: F(1)}
    assert E.cg(F(1), F(1), h, -h, h, h).t == {6: F(1, 3)}           # sqrt(2/3)


def test_one_plaquette_and_recoupling_split():
    h, z, o = F(1, 2), F(0), F(1)
    assert E.ladder_eval([h], [h], [h, h]).rational() == F(1, 4)
    assert E.ladder_eval([o], [o], [o, o]).rational() == F(1, 9)
    n0 = E.ladder_eval([h, h], [h, h], [h, z, h]).rational()
    n1 = E.ladder_eval([h, h], [h, h], [h, o, h]).rational()
    assert (n0, 3 * n1) == (F(1, 16), F(3, 16))


def test_surd_cancellation_is_enforced():
    s = E.Surd.sqrt_of(F(2)) * E.Surd.sqrt_of(F(2))
    assert s.is_rational() and s.rational() == 2
