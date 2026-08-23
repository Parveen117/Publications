import os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym35_commutator_remainder as ym35  # noqa: E402


def test_adjacent_weights_and_commuting_model():
    w = ym35.site_content_weights([F(1, 2), F(1, 2)], [F(1, 2), F(1, 2)], site=1, contents=[F(0), F(1)])
    assert w[F(0)] == F(1, 4) and w[F(1)] == F(3, 4)
    V, Xs = ym35.t3_commuting_model(F(1, 8), F(1, 3), 5)
    assert V == (1 + F(1, 64) * F(1, 3)) ** 4 and len(set(X / V for X in Xs)) == 1


def test_verdict_and_pin():
    cert = ym35.run()
    assert cert["verdict"] == "PASS"
    assert cert["T1_commutator_locality_nonadjacent_no_mixing"]
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM35.sha256")).read().strip()
    assert ym35.canonical_sha(cert) == pin
