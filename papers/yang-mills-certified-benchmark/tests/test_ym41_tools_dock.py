import os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym41_tools_dock as ym41  # noqa: E402


def test_characters_native_recurrence():
    assert ym41.chi_of_u(F(3, 2)) == [F(0), F(-2), F(0), F(1)]
    assert ym41.chi_of_u(F(2)) == [F(1), F(0), F(-3), F(0), F(1)]


def test_char_component_reassembly_spot():
    m = (2, 0, 0, 0)
    total = {}
    for c in (F(0), F(1, 2), F(1)):
        for e, v in ym41.char_component(m, c).items():
            total[e] = total.get(e, F(0)) + v
    # compare as functions against probes
    import ym37_space_transfer as y37
    for p in ((0, 0, 0, 0), (2, 0, 0, 0), (0, 2, 0, 0)):
        lhs = y37.moment4(tuple(a + b for a, b in zip(m, p)))
        rhs = sum(v * y37.moment4(tuple(a + b for a, b in zip(e, p))) for e, v in total.items())
        assert lhs == rhs


def test_verdict_and_pin():
    cert = ym41.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM41.sha256")).read().strip()
    assert ym41.canonical_sha(cert) == pin
    for b in cert["meter_column_pencil"]["flips"]:
        assert (abs(b["jump"]) % 2 == 1) == b["det_flip"]
