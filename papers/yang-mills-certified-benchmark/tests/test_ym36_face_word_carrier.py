import os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym36_face_word_carrier as ym36  # noqa: E402


def test_small_matrix_exact_and_symmetric():
    words, T = ym36.word_matrix(F(1, 8), 3)
    assert T[(frozenset(), frozenset())] == 1
    assert all(T[(S, Sp)] == T[(Sp, S)] for S in words for Sp in words)


def test_verdict_and_pin():
    cert = ym36.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM36.sha256")).read().strip()
    assert ym36.canonical_sha(cert) == pin
