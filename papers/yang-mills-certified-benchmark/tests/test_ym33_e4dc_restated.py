import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym33_e4dc_restated as ym33  # noqa: E402


def test_verdict_and_pin():
    cert = ym33.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM33.sha256")).read().strip()
    assert ym33.canonical_sha(cert) == pin
