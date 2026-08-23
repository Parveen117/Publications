import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ng1_native_gap_statement as ng1  # noqa: E402


def test_verdict_and_pin():
    cert = ng1.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_NG1.sha256")).read().strip()
    assert ng1.canonical_sha(cert) == pin


def test_ng_is_not_claimed_to_be_clay():
    cert = ng1.run()
    assert "not_the_clay_statement" in cert["claim_status"].lower().replace("ng_is_", "")
    assert cert["obligations"]["DICT_translation_to_clay"].startswith("OPEN")
    assert cert["obligations"]["O4_time_direction_t"].startswith("OPEN")


def test_every_tooth_is_live():
    cert = ng1.run()
    assert len(cert["refutation_conditions"]) == 5
    for k, v in cert["refutation_conditions"].items():
        assert ("live" in v) or ("in force" in v)


def test_control_can_fail():
    # C1 must be capable of failing: an unknown tag makes it report missing
    saved = list(ng1.CITED)
    try:
        ng1.CITED.append("YM999")
        cert = ng1.run()
        assert cert["verdict"] == "FAIL"
        assert "YM999" in cert["controls"]["C1_missing"]
    finally:
        ng1.CITED[:] = saved
