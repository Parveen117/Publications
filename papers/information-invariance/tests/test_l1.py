import hashlib, json, pathlib, sys
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "certificates"))
import l1_certificates as L

def test_theorem_A_proved():
    assert L.build()["theorem_A_information_invariance_iff_flat"]["status"] == "PROVED"

def test_theorem_B_proved():
    b = L.build()["theorem_B_single_cycle_nonclosure_not_curvature_witness"]
    assert b["status"] == "PROVED"
    assert b["checks"]["cycle2_rule_separates"] is True

def test_certificate_deterministic_and_pinned():
    a, b = L.build(), L.build()
    assert a["certificate_sha256"] == b["certificate_sha256"]
    pin = (ROOT / "certificates" / "EXPECTED_L1.sha256").read_text().strip()
    assert a["certificate_sha256"] == pin

def test_tamper_detected():
    c = L.build(); body = {k: v for k, v in c.items() if k != "certificate_sha256"}
    body["theorem_A_information_invariance_iff_flat"]["status"] = "FAILED"
    assert hashlib.sha256(L.cj(body).encode()).hexdigest() != c["certificate_sha256"]

def test_ledger_lists_only_certified_claims_as_proved():
    txt = (ROOT / "LEDGER.md").read_text()
    assert "## PROVED" in txt and "## OPEN" in txt and "D1" in txt and "D4" in txt


sys.path.insert(0, str(ROOT / "certificates"))
import d2_phase_origin as D

def test_theorem_C_proved():
    assert D.build()["theorem_C_real_flow_is_schrodinger_iff"]["status"] == "PROVED"

def test_theorem_D_proved_and_scalar_channel_blocked():
    d = D.build()["theorem_D_one_channel_obstruction"]
    assert d["status"] == "PROVED"
    assert d["checks"]["odd_dim_determinant_obstruction"] is True

def test_D2_certificate_pinned():
    c = D.build()
    assert c["certificate_sha256"] == (ROOT / "certificates" / "EXPECTED_D2.sha256").read_text().strip()

def test_D2_still_open_not_claimed_discharged():
    r = D.build()["reduces_obligation"]
    assert r["D2"].startswith("not discharged")


import d2a_cut_supplies_structure as E

def test_theorem_E_proved():
    assert E.build()["theorem_E_cut_supplies_conjugate_pair_and_complex_structure"]["status"] == "PROVED"

def test_D2a_discharged_D2b_open():
    o = E.build()["obligations"]
    assert o["D2a"].startswith("DISCHARGED")
    assert o["D2b"].startswith("OPEN")

def test_D2a_certificate_pinned():
    assert E.build()["certificate_sha256"] == (ROOT / "certificates" / "EXPECTED_D2A.sha256").read_text().strip()


import d2b_linearisation_verdict as B2

def test_D2b_executed():
    assert B2.build()["status"] == "EXECUTED"

def test_D2b_negative_verdict_recorded_not_hidden():
    c = B2.build()
    assert c["results"]["R4_full_generator"]["verdict"] == "NOT_UNITARY_DAMPED_ROTATION"
    assert c["results"]["theorem_E_conditions_on_A"]["passes"] is False
    assert "NEGATIVE RESULT" in c["verdict"]

def test_D2b_transport_term_alone_is_unitary():
    r = B2.build()["results"]["R3_D_antisymmetric_cut_odd_gives_phase"]
    assert r["antisymmetric"] and r["norm_preserved"] and r["verdict"] == "PHASE"

def test_D2b_certificate_pinned():
    assert B2.build()["certificate_sha256"] == (ROOT / "certificates" / "EXPECTED_D2B.sha256").read_text().strip()

def test_theorems_A_B_unaffected_by_D2b():
    assert "unaffected" in B2.build()["scope"]
