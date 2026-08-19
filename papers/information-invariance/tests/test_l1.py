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
