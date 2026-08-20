import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym12_square_sourced as ym12  # noqa: E402


def test_verdict_and_pin():
    cert = ym12.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM12.sha256")).read().strip()
    assert ym12.canonical_sha(cert) == pin


def test_square_exponent_identities():
    assert ym12.exponent_square_identities(F(1, 4), F(1))
    assert ym12.exponent_square_identities(F(3, 7), F(2, 5))


def test_halving_tamper_separates():
    D = ym12.square_defect(F(1, 4))
    Dw = ym12.square_defect(F(1, 4), wrong_half=True)
    assert ym12.matrices_separated(D, Dw)
    assert ym12.defect_diag_nonneg(D)


def test_chain_ratio_invariance_exact():
    for m in (1, 2, 3, 5, 8):
        assert ym12.chain_ratio_invariant(m, F(1, 4))


def test_volume_uniform_bound_three_eighths():
    assert ym12.chain_uniform_bound(F(1, 16)) == F(3, 8)
    for a in (F(1), F(1, 100)):
        assert ym12.interval_route_chain(a, F(1, 16))


def test_planted_shared_face_contrast():
    from ym11_gate_verdicts import sandwich_volume_bound
    assert sandwich_volume_bound(4, F(1, 16)) == 0        # sharing kills sup
    assert ym12.chain_uniform_bound(F(1, 16)) == F(3, 8)  # chain survives


def test_governance_complete_and_corrections_recorded():
    g = ym12.GOVERNANCE
    assert set(g) == {f"YM{i}" for i in range(1, 13)}
    assert "standing_correction" in g["YM9"]
    assert "DEMOTED" in g["YM8"]["standing_correction"]
    assert g["YM12"]["consumers"] == ["TERMINAL"]
    assert len(ym12.DO_NOT_REOPEN) >= 3
