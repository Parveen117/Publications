import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym10_blindness_ledger as ym10  # noqa: E402
import ym9_uniform_heat_kernel as ym9  # noqa: E402


def test_verdict_and_pin():
    cert = ym10.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM10.sha256")).read().strip()
    assert ym10.canonical_sha(cert) == pin


def test_multiplicity_law_reproduces_carriers():
    assert ym10.carrier_dimension(ym10.V5) == 5
    assert ym10.carrier_dimension(ym10.V7) == 7
    assert ym10.multiplicity(1, 1) == 2      # (1/2,1/2) sector, YM-6 claim
    assert ym10.multiplicity(0, 0) == 1
    assert ym10.multiplicity(3, 1) == 2


def test_ym9_amendment_is_material_and_both_a_independent():
    assert ym10.uniform_count_contents(F(2)) == 4
    assert ym10.uniform_count_with_multiplicity(F(2)) == 5
    # YM-9's amended function reports both, and neither depends on a
    assert ym9.seam_count(F(2)) == 4
    assert ym9.seam_count(F(2), with_multiplicity=True) == 5


def test_v5_is_blindness_free_in_its_window():
    for s in [F(1), F(2)]:
        assert ym10.blind_dimension(ym10.V5, s) == 0
        assert ym10.probe_count(ym10.V5, s) == 0


def test_v7_strictly_less_blind_than_v5_where_it_matters():
    assert ym10.blind_dimension(ym10.V7, F(3)) < ym10.blind_dimension(ym10.V5, F(3))
    assert ym10.blind_dimension(ym10.V9, F(3)) < ym10.blind_dimension(ym10.V7, F(3))


def test_composition_increments_nonnegative():
    for s in ym10.S_GRID:
        assert (ym10.blind_dimension(ym10.V5, s)
                - ym10.blind_dimension(ym10.V7, s)) >= 0
        assert (ym10.blind_dimension(ym10.V7, s)
                - ym10.blind_dimension(ym10.V9, s)) >= 0


def test_multiplicity_tamper_changes_ledger():
    assert (ym10.blind_dimension(ym10.V5, F(4), tamper=True)
            != ym10.blind_dimension(ym10.V5, F(4)))


def test_honest_remainder_scopes_to_free_target():
    hr = ym10.run()["honest_remainder"]
    assert "LOWER bound" in hr["free_target_only"]
    assert hr["clay_predicate"] == "OPEN"
