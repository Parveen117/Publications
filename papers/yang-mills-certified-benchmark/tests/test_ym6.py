import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym6_seam_integer_dock as ym6  # noqa: E402


def test_verdict_and_pin():
    cert = ym6.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM6.sha256")).read().strip()
    assert ym6.canonical_sha(cert) == pin


def test_gram_exact_structure():
    G = ym6.gram5()
    for i in range(3):
        for j in range(5):
            assert G[i][j] == (F(1) if i == j else F(0))
    assert G[3][3] == 1 and G[4][4] == 1
    assert G[3][4] == F(1, 2) and G[4][3] == F(1, 2)


def test_exact_count_beyond_ym5_reach():
    r = ym6.dock_cell(F(1, 2), F(1))
    assert r["status"] == "OK" and r["count_exact"] == 1


def test_low_mu_counts_three_or_more():
    r = ym6.dock_cell(F(1, 8), F(3, 10))
    assert r["count_lower"] >= 3


def test_refusal_recorded_fail_closed():
    r = ym6.dock_cell(F(7, 10), F(5, 4))
    assert r["status"] == "OK" and not r["certified"]
    r2 = ym6.dock_cell(F(1, 2), F(1, 2))     # mu below |C| bound
    assert r2["status"].startswith("REFUSED")


def test_inertia_engine_known_cases():
    def ivm(rows):
        return [[ym6.Iv(F(x)) for x in row] for row in rows]
    assert ym6.ldl_inertia(ivm([[5]])) == (1, 0)
    assert ym6.ldl_inertia(ivm([[1, 2], [2, 1]])) == (1, 1)
    assert ym6.ldl_inertia(ivm([[2, 1, 0], [1, 2, 1], [0, 1, 2]])) == (3, 0)
    assert ym6.ldl_inertia([[ym6.Iv(F(-1), F(1))]]) is None


def test_general_pairing_formula():
    # Int chi_a chi_b chi_c(AB^-1) = delta/d via the slot machinery at k=0
    G = ym6.gram5()
    # <phi_4, phi_5> = Int chi12(A) chi12(B) chi12(AB^-1) = 1/2 (YM-3 L4)
    assert G[3][4] == F(1, 2)
