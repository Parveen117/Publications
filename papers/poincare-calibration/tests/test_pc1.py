import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import sympy as sp  # noqa: E402
import pc1_lambda_monotonicity as pc1  # noqa: E402


def test_verdict_and_pin():
    cert = pc1.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_PC1.sha256")).read().strip()
    assert pc1.canonical_sha(cert) == pin


def test_monotonicity_identity_is_exactly_zero():
    rem, _ = pc1.monotonicity_remainder()
    assert sp.simplify(rem) == 0


def test_tampers_break_the_identity():
    assert sp.simplify(pc1.monotonicity_remainder(tamper_R=True)[0]) != 0
    assert sp.simplify(pc1.monotonicity_remainder(tamper_flow=True)[0]) != 0


def test_round_sphere_exact_extinction_and_growth():
    a0 = F(1)
    a, lam = pc1.round_solution(a0, F(1, 4))
    assert a == F(1, 2) and lam == F(6)
    _, lam2 = pc1.round_solution(a0, F(2, 5))
    assert lam2 > lam
    assert pc1.round_solution(a0, F(1, 2))[0] == 0


def test_R_matches_closed_form_exactly():
    for (a, b, c) in [(F(1), F(1), F(1)), (F(3), F(2), F(1)),
                      (F(7, 3), F(1, 2), F(5))]:
        ref = (2 * (a * b + b * c + c * a)
               - (a * a + b * b + c * c)) / (a * b * c)
        assert pc1.scalar_R(a, b, c) == ref


def test_berger_sign_law_exact_both_sides():
    du_expr, (x, z) = pc1.berger_anisotropy_identity()
    assert sp.simplify(du_expr - (-4 * (x - z) / (x * z))) == 0
    assert (-4 * (F(2) - F(1)) / (F(2) * F(1))) < 0
    assert (-4 * (F(1) - F(2)) / (F(1) * F(2))) > 0


def test_ric_norm_strictly_positive():
    for s in [(F(1), F(1), F(1)), (F(3), F(2), F(1)), (F(1), F(1), F(1, 4))]:
        assert pc1.ric_norm_sq(*s) > 0


def test_honest_remainder_declared():
    cert = pc1.run()
    hr = cert["honest_remainder"]
    assert hr["poincare_conjecture"] == "NOT CLAIMED"
    assert any("surgery" in x for x in hr["obstruction_tower_untouched"])
    assert "never rederived" in cert["anchor_pinned"]
