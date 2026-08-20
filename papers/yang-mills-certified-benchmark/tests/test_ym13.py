import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym13_beta_audit as ym13  # noqa: E402
import ym13_independent_witnesses as W  # noqa: E402


def test_verdict_and_pin():
    cert = ym13.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM13.sha256")).read().strip()
    assert ym13.canonical_sha(cert) == pin


def test_beta_computed_from_gram():
    _, _, bb = ym13.gram_and_beta(ym13.H_BEFORE_OFF)
    _, _, ba = ym13.gram_and_beta(ym13.H_AFTER_OFF)
    assert bb == F(4, 5) and ba == F(3, 5)
    assert ym13.RH_BETA == F(1, 5) and ym13.RH_BETA < ba < bb


def test_cf_bessel_matches_series_independently():
    from ym1_certified_gap import bessel_I, TERMS
    lam = bessel_I(2, F(2), TERMS) / bessel_I(1, F(2), TERMS)
    lo, hi = W.bessel_ratio_cf(1, F(2))
    assert lo <= lam.hi and lam.lo <= hi
    assert hi - lo < F(1, 10 ** 50)


def test_witness_module_import_isolated():
    import ast
    src = open(os.path.join(HERE, "..", "certificates",
                            "ym13_independent_witnesses.py")).read()
    for node in ast.walk(ast.parse(src)):
        if isinstance(node, ast.ImportFrom):
            assert not (node.module or "").startswith("ym")
        if isinstance(node, ast.Import):
            assert all(not a.name.startswith("ym") for a in node.names)


def test_compound_exp_sound():
    import math
    for x in [F(3, 8), F(-1, 2)]:
        lo, hi = W.exp_bounds_compound(x)
        assert float(lo) <= math.exp(float(x)) <= float(hi)


def test_retained_overlaps_have_refusal_notes():
    cert = ym13.run()
    assert len(cert["retained_overlaps"]) == 3
    assert "H_23=1" in cert["retained_overlaps"]
