import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym11_gate_verdicts as ym11  # noqa: E402


def test_verdict_and_pin():
    cert = ym11.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM11.sha256")).read().strip()
    assert ym11.canonical_sha(cert) == pin


def test_gauge_counting_reproduces_theta_graph():
    assert ym11.betti_one(3) == 2 and ym11.n_faces(3) == 3
    for n in range(2, 10):
        assert ym11.betti_one(n) == n - 1
        assert ym11.n_faces(n) == n * (n - 1) // 2


def test_free_gap_volume_and_cutoff_uniform():
    assert ym11.free_gap_any_volume() == F(3, 4)


def test_sandwich_critical_volume_exact():
    assert ym11.sandwich_volume_bound(3, F(1, 16)) == F(3, 8)
    assert ym11.sandwich_volume_bound(4, F(1, 16)) == F(0)
    assert ym11.sandwich_survives(3, F(1, 16))
    assert not ym11.sandwich_survives(4, F(1, 16))
    assert ym11.critical_volume(F(1, 16)) == 4


def test_smaller_theta_pushes_critical_volume_out():
    n16 = ym11.critical_volume(F(1, 16))
    n64 = ym11.critical_volume(F(1, 64))
    assert n64 > n16          # weaker coupling survives more volume


def test_counting_layer_regulator_independent():
    a = ym11.ledger_fingerprint()
    b = ym11.ledger_fingerprint()
    assert a == b
    # heat and Wilson index the same content set
    assert set(ym11.heat_coeffs(F(1))) == set(ym11.wilson_coeffs(F(2)))


def test_metric_layer_actually_differs():
    hg = ym11.heat_reduced_gap(F(1))
    wg = ym11.wilson_reduced_gap(F(1))
    assert wg.lo > hg or wg.hi < hg


def test_positive_coefficient_class_enforced():
    assert ym11.coeffs_positive(ym11.wilson_coeffs(F(2)))


def test_scope_warning_and_open_list_present():
    cert = ym11.run()
    assert "not claimed" in cert["scope_warning"].lower()
    assert "Clay predicate" in cert["still_open"]
    assert "interacting half of gate 2" in cert["still_open"]
