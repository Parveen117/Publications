import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym9_uniform_heat_kernel as ym9  # noqa: E402


def test_verdict_and_pin():
    cert = ym9.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM9.sha256")).read().strip()
    assert ym9.canonical_sha(cert) == pin


def test_free_gap_exactly_three_quarters_all_a():
    assert ym9.free_reduced_gap_exact() == F(3, 4)
    assert ym9.casimir(1) == F(3, 4) and ym9.casimir(2) == F(2)


def test_uniform_bound_exact_and_threshold():
    assert ym9.uniform_interacting_gap_bound(F(1, 16)) == F(3, 8)
    assert ym9.certifies_uniform(F(1, 16))
    assert not ym9.certifies_uniform(F(1, 8))       # exactly critical
    assert not ym9.certifies_uniform(F(1, 4))


def test_semigroup_and_subdivision_exact():
    assert ym9.semigroup_exponent_identity(F(1, 3), F(2, 7))
    for n in (2, 3, 10, 97):
        assert ym9.subdivision_exact(F(3, 5), n)


def test_interval_route_agrees_across_six_decades():
    for a in [F(1), F(1, 100), F(1, 10 ** 6)]:
        gap, agree = ym9.interval_route_check(a, F(1, 16))
        assert agree
        assert gap.lo <= F(3, 8) <= gap.hi


def test_seam_count_is_a_independent_combinatorial():
    # count depends only on s, never on a
    assert ym9.seam_count(F(1)) == 3
    assert ym9.seam_count(F(2)) == 4
    assert ym9.seam_count(F(4)) == 10
    assert ym9.seam_count(F(2), casimir_tamper=True) != ym9.seam_count(F(2))


def test_wilson_diverges_as_spacing_shrinks():
    w1 = ym9.wilson_reduced_gap(F(1))
    w4 = ym9.wilson_reduced_gap(F(1, 4))
    assert w4.lo > 3 * w1.hi


def test_honest_remainder_present():
    hr = ym9.run()["honest_remainder"]
    assert hr["clay_predicate"] == "OPEN"
    assert "DECLARED" in hr["trajectory"]
    assert "not the physical mass gap" in hr["object"]
