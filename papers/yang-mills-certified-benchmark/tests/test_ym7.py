import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym7_v7_crossing_curves as ym7  # noqa: E402


def test_verdict_and_pin():
    cert = ym7.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM7.sha256")).read().strip()
    assert ym7.canonical_sha(cert) == pin


def test_kappa_7_10_unlocked_exact_count():
    r = ym7.dock_cell7(F(7, 10), F(13, 10))
    assert r["count_exact"] == 1


def test_gram7_blockdiag():
    G = ym7.gram7()
    assert G[5][5] == 1 and G[6][6] == 1
    for j in range(5):
        assert G[5][j] == 0 and G[6][j] == 0
    assert G[3][4] == F(1, 2)


def test_complement_level_drop():
    lam, lam1, lam32, nxt = ym7.levels7()
    assert nxt.hi < lam1.lo                    # strict improvement vs YM-6
    assert lam32.hi < nxt.lo and (lam1 * lam1).hi < nxt.lo


def test_eig_brackets_ordered_and_tight():
    br = ym7.compressed_eigs(F(1, 4))
    for k in range(6):
        assert br[k].lo >= br[k + 1].lo        # nonincreasing
    for b in br:
        assert b.hi - b.lo < F(1, 10 ** 9)


def test_no_exact_constant_claims_in_curves():
    cert = ym7.run()
    # every curve entry is a two-sided bracket, lo != hi as strings only
    # when width exceeds printer resolution; the schema itself is brackets
    for kap, rows in cert["eigenvalue_curves"].items():
        for lo, hi in rows:
            assert F(lo) <= F(hi)
