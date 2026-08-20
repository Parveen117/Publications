import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ym4_symmetry_protected as ym4  # noqa: E402


def test_verdict_and_pin():
    cert = ym4.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YM4.sha256")).read().strip()
    assert ym4.canonical_sha(cert) == pin


def test_pairing_tensor_exact():
    assert ym4.pairing_g(1, 1, 1) == F(1, 2)
    assert ym4.pairing_g(2, 2, 2) == F(1, 3)
    assert ym4.pairing_g(1, 1, 0) == 0
    assert ym4.pairing_g(1, 0, 1) == 0
    # p<->q symmetry (T1 machine content)
    for p in range(4):
        for q in range(4):
            for r in range(4):
                assert ym4.pairing_g(p, q, r) == ym4.pairing_g(q, p, r)


def test_swap_commutes_on_all_grid():
    for kap in ym4.GRID:
        r = ym4.sector_spectra(kap)
        assert r["comm_zero"]


def test_sectors_split_with_even_rising():
    lo = ym4.sector_spectra(F(1, 8))
    hi = ym4.sector_spectra(F(1, 2))
    # even (crossing) branch rises, odd branch falls with coupling
    assert hi["even_second"].lo > lo["even_second"].hi - F(1, 100)
    assert hi["odd"].hi < lo["odd"].lo + F(1, 100)
    # strict separation of sectors at kappa=1/2
    assert hi["even_second"].lo > hi["odd"].hi


def test_character_ring():
    assert ym4.chi_mul(1, 1) == [0, 2]
    prod = ym4.ring_mul({1: F(1)}, {1: F(1)})
    assert prod == {0: F(1), 2: F(1)}


def test_truncation_error_small():
    assert ym4.tail_E(F(1, 2), ym4.P_CUT) < F(1, 10 ** 8)


def test_iv_sqrt_bracket():
    s = ym4.iv_sqrt(ym4.Iv(F(2)))
    a_lo = F("1.41421356237309504880168872420")   # < sqrt 2
    a_hi = F("1.41421356237309504880168872421")   # > sqrt 2
    assert s.lo < a_hi and s.hi > a_lo            # brackets overlap
    assert s.lo * s.lo <= 2 <= s.hi * s.hi        # true containment
