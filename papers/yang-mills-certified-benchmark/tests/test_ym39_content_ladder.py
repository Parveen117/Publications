import os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym39_content_ladder as ym39  # noqa: E402
import ym37_space_transfer as ym37  # noqa: E402


def test_characters_exact():
    # chi_1 = chi_half^2 - 1 and chi_3/2 = chi_half^3 - 2 chi_half as polynomials
    t = 2
    u = ym39.chi_poly(F(1, 2), 0, t)
    c1 = ym39.chi_poly(F(1), 0, t)
    c32 = ym39.chi_poly(F(3, 2), 0, t)
    assert c1 == ym37.padd(ym37.pmul(u, u), ym37.pconst(1, t), F(-1))
    assert c32 == ym37.padd(ym37.pmul(ym37.pmul(u, u), u), u, F(-2))


def test_mixed_matrix_vanishes_same_parity_does_not():
    from ym31_uniform_floor import fj, lam, rnd_down
    Z, H, O = F(0), F(1, 2), F(1)
    kap = F(1, 8)
    fpt = {Z: rnd_down(fj(0, kap).lo), H: rnd_down(fj(1, kap).lo), O: rnd_down(fj(2, kap).lo)}
    klo = {Z: F(1), H: rnd_down(lam(1).lo), O: rnd_down(lam(2).lo)}
    basis, M, N, B, e0 = ym37.build_pencil(2, fpt, klo)
    Nmix = ym39.build_N(basis, klo, H, O)
    assert all(x == 0 for row in Nmix for x in row)
    Nsame = ym39.build_N(basis, klo, H, F(3, 2))
    assert any(x != 0 for row in Nsame for x in row)


def test_verdict_and_pin():
    cert = ym39.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM39.sha256")).read().strip()
    assert ym39.canonical_sha(cert) == pin
