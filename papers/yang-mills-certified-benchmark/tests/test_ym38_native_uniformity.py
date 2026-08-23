import os, sys
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))
import ym38_native_uniformity as ym38  # noqa: E402
import ym37_space_transfer as ym37  # noqa: E402


def test_sqrt_bounds_certified():
    q = F(2)
    assert ym38.sqrt_up(q) ** 2 >= q >= ym38.sqrt_lo(q) ** 2


def test_native_count_is_t53_elimination():
    # count_pos_native must be the exact symmetric elimination, not Bareiss
    W = [[F(1, 2), F(1, 2)], [F(1, 2), F(1, 2)]]
    assert ym38.count_pos_native(W) == 1
    assert ym38.count_pos_native([[F(-1), F(0)], [F(0), F(-2)]]) == 0


def test_defect_B_orthogonality_exact():
    from ym31_uniform_floor import fj, lam, rnd_down
    Z, H, O = F(0), F(1, 2), F(1)
    kap = F(1, 8)
    fpt = {Z: rnd_down(fj(0, kap).lo), H: rnd_down(fj(1, kap).lo), O: rnd_down(fj(2, kap).lo)}
    klo = {Z: F(1), H: rnd_down(lam(1).lo), O: rnd_down(lam(2).lo)}
    basis, M, N, B, e0i = ym37.build_pencil(2, fpt, klo)
    n = len(basis)
    e0 = [F(0)] * n
    e0[e0i] = F(1)
    x = e0
    for _ in range(4):
        x = ym37.tau_apply(M, B, x)
    W_w = ym38.bq(B, x, x)
    th = ym38.bq(B, x, ym37.tau_apply(M, B, x)) / W_w
    tw = ym37.tau_apply(M, B, x)
    r = [tw[k] - th * x[k] for k in range(n)]
    assert ym38.bq(B, r, x) == 0


def test_verdict_and_pin():
    cert = ym38.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates", "EXPECTED_YM38.sha256")).read().strip()
    assert ym38.canonical_sha(cert) == pin
