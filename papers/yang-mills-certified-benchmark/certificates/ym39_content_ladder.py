"""YM-39: THE CONTENT LADDER ON THE SPACE TRANSFER — centre superselection
(exact) and the certified ordering of excitation channels, m-uniform.

Session item (b): contents >= 1 via the r_j ladder. YM-36's carrier held
only content-1/2 insertions; the honest worry was a lighter channel at
content >= 1. This capsule answers it at the carrier level with the
YM-38 native machinery — no new verdict types: T53 elimination sign
checks, the cut-square inequality, theorum/28's outward certificate.

 (T1) CENTRE SUPERSELECTION, EXACT. The centre Z2 of SU(2) acts on every
      rail as the total flip x -> -x (YM-37/38's parity). Faces and rungs
      are class functions of PRODUCTS of group elements — centre-blind —
      while chi_j picks the sign (-1)^{2j}. Hence every fabric quantity
      with one half-integer and one integer insertion vanishes EXACTLY:
      the mixed matrix N_mix (chi_half at one end, chi_1 at the other)
      is the ZERO matrix on the sector (checked entry by entry, exact
      rationals), and the YM-33 engine's mixed fabric partitions are 0
      for every tested (m, p). Half-integer and integer content are
      superselected on the fabric at every coupling, every m, every
      time order — a grading theorem, not an estimate. (This is the
      fabric instance of the odd/even split the operator ladder keeps
      meeting: the commuting zero-residue sector of theorum/51's
      exchange law is centre-even.)

 (T2) CHANNEL BULK RATIOS, NATIVE UNIFORMITY. For content c in
      {1/2, 1, 3/2}: N_c built from the insertion chi_c(x) chi_c(y)
      (chi_half = u, chi_1 = u^2 - 1, chi_3/2 = u^3 - 2u, u = 2 x0);
      rho_c(m,p) := (x_{m-p}' N_c x_{p-1}) / (x_0' B x_m). The YM-38
      derivation never used anything about N beyond its pencil norm and
      the positivity of the reference functional, so it applies verbatim
      per channel: |rho_c(m,p) - rho_c^ref| <= A_c (L^{j-p0} + L^{jc-p0})
      for ALL m, p with j = min(p-1, m-p) >= p0 = 2, constants certified
      per channel (same L — the transfer is the same; only A_c changes).

 (T3) LADDER VERDICT. Certified strict ordering, uniform in m:
          rho_{1/2} - err_{1/2}  >  rho_1 + err_1  >  rho_{3/2} + ...
      i.e. the content-1/2 excitation is the LIGHTEST (largest ratio,
      smallest gap) of the three channels on this carrier at every grid
      kappa, for every m and interior p. Combined with T1 (no mixing),
      the YM-36 concern is closed at the carrier level: no channel at
      content >= 1 undercuts the J = 1/2 sector, and none can mix into it.

Boundary (honest): carrier-level statements on the declared truncated
faces and rung kernel (contents <= 1 in both) — chi_3/2 is therefore a
COMPOSITE channel here (its rung propagation goes through fusion), its
row is data, not a fundamental-channel claim; compressed quantities
only; no operator upper bound; E4D-C stays OPEN. The heat-kernel
trajectory tie (YM-22: r_j/a -> 0 for j >= 1, so the ladder collapses
onto content 1/2 in the cutoff limit) is NAMED as the next quantified
step, not claimed here.

Controls:
  C1  N_mix identically zero (exact) AND engine mixed partitions zero
      for m = 2, 3, all p; a same-parity "tamper" (chi_half x chi_3/2,
      both half-integer) is NONZERO — the vanishing bites on parity,
      not on some accidental orthogonality.
  C2  per-channel positivity F_lo, Z_lo > 0 and ball/contract checks
      (theorum/28 outward certificate) as in YM-38.
  C3  per-channel exact tables: every (m,p) row inside its bound.
  C4  T3 ordering inequalities as exact rationals, with the uniformity
      margins included on both sides.
  C5  sanity anchor: rho_1 bulk within 25 percent of the declared
      lambda_1 = 9/64 scale at kappa = 1/8 (dressing-sized deviation),
      and rho_{1/2} matches YM-38's rho_c exactly (same construction).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(1000000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import _dec, canonical_sha  # noqa: E402
from ym31_uniform_floor import fj, lam, rnd_down  # noqa: E402
import ym37_space_transfer as y37  # noqa: E402
import ym38_native_uniformity as y38  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
ZERO, HALF, ONE = F(0), F(1, 2), F(1)
P0 = 2
K_EXACT = 34


def chi_poly(c, rail, t):
    """chi_c of rail's quaternion as polynomial: u = 2 x0."""
    u = y37.padd({}, y37.var(4 * rail, t), F(2))
    if c == HALF:
        return u
    u2 = y37.pmul(u, u)
    if c == ONE:
        return y37.padd(u2, y37.pconst(1, t), F(-1))
    if c == F(3, 2):
        u3 = y37.pmul(u2, u)
        return y37.padd(u3, u, F(-2))
    raise ValueError(c)


def build_N(basis, kcoef, ca, cb, t=2):
    K = y37.class_kernel(kcoef, y37.dot(0, 1, t), t)
    ins = y37.pmul(chi_poly(ca, 0, t), chi_poly(cb, t - 1, t))
    KI = y37.pmul(K, ins)
    n = len(basis)
    N = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            N[i][j] = N[j][i] = y37.integrate(y37.pmul(y37.pmul(basis[i], basis[j]), KI), t)
    return N


def channel_block(M, B, N, xs, w, W_w, theta, r_w, Rhat, sig2, tbar, L, Dsum, nsq0, b0):
    """the YM-38 assembly for one insertion matrix N (identical derivation)."""
    n = len(M)
    Wup = y38.sqrt_up(W_w)
    nu = y37.pencil_norm_hi(N, B, F(64), steps=40)
    Ntil = y38.dotv(w, y37.matvec(N, w))
    theta_lo = theta - Rhat * tbar / W_w
    K1 = (Rhat / (W_w * theta_lo)) * Dsum / (1 - L)
    KN = nu * 2 * (Wup + tbar) * Dsum / (1 - L)
    KZ = (y38.sqrt_up(nsq0) / b0 + tbar) * Dsum / (1 - L)
    F_lo = Ntil - nu * (2 * Wup * tbar + tbar * tbar) - KN
    Z_lo = W_w - (y38.sqrt_up(nsq0) / b0) * tbar - KZ
    ok = F_lo > 0 and Z_lo > 0 and K1 <= F(1, 4)
    if not ok:
        return None

    def rho(m, p):
        a = y38.dotv(xs[m - p], y37.matvec(N, xs[p - 1]))
        z = y38.dotv(xs[0], y37.matvec(B, xs[m]))
        return a / z
    rho_c = rho(26, 13)
    jc = 12
    rho_hi = (Ntil + nu * (2 * Wup * tbar + tbar * tbar) + KN) / (theta_lo * Z_lo)
    A = rho_hi * (4 * K1 + KN / F_lo + KZ / Z_lo) * 2
    table = {}
    ok_rows = True
    for (m, p) in ((26, 13), (30, 15), (34, 17), (30, 10), (16, 8), (13, 3)):
        j = min(p - 1, m - p)
        err = A * (L ** (j - P0) + L ** (jc - P0))
        val = rho(m, p)
        inside = abs(val - rho_c) <= err
        ok_rows = ok_rows and inside
        table[f"{m},{p}"] = {"rho": _dec(val, 16), "bound": _dec(err, 16), "inside": bool(inside)}
    # uniform two-sided band over ALL m, p with j >= P0
    band = A * (1 + L ** (jc - P0))
    return {"rho_ref_26_13": rho_c, "A": A, "band_all_interior": band,
            "table": table, "ok_rows": ok_rows}


def run():
    grid = {}
    C = [True] * 5
    for kap in GRID:
        fpt = {ZERO: rnd_down(fj(0, kap).lo), HALF: rnd_down(fj(1, kap).lo), ONE: rnd_down(fj(2, kap).lo)}
        klo = {ZERO: F(1), HALF: rnd_down(lam(1).lo), ONE: rnd_down(lam(2).lo)}
        basis, M, Nhalf_builtin, B, e0i = y37.build_pencil(2, fpt, klo)
        n = len(basis)
        e0 = [F(0)] * n
        e0[e0i] = F(1)
        xs = [e0]
        for _ in range(K_EXACT):
            xs.append(y37.tau_apply(M, B, xs[-1]))
        w = xs[12]
        W_w = y38.bq(B, w, w)
        tw = y37.tau_apply(M, B, w)
        theta = y38.bq(B, w, tw) / W_w
        r_w = [tw[k] - theta * w[k] for k in range(n)]
        assert y38.bq(B, r_w, w) == 0
        Rhat = y38.sqrt_up(y38.bq(B, r_w, r_w))
        U, Mp, Bp = y38.restricted(M, B, w, W_w)
        sig2 = y38.doubled_ceiling(Mp, Bp, hi=F(1))
        bs, nsqs, t_up = [], [], []
        for x in xs:
            b = y38.bq(B, w, x) / W_w
            yv = [x[k] - b * w[k] for k in range(n)]
            bs.append(b)
            nsqs.append(y38.bq(B, yv, yv))
            t_up.append(y38.sqrt_up(nsqs[-1]) / b)
        tbar = max(t_up[P0:]) * F(3, 2)
        theta_lo = theta - Rhat * tbar / W_w
        theta_hi = theta + Rhat * tbar / W_w
        inv_ok = theta_lo > 0 and (Rhat + sig2 * tbar) / theta_lo <= tbar
        L = sig2 / theta_lo + theta_hi * tbar * Rhat / (W_w * theta_lo * theta_lo)
        c2 = inv_ok and L < 1 and all(t <= tbar for t in t_up[P0:])
        yh = [[(xs[p][k] - bs[p] * w[k]) / bs[p] for k in range(n)] for p in range(K_EXACT + 1)]
        d0v = [yh[P0 + 1][k] - yh[P0][k] for k in range(n)]
        Dsum = y38.sqrt_up(y38.bq(B, d0v, d0v)) / (1 - L)
        # ---- T1: centre superselection ----
        Nmix = build_N(basis, klo, HALF, ONE)
        c1 = all(x == 0 for row in Nmix for x in row)
        Nsame = build_N(basis, klo, HALF, F(3, 2))       # both half-integer: must NOT vanish
        c1 = c1 and any(x != 0 for row in Nsame for x in row)
        from ym33_t_row_engine import fabric_partition
        for m in (2, 3):
            for p in range(1, m + 1):
                if fabric_partition([fpt, fpt], klo, m, 2, {(0, p): HALF, (1, p): ONE}) != 0:
                    c1 = False
        # ---- T2: channels ----
        chans = {}
        okpos = True
        for c in (HALF, ONE, F(3, 2)):
            Nc = build_N(basis, klo, c, c)
            if c == HALF and Nc != Nhalf_builtin:
                okpos = False
            blk = channel_block(M, B, Nc, xs, w, W_w, theta, r_w, Rhat, sig2, tbar, L, Dsum, nsqs[0], bs[0])
            if blk is None:
                okpos = False
                continue
            if not blk["ok_rows"]:
                C[2] = False
            chans[str(c)] = blk
        C[1] = C[1] and c2 and okpos
        C[0] = C[0] and c1
        # ---- T3: ordering with uniformity margins ----
        c4 = True
        order = [HALF, ONE, F(3, 2)]
        for a, b in zip(order, order[1:]):
            ra, rb = chans[str(a)], chans[str(b)]
            if not (ra["rho_ref_26_13"] - ra["band_all_interior"] > rb["rho_ref_26_13"] + rb["band_all_interior"]):
                c4 = False
        C[3] = C[3] and c4
        # ---- C5 sanity anchors ----
        c5 = True
        if kap == F(1, 8):
            r1 = chans[str(ONE)]["rho_ref_26_13"]
            c5 = abs(r1 - F(9, 64)) <= F(9, 64) / 4
        y38r = json.load(open(os.path.join(HERE, "YM38_RESULT.json")))
        rc38 = F(y38r["grid"][str(kap)]["j_ge_2"]["rho_c_exact_26_13"])
        slack = F(1, 10 ** 13)
        c5 = c5 and abs(chans[str(HALF)]["rho_ref_26_13"] - rc38) <= slack
        C[4] = C[4] and c5
        grid[str(kap)] = {
            "L_beta": _dec(L, 12),
            "channels": {c: {"rho_ref": _dec(v["rho_ref_26_13"], 16),
                             "uniform_band": _dec(v["band_all_interior"], 20),
                             "table": v["table"]} for c, v in chans.items()},
            "ordering_half_gt_one_gt_threehalf_with_margins": bool(c4),
        }
    ok = all(C)
    cert = {
        "certificate_type": "YM39_CONTENT_LADDER_SUPERSELECTION_AND_ORDERING",
        "claim_status": "centre_superselection_EXACT__channel_ratios_native_m_uniform__"
                        "J_half_lightest_certified_on_carrier__composite_3half_data_only__E4DC_OPEN",
        "grid": grid,
        "controls": {"C1_mixed_channel_vanishes_exact_and_parity_tamper_nonzero": bool(C[0]),
                     "C2_outward_certificate_and_channel_positivity": bool(C[1]),
                     "C3_all_table_rows_inside": bool(C[2]),
                     "C4_ordering_with_uniform_margins": bool(C[3]),
                     "C5_anchors_lambda1_scale_and_ym38_match": bool(C[4])},
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM39_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    with open(os.path.join(HERE, "EXPECTED_YM39.sha256"), "w") as f:
        f.write(canonical_sha(cert) + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(k, "L", v["L_beta"], {c: (d["rho_ref"], d["uniform_band"]) for c, d in v["channels"].items()},
              "order", v["ordering_half_gt_one_gt_threehalf_with_margins"])
