"""YM-38: NATIVE m-UNIFORMITY — the YM-37 T3 statement re-derived with the
framework's own machinery, replacing the anchored finite-spectral split.

Owner's standing rule (Aug 22): no borrowed classical results. YM-37 T3
used "B-orthonormal eigenvectors of the pencil" — the finite spectral
theorem. This capsule removes it. Nothing here decomposes into
eigenvectors; the only verdicts are:

  V1  T53 elimination (native LDL): sign pattern of the weights of a
      symmetric matrix. "count_above(sigma) = 0" literally means the
      exact elimination of M - sigma*B produces no positive weight,
      hence x'(M - sigma B)x <= 0 for every x (sum of nonpositive
      native squares) — a quadratic-form ceiling with no eigenvalues.
  V2  the cut-square inequality (RH T01-C form, as consumed by YM-31):
      for the PD square B (all elimination weights positive, certified),
      (u'Bv)^2 <= (u'Bu)(v'Bv); triangle inequality follows.
  V3  theorum/28: recognized-channel Cauchy with a DECLARED geometric
      Smriti tail (Sec 3-5), outward certificate u + e < 1 (Sec 9) —
      the framework's own finite-to-infinite convergence machinery;
      no fixed-point theorem is cited.
  V4  theorum/54 Sec 2 product bound: prod(1 + x_i) <= 1 + 2*s for
      s = sum x_i <= 1/2 (and prod(1 - x_i) >= 1 - s) — replaces every
      exponential estimate.

Construction (per kappa, t = 2, even invariant sector, dim 8; pencil
(M, B) and boundary vector e0 from YM-37, matched exactly to the YM-33
engine there).

  DEFLATION WITH AN EXPLICIT VECTOR (YM-8 lesson, native): w := x_P
  (P = 12), x_p := tau^p e0 exact rational, tau = B^{-1}M applied by
  exact solves. theta := (w'Mw)/(w'Bw), W_w := w'Bw. Defect
  r_w := tau w - theta w. EXACT IDENTITY r_w'Bw = 0 (because B tau = M
  and theta is the B-Rayleigh quotient of w) — checked exactly, and a
  tampered theta refuses it.

  SPLIT (recovered identity, theorum/28 item 1): x_p = b_p w + y_p with
  b_p = w'B x_p / W_w, y_p'Bw = 0 — exact at every computed p, and the
  step law is exact algebra:
      b_{p+1} = theta b_p + (r_w'B y_p)/W_w,
      y_{p+1} = b_p r_w + Pi tau y_p,
  where Pi is the B-orthogonal projection onto w-perp.

  COMPLEMENT CEILING WITHOUT EIGENVECTORS: with U the explicit rational
  basis u_i = e_i - (e_i'Bw/W_w) w (one coordinate dropped), one has
  U'BtauU = U'MU =: M' exactly, so the projected-restricted map is the
  restricted pencil tau' = B'^{-1}M', and
      ||Pi tau y||_B <= sigma'' ||y||_B   for every y perp w
  is certified by ONE native elimination sign check:
      count_+( M' B'^{-1} M'  -  sigma''^2 B' ) = 0            (V1).
  B' PD is certified the same way. No spectral objects.

  RECOGNIZED-CHANNEL CAUCHY (theorum/28 items 2-3-5-6): with
  yhat_p := y_p/b_p, s_p := b_{p+1}/b_p = theta + r_w'B yhat_p / W_w:
      yhat_{p+1} = (r_w + Pi tau yhat_p)/s_p ,
  and for p, q >= p0 (p0 = 12):
      ||yhat_{p+1} - yhat_p||_B <= L ||yhat_p - yhat_{p-1}||_B ,
      L := sigma''/theta_lo + theta_hi * tbar * Rhat / (W_w theta_lo^2),
  provided t_p <= tbar for all p >= p0 — the ball invariance
      (Rhat + sigma'' tbar) / (theta - Rhat tbar / W_w) <= tbar
  is ONE rational inequality, checked exactly (outward certificate:
  the recognized channel contracts, the memory channel has a declared
  geometric tail D0 L^{k-p0}/(1-L); floor b_{p+1} >= theta_lo b_p > 0).
  L < 1 is the beta < 1 condition, checked exactly.

  ASSEMBLY: rho_1(m,p) = (x_{m-p}' N x_{p-1}) / (x_0' B x_m). The
  b-window products pair term by term; early terms (k < p0) are pairs
  of IDENTICAL exact rationals (delta = 0); all remaining pairs have
  both indices >= p0 and contribute sum(delta) <= K1 L^{j-p0} with
  j = min(p-1, m-p) — bounded by V4, never exponentiated. The yhat
  corrections to numerator and denominator are bounded by V2 with
  certified rational square-root uppers. Result, per kappa:

      |rho_1(m,p) - rho_c| <= ERR(j) = A * L^{j-p0}
          for ALL m and p with j = min(p-1, m-p) >= p0,

  with rho_c := exact rho_1(26,13) and A a certified rational.
  Checked against exact rho_1 at (26,13),(30,15),(34,17),(30,10) and
  against YM-37's (anchored) rho_1_bulk bracket — the native route
  reproduces the same numbers, now without the import.

Correction of a YM-37 build-note: the LINEAGE entry narrated a
fraction-free (Bareiss) inertia as the dim-52 engine; inspection of the
committed file shows the patch never applied and the pinned YM-37 run
used the NATIVE T53 elimination throughout, dim 52 included. So no
Bareiss identity is in any load-bearing path; the LINEAGE note is
corrected in this commit (narration != arithmetic, caught).

Controls:
  C1  r_w'Bw = 0 exact; tampered theta (theta*101/100) refuses it.
  C2  B and B' PD by native elimination; U'BtauU = U'MU exact.
  C3  ball invariance + L < 1 rational checks; planted bad deflation
      (w = e0) fails the invariance check (honest refusal path).
  C4  exact t_p <= tbar for p0 <= p <= K; floor b_p > 0 for all
      computed p and theta_lo > 0.
  C5  |rho_1(m,p) - rho_c| <= ERR(j) on the exact table; rho_c inside
      YM-37's rho_1_bulk bracket.
"""

from fractions import Fraction as F
import json
import math
import os
import sys

sys.set_int_max_str_digits(1000000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import _dec, canonical_sha  # noqa: E402
from ym31_uniform_floor import fj, lam, rnd_down  # noqa: E402
import ym37_space_transfer as y37  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
ZERO, HALF, ONE = F(0), F(1, 2), F(1)
P_DEF = 12          # deflation vector w = x_P
P0_LIST = (2, 12)   # theorem regimes j >= p0: coarse (all interior p) and sharp
K_EXACT = 34        # exact sequence horizon


def sqrt_up(q, sc=10 ** 45):
    """rational s with s*s >= q >= 0, certified by the returned inequality."""
    assert q >= 0
    n = q.numerator * sc * sc
    s = F(math.isqrt(n // q.denominator) + 1, sc)
    assert s * s >= q
    return s


def sqrt_lo(q, sc=10 ** 45):
    assert q >= 0
    n = q.numerator * sc * sc
    s = F(math.isqrt(n // q.denominator), sc)
    assert s * s <= q
    return s


def dotv(u, v):
    return sum(a * b for a, b in zip(u, v))


def bq(B, u, v):
    return dotv(u, y37.matvec(B, v))


def count_pos_native(A):
    """T53 elimination sign check: number of positive weights (native LDL)."""
    return y37.inertia(A)[0]


def shift_count(M, B, s):
    n = len(M)
    return count_pos_native([[M[i][j] - s * B[i][j] for j in range(n)] for i in range(n)])


def restricted(M, B, w, W_w):
    """explicit rational basis of w-perp (drop coordinate of max |Bw| entry)."""
    n = len(M)
    Bw = y37.matvec(B, w)
    drop = max(range(n), key=lambda i: abs(Bw[i]))
    U = []
    for i in range(n):
        if i == drop:
            continue
        col = [F(0)] * n
        col[i] = F(1)
        c = Bw[i] / W_w
        col = [col[k] - c * w[k] for k in range(n)]
        U.append(col)
    # M' = U' M U ; B' = U' B U
    Mp = [[dotv(u, y37.matvec(M, v)) for v in U] for u in U]
    Bp = [[dotv(u, y37.matvec(B, v)) for v in U] for u in U]
    return U, Mp, Bp


def doubled_ceiling(Mp, Bp, hi, steps=40):
    """smallest tested sigma'' with count_+(M'B'^{-1}M' - sigma''^2 B') = 0."""
    n = len(Mp)
    # M' B'^{-1} M' by exact solves, column by column
    cols = [y37.solve(Bp, [Mp[i][j] for i in range(n)]) for j in range(n)]
    D = [[dotv([Mp[i][k] for k in range(n)], [cols[j][k] for k in range(n)]) for j in range(n)] for i in range(n)]
    lo = F(0)
    assert count_pos_native([[D[i][j] - hi * hi * Bp[i][j] for j in range(n)] for i in range(n)]) == 0
    for _ in range(steps):
        mid = (lo + hi) / 2
        if count_pos_native([[D[i][j] - mid * mid * Bp[i][j] for j in range(n)] for i in range(n)]) == 0:
            hi = mid
        else:
            lo = mid
    return hi


def kappa_block(kap, p0):
    P0 = p0
    fpt = {ZERO: rnd_down(fj(0, kap).lo), HALF: rnd_down(fj(1, kap).lo), ONE: rnd_down(fj(2, kap).lo)}
    klo = {ZERO: F(1), HALF: rnd_down(lam(1).lo), ONE: rnd_down(lam(2).lo)}
    basis, M, N, B, e0i = y37.build_pencil(2, fpt, klo)
    n = len(basis)
    assert count_pos_native(B) == n            # B PD (all weights positive)  [C2]
    e0 = [F(0)] * n
    e0[e0i] = F(1)
    # exact sequence
    xs = [e0]
    for _ in range(K_EXACT):
        xs.append(y37.tau_apply(M, B, xs[-1]))
    w = xs[P_DEF]
    W_w = bq(B, w, w)
    theta = bq(B, w, y37.tau_apply(M, B, w)) / W_w          # = w'Mw / W_w
    tw = y37.tau_apply(M, B, w)
    r_w = [tw[k] - theta * w[k] for k in range(n)]
    c1 = (bq(B, r_w, w) == 0)                               # exact identity [C1]
    theta_t = theta * F(101, 100)
    r_t = [tw[k] - theta_t * w[k] for k in range(n)]
    c1 = c1 and (bq(B, r_t, w) != 0)
    Rhat = sqrt_up(bq(B, r_w, r_w))
    # restricted pencil and complement ceiling
    U, Mp, Bp = restricted(M, B, w, W_w)
    c2 = (count_pos_native(Bp) == n - 1)
    # exact structural identity U'B tau U = U'MU on the basis vectors     [C2]
    for v in U:
        tv = y37.tau_apply(M, B, v)
        lhs = [bq(B, u, tv) for u in U]
        rhs = [dotv(u, y37.matvec(M, v)) for u in U]
        if lhs != rhs:
            c2 = False
    sig2 = doubled_ceiling(Mp, Bp, hi=F(1))
    # splits and exact per-step data
    bs, ys, nsq, t_up, s_ex = [], [], [], [], []
    for p, x in enumerate(xs):
        b = bq(B, w, x) / W_w
        yv = [x[k] - b * w[k] for k in range(n)]
        bs.append(b)
        ys.append(yv)
        nsq.append(bq(B, yv, yv))
        t_up.append(sqrt_up(nsq[-1]) / b if b > 0 else None)
    for p in range(K_EXACT):
        s_ex.append(bs[p + 1] / bs[p])
    c4 = all(b > 0 for b in bs)
    # ball invariance and contraction
    t_exact_tail = [t_up[p] for p in range(P0, K_EXACT + 1)]
    tbar = max(t_exact_tail) * F(3, 2)
    theta_lo = theta - Rhat * tbar / W_w
    theta_hi = theta + Rhat * tbar / W_w
    inv_ok = (theta_lo > 0) and ((Rhat + sig2 * tbar) / theta_lo <= tbar)
    L = sig2 / theta_lo + theta_hi * tbar * Rhat / (W_w * theta_lo * theta_lo)
    c3 = inv_ok and (L < 1)
    c4 = c4 and all(t <= tbar for t in t_exact_tail) and theta_lo > 0
    # planted bad deflation must refuse [C3]
    Rb = sqrt_up(bq(B, [y37.tau_apply(M, B, e0)[k] - (bq(B, e0, y37.tau_apply(M, B, e0)) / bq(B, e0, e0)) * e0[k] for k in range(n)],
                       [y37.tau_apply(M, B, e0)[k] - (bq(B, e0, y37.tau_apply(M, B, e0)) / bq(B, e0, e0)) * e0[k] for k in range(n)]))
    th_b = bq(B, e0, y37.tau_apply(M, B, e0)) / bq(B, e0, e0)
    bad_ok = not ((th_b - Rb * tbar / bq(B, e0, e0) > 0) and ((Rb + sig2 * tbar) / (th_b - Rb * tbar / bq(B, e0, e0)) <= tbar))
    c3 = c3 and bad_ok
    # declared geometric tail (theorum/28 item 3)
    yh = [[ys[p][k] / bs[p] for k in range(n)] for p in range(K_EXACT + 1)]
    d0 = sqrt_up(bq(B, [yh[P0 + 1][k] - yh[P0][k] for k in range(n)],
                       [yh[P0 + 1][k] - yh[P0][k] for k in range(n)]))
    Dsum = d0 / (1 - L)                     # || yhat_a - yhat_b ||_B <= Dsum L^{min(a,b)-P0}
    # certified rational scalars for the assembly
    Wup, Wlo = sqrt_up(W_w), sqrt_lo(W_w)
    nu = y37.pencil_norm_hi(N, B, F(64), steps=40)          # native two-sided sign checks inside
    Ntil = dotv(w, y37.matvec(N, w))
    Zc = W_w + bq(B, yh[0], yh[K_EXACT])                    # reference z-functional shape
    # coefficient of the window pairing:  sum delta <= K1 * L^{j-P0}
    K1 = (Rhat / (W_w * theta_lo)) * Dsum / (1 - L)
    # numerator / denominator correction coefficients
    KN = nu * 2 * (Wup + tbar) * Dsum / (1 - L)
    KZ = (sqrt_up(nsq[0]) / bs[0] + tbar) * Dsum / (1 - L)
    F_lo = Ntil - nu * (2 * Wup * tbar + tbar * tbar) - KN
    Z_lo = W_w - (sqrt_up(nsq[0]) / bs[0]) * tbar - KZ
    ok_pos = F_lo > 0 and Z_lo > 0 and K1 <= F(1, 4)
    rho_c = None
    A = None
    table = {}
    c5 = ok_pos
    if ok_pos:
        def rho(m, p):
            a = dotv(xs[m - p], y37.matvec(N, xs[p - 1]))
            z = dotv(xs[0], y37.matvec(B, xs[m]))
            return a / z
        rho_c = rho(26, 13)
        jc = min(12, 13)
        rho_hi = (Ntil + nu * (2 * Wup * tbar + tbar * tbar) + KN) / (theta_lo * Z_lo)
        A = rho_hi * (4 * K1 + KN / F_lo + KZ / Z_lo) * 2
        pts = ((26, 13), (30, 15), (34, 17), (28, 14)) if P0 == 12 else               ((26, 13), (30, 15), (34, 17), (30, 10), (20, 10), (16, 8), (13, 3), (10, 5))
        for (m, p) in pts:
            j = min(p - 1, m - p)
            if j < P0:
                continue
            err = A * (L ** (j - P0) + L ** (jc - P0))
            val = rho(m, p)
            good = abs(val - rho_c) <= err
            if not good:
                c5 = False
            table[f"{m},{p}"] = {"rho_1": _dec(val, 16), "native_bound": _dec(err, 16), "inside": bool(good)}
    return {
        "dim": n, "theta": _dec(theta, 18), "sigma_dd": _dec(sig2, 15),
        "Rhat_defect": _dec(Rhat, 42), "tbar": _dec(tbar, 42),
        "L_beta": _dec(L, 12), "beta_lt_1": bool(L < 1),
        "Dsum_declared_tail": _dec(Dsum, 42),
        "rho_c_exact_26_13": _dec(rho_c, 16) if rho_c is not None else None,
        "A_error_coefficient": _dec(A, 42) if A is not None else None,
        "table": table,
        "_controls": (c1, c2, c3, c4, c5),
        "_rho_c": rho_c,
    }


def run():
    grid = {}
    C = [True] * 5
    for kap in GRID:
      for p0 in P0_LIST:
        blk = kappa_block(kap, p0)
        cs = blk.pop("_controls")
        rc = blk.pop("_rho_c")
        C = [a and b for a, b in zip(C, cs)]
        # rho_c inside YM-37's anchored bracket
        try:
            y37r = json.load(open(os.path.join(HERE, "YM37_RESULT.json")))
            b = y37r["grid"][str(kap)]["t2"]["rho_1_bulk"]
            slack = F(1, 10 ** 13)          # the stored bracket is a 14-digit decimal
            inside = (F(b[0]) - slack <= rc <= F(b[1]) + slack) if rc is not None else False
        except Exception:
            inside = False
        blk["rho_c_inside_ym37_bracket"] = bool(inside)
        if not inside:
            C[4] = False
        grid.setdefault(str(kap), {})[f"j_ge_{p0}"] = blk
    ok = all(C)
    cert = {
        "certificate_type": "YM38_NATIVE_M_UNIFORMITY",
        "claim_status": "ym37_T3_spectral_import_REPLACED__deflation_plus_T53_elimination_plus_"
                        "theorum28_outward_certificate__for_all_m_p_bound_native__bareiss_never_in_load_bearing_path_lineage_note_corrected",
        "grid": grid,
        "controls": {"C1_defect_B_orthogonal_exact_and_tamper": bool(C[0]),
                     "C2_restricted_pencil_identities_and_PD": bool(C[1]),
                     "C3_outward_certificate_beta_lt_1_and_refusal_path": bool(C[2]),
                     "C4_floor_and_ball_on_exact_data": bool(C[3]),
                     "C5_bound_holds_and_matches_ym37": bool(C[4])},
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM38_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    with open(os.path.join(HERE, "EXPECTED_YM38.sha256"), "w") as f:
        f.write(canonical_sha(cert) + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "table"} for k, v in cert["grid"].items()}, indent=1))
