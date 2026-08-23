"""YM-40: THE CONTENT LADDER ON THE HEAT-KERNEL TRAJECTORY — the channel
gaps become the Casimir ladder in the cutoff direction, with certified
brackets uniform in m.

YM-39's named next step. Setup = YM-22's family: heat-kernel TIME kernel
lambda_c(a) = Exp_Sigma(-a C_c) (C_1/2 = 3/4, C_1 = 2; rounded-down
rationals via YM-13's compound-interest bounds, the audited native
exp), Wilson space faces at kappa = theta a, theta = 1/16, over
a in {1, 1/4, 1/16, 1/64}. Machinery = YM-37/38/39 verbatim per a: the
same even-sector pencil, the same deflation, the same outward
certificate; channel matrices from ym39.build_N. Native log throughout
(F00-G odd series, ym1.log_iv). No new verdict types.

 (T1) Superselection persists on the whole family: the mixed
      half-integer/integer insertion matrix is EXACTLY zero at every a
      (the centre argument is a-independent; re-verified entry by entry).

 (T2) Per-a, per-channel bulk ratios rho_c(a) with the YM-38 m-uniform
      band, and certified physical rates
          gamma_c(a) := -Log(rho_c(a) -+ band) / a,
      two-sided, valid for EVERY m and interior p.

 (T3) CASIMIR LADDER. Certified for every grid a:
          Delta_gamma(a) := gamma_1(a) - gamma_1/2(a) >= 1,
      and a two-sided PINCH onto the Casimir difference:
          |Delta_gamma(a) - 5/4| <= eps(a),
      with eps(a) certified and strictly decreasing along the grid
      (6.4e-4 -> 1.6e-4 -> 7.0e-5 -> 2.0e-5). So in the cutoff
      direction the content-1 channel does not merely stay heavier
      (YM-39): its EXCESS rate is pinned onto the Casimir gap. Build
      note (draft refused by the arithmetic, YM-22's lesson again in
      mirror image): the first draft claimed monotone approach from
      BELOW; the certified brackets cross 5/4 after a = 1 and approach
      from above — the trajectory correction is signed, not one-sided.
      Only the pinch and the floor are claimed. Likewise gamma_1/2(a)
      itself crosses C_1/2 = 3/4 (from below at a = 1 to slightly above
      at a = 1/64): the dressing term changes sign along the family.

 (T4) Ladder collapse tie, quantified on the same grid: certified
      brackets for r_1(theta a)/r_1/2(theta a) decreasing toward 0 —
      the face-side statement (YM-22) next to the spectrum-side
      statement (T3), same a's.

Boundary (honest): all statements are about the order-1 dressed
two-point ratio on the carrier (compressed, lower-bound side), NOT the
operator gap; time kernel truncated at contents <= 1, so no chi_3/2
Casimir claim (C_3/2 = 15/4 is outside the declared kernel; the 3/2 row
is omitted here rather than reported as if fundamental); trajectory
DECLARED (theta = 1/16), not derived; E4D-C OPEN.

Controls:
  C1  mixed matrix zero at every a; same-parity pair nonzero.
  C2  outward certificate (ball + beta = L < 1) at every a; channel
      positivity F_lo, Z_lo > 0.
  C3  every exact table row inside its band (spot m,p per a).
  C4  gamma brackets: lower <= upper, and log_iv consistency
      (Exp_Sigma(-a gamma_hi) <= rho - band <= rho + band <=
       Exp_Sigma(-a gamma_lo)) checked with the independent exp route.
  C5  T3 inequalities exact: floor Delta_gamma >= 1 at every a, the
      pinch |Delta_gamma - 5/4| <= eps(a) with eps strictly decreasing;
      T4 ratio brackets strictly decreasing.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(1000000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import _dec, canonical_sha, log_iv, Iv  # noqa: E402
from ym19_dobrushin_dock import exp_neg  # noqa: E402
from ym31_uniform_floor import fj, rnd_down  # noqa: E402
import ym37_space_transfer as y37  # noqa: E402
import ym38_native_uniformity as y38  # noqa: E402
import ym39_content_ladder as y39  # noqa: E402

ZERO, HALF, ONE = F(0), F(1, 2), F(1)
THETA = F(1, 16)
A_GRID = [F(1), F(1, 4), F(1, 16), F(1, 64)]
C_HALF, C_ONE = F(3, 4), F(2)
K_EXACT = 34
LOG_TERMS = 60


def rdn(x, d=10 ** 30):
    return F((x.numerator * d) // x.denominator, d)


def rup(x, d=10 ** 30):
    return F((x.numerator * d) // x.denominator + 1, d)


def gamma_bracket(rho_ref, band, a):
    # YM-23 perf lesson: round outward BEFORE odd-series logs
    lo_r, hi_r = rdn(rho_ref - band), rup(rho_ref + band)
    assert 0 < lo_r and hi_r < 1
    lg = log_iv(Iv(lo_r, hi_r), LOG_TERMS)
    g_lo, g_hi = rdn(-lg.hi / a, 10 ** 12), rup(-lg.lo / a, 10 ** 12)
    # C4: independent exp route must re-enclose the ratio interval
    e_hi = exp_neg(a * g_lo)     # >= true e^{-a g_lo} >= hi_r
    e_lo = exp_neg(a * g_hi)     # <= ... <= lo_r
    assert e_lo.lo <= lo_r and hi_r <= e_hi.hi
    return g_lo, g_hi


def run():
    grid = {}
    C = [True] * 5
    prev_dlo = None
    prev_ratio_hi = None
    for a in A_GRID:
        kap = THETA * a
        fpt = {ZERO: rnd_down(fj(0, kap).lo), HALF: rnd_down(fj(1, kap).lo), ONE: rnd_down(fj(2, kap).lo)}
        klo = {ZERO: F(1), HALF: rnd_down(exp_neg(C_HALF * a).lo), ONE: rnd_down(exp_neg(C_ONE * a).lo)}
        basis, M, Nhalf, B, e0i = y37.build_pencil(2, fpt, klo)
        n = len(basis)
        e0 = [F(0)] * n
        e0[e0i] = F(1)
        xs = [e0]
        for _ in range(K_EXACT):
            xs.append(y37.tau_apply(M, B, xs[-1]))
        w = xs[12]
        W_w = y38.bq(B, w, w)
        tw = y37.tau_apply(M, B, w)
        theta_r = y38.bq(B, w, tw) / W_w
        r_w = [tw[k] - theta_r * w[k] for k in range(n)]
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
        tbar = max(t_up[y39.P0:]) * F(3, 2)
        th_lo = theta_r - Rhat * tbar / W_w
        th_hi = theta_r + Rhat * tbar / W_w
        inv_ok = th_lo > 0 and (Rhat + sig2 * tbar) / th_lo <= tbar
        L = sig2 / th_lo + th_hi * tbar * Rhat / (W_w * th_lo * th_lo)
        C[1] = C[1] and inv_ok and L < 1 and all(t <= tbar for t in t_up[y39.P0:])
        yhd = [[(xs[y39.P0 + 1][k] - bs[y39.P0 + 1] * w[k]) / bs[y39.P0 + 1]
                - (xs[y39.P0][k] - bs[y39.P0] * w[k]) / bs[y39.P0] for k in range(n)]][0]
        Dsum = y38.sqrt_up(y38.bq(B, yhd, yhd)) / (1 - L)
        # T1 superselection at this a
        Nmix = y39.build_N(basis, klo, HALF, ONE)
        C[0] = C[0] and all(x == 0 for row in Nmix for x in row)
        Nsame = y39.build_N(basis, klo, HALF, F(3, 2))
        C[0] = C[0] and any(x != 0 for row in Nsame for x in row)
        # T2 channels 1/2 and 1
        entry = {"kappa": str(kap)}
        gammas = {}
        for c, name in ((HALF, "half"), (ONE, "one")):
            Nc = Nhalf if c == HALF else y39.build_N(basis, klo, c, c)
            blk = y39.channel_block(M, B, Nc, xs, w, W_w, theta_r, r_w, Rhat, sig2, tbar, L, Dsum, nsqs[0], bs[0])
            if blk is None or not blk["ok_rows"]:
                C[1] = C[1] and blk is not None
                C[2] = C[2] and (blk is not None and blk["ok_rows"])
                continue
            g_lo, g_hi = gamma_bracket(blk["rho_ref_26_13"], blk["band_all_interior"], a)
            C[3] = C[3] and g_lo <= g_hi
            gammas[name] = (g_lo, g_hi)
            entry[f"rho_{name}"] = _dec(blk["rho_ref_26_13"], 16)
            entry[f"band_{name}"] = _dec(blk["band_all_interior"], 20)
            entry[f"gamma_{name}"] = [_dec(g_lo, 12), _dec(g_hi, 12)]
        # T3 Casimir ladder
        d_lo = gammas["one"][0] - gammas["half"][1]
        d_hi = gammas["one"][1] - gammas["half"][0]
        floor_ok = d_lo >= 1
        eps = max(abs(d_lo - F(5, 4)), abs(d_hi - F(5, 4)))
        pinch_ok = True if prev_dlo is None else eps < prev_dlo
        C[4] = C[4] and floor_ok and pinch_ok
        entry["Delta_gamma"] = [_dec(d_lo, 12), _dec(d_hi, 12)]
        entry["Delta_floor_ge_1"] = bool(floor_ok)
        entry["pinch_eps_to_5_4"] = _dec(eps, 12)
        entry["pinch_eps_strictly_decreasing"] = bool(pinch_ok)
        prev_dlo = eps
        # T4 face-ladder collapse on the same grid
        rr = (fj(2, kap) / fj(1, kap))
        entry["r1_over_rhalf"] = [_dec(rr.lo, 12), _dec(rr.hi, 12)]
        if prev_ratio_hi is not None:
            C[4] = C[4] and rr.hi < prev_ratio_hi
        prev_ratio_hi = rr.hi
        entry["L_beta"] = _dec(L, 14)
        grid[str(a)] = entry
    ok = all(C)
    cert = {
        "certificate_type": "YM40_HEAT_KERNEL_CONTENT_LADDER_CASIMIR",
        "claim_status": "superselection_on_family_EXACT__channel_rates_bracketed_m_uniform__"
                        "Delta_gamma_floor_1_and_two_sided_pinch_onto_Casimir_5_4__faces_collapse__"
                        "compressed_order1_only__trajectory_declared__E4DC_OPEN",
        "theta": str(THETA),
        "grid": grid,
        "controls": {"C1_superselection_every_a": bool(C[0]),
                     "C2_outward_certificate_every_a": bool(C[1]),
                     "C3_tables_inside": bool(C[2]),
                     "C4_gamma_brackets_and_exp_crosscheck": bool(C[3]),
                     "C5_ladder_floor_pinch_and_face_collapse": bool(C[4])},
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM40_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    with open(os.path.join(HERE, "EXPECTED_YM40.sha256"), "w") as f:
        f.write(canonical_sha(cert) + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print("a =", k, "gamma_half", v.get("gamma_half"), "gamma_one", v.get("gamma_one"),
              "Delta", v.get("Delta_gamma"), "r1/rhalf", v.get("r1_over_rhalf"))
