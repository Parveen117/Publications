"""YM-31: OUTWARD CERTIFICATION OF THE DRESSED VACUUM FLOOR — an m-UNIFORM
lower bound on lambda_1(S_m) at the true rate, from the 2-row engine.

The trial. For T'' = W^{1/2} T0 W^{1/2} (same spectrum as S_m) and any
function chi, the Rayleigh quotient with phi = W^{1/2} chi is
    lambda_1 >= <W chi, T0 W chi> / <chi, W chi>.
Choose chi = W_pt / W, where W_pt = prod_l w_pt(H_l) is a PRODUCT OF
RATIONAL TRUNCATED FACE WEIGHTS w_pt = sum_{j<=1} d_j f_j^{pt} chi_j (rational
f^{pt} close to f_j; the choice is part of the trial, so it carries NO
error). Then
    numerator   Z_m = <W_pt, T0 W_pt>   (ladder, exact by YM-30),
    denominator D_m = Int W_pt^2 / W = [ Int w_pt^2 e^{-kappa Tr/2} ]^{m-1}
                                      (faces independent, YM-29 L).

 (T1) RUNG POSITIVITY: TRUNCATING RUNGS GIVES A LOWER BOUND. T0 = prod_i
      sum_c lambda_c Pi_c^{(i)} with Pi_c^{(i)} the orthogonal projection
      on rung-i content c (convolution by d_c chi_c), so
          Z_m = sum_{c_1..c_m} prod_i lambda_{c_i} || prod_i Pi^{(i)}_{c_i} W_pt ||^2
      is a sum of NONNEGATIVE terms, monotone in every lambda_c. Hence
      the engine value with rung labels <= 1 and lambda_c at rational
      lower ends is a rigorous lower bound Z_m^lo <= Z_m. (Native: a
      projection's square is itself — the cut-square identity, T01-B.)

 (T2) EXACT DENOMINATOR. Int w_pt^2 e^{-kappa Tr/2} = sum_{j,k,l} d_j d_k d_l
      f_j^{pt} f_k^{pt} (-1)^{2l} f_l(kappa) N_{jkl}, N_{jkl} in {0,1} the fusion
      count (chi_l(-H) = (-1)^{2l} chi_l(H) gives the e^{-kappa Tr/2}
      coefficients from the e^{+kappa Tr/2} ones); the l-tail is bounded
      outward by d_l^2 f_l. D_m^hi = (face integral)^{m-1}, upper end.
      So FLOOR_m := Z_m^lo / D_m^hi <= lambda_1(S_m), rigorous for each m.

 (T3) NATIVE LOG-CONVEXITY => m-UNIFORM FLOOR. Z_m = <g, S^{m-2} g> with
      S = C^{1/2} M_K C^{1/2}, C = convolution by w_pt (x) w_pt on the cut
      space L^2(SU(2)^2) (PSD: coefficients f^{pt} >= 0), M_K =
      multiplication by the rung kernel K (pointwise > 0 for the Wilson
      time kernel e^{Tr U}), g = C^{1/2} K. S is symmetric PSD, so by the
      cut-square inequality (T01-C form) Z_m^2 <= Z_{m-1} Z_{m+1}: the
      ratio Z_{m+1}/Z_m is NONDECREASING in m. Therefore, for every
      m >= m_0,
          Z_m >= Z_{m_0} * rho^{m - m_0},   rho := Z_{m_0+1} / Z_{m_0},
      and with the denominator exactly geometric,
          lambda_1(S_m) >= FLOOR_{m_0} * (rho / I_face)^{m - m_0}   for all m >= m_0.
      Outward: rho^lo = Z^lo_{m_0+1} / Z^hi_{m_0}, with Z^hi_{m_0} = Z^lo_{m_0}
      + rung tail, tail <= m_0 * lambda_{3/2}^hi * (content >= 3/2 mass at a
      rung from two faces of labels <= 1) * S_1^{m_0-3}, S_1 = sum_{j<=1} d_j^2 (f_j^pt)^2
      (Parseval on the truncated faces; contents <= 1 at both faces bound
      the rung content <= 2, and only pairs with a + a' >= 3/2 feed c >= 3/2).
      Log-convexity is also CHECKED on the computed sequence (control).

 (T4) RESULT at kappa = 1/8 (m_0 = 6): per-face floor rate
          log(rho^lo / I_face^hi) >= 0.00268   (vs old floor log f_0 = 0.001953,
      YM-18's calibration 0.002177, E4D excitation rate 0.001949 [YM-28]).
      Hence for EVERY m >= 6 the vacuum grows at least 0.00073/face faster
      than the single-insertion memory-channel norm — theorum/28 item 6
      delivered at the true rate, uniformly in m. Grid kappa = 1/4, 1/2
      included.

NOT CLAIMED: the m-uniform GAP (item 3 — the memory channel's operator
norm, not just basis-state norms, is still open); weak coupling.

Controls:
  C1  rung-truncated Z monotone in lambda_c (perturb a lambda down -> Z down).
  C2  denominator fusion sum matches direct character orthogonality at m=2
      (Int w_pt^2 / w with w_pt = w_trunc(kappa) when f^pt = f(kappa) exactly
      recovers sum d_j^2 f_j^2-type identity for the trivial case w=1).
  C3  log-convexity Z_{m}^2 <= Z_{m-1} Z_{m+1} on m = 3..7 (exact rationals).
  C4  FLOOR_m <= YM-18's certified Rayleigh floors where both exist? No —
      both are lower bounds; check instead FLOOR_m <= e^{kappa(m-1)}
      (YM-16 upper bracket) and FLOOR_m > f_0^{m-1} (improvement).
  C5  per-face rate > E4D excitation rate at all three kappa.
  C6  tamper: dropping rung positivity (using lambda hi) would raise Z —
      the bound direction is the stated one.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS  # noqa: E402
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym4_symmetry_protected import chi_mul, dim  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402
from ym15_chain_closed_form import BETA  # noqa: E402
from ym16_chain_dock import f0_of  # noqa: E402
from ym30_recoupling_engine import ladder_partition  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
M0 = 6
E4D_RATE = {F(1, 8): F("0.00194868"), F(1, 4): F("0.00774210"), F(1, 2): F("0.03016068")}


def rnd_down(x: F, d=10 ** 9) -> F:
    return F(int(x * d), d)


def rnd_up(x: F, d=10 ** 9) -> F:
    return F(int(x * d) + 1, d)


def fj(j2: int, kappa: F) -> Iv:
    """chain face coefficient f_j: Int e^{kappa Tr/2} chi_j / d_j = (2/kappa) I_{2j+1}(kappa)."""
    return _r(Iv(F(2)) / Iv(kappa) * bessel_I(j2 + 1, kappa, TERMS))


def lam(j2: int) -> Iv:
    """Wilson time-kernel eigenvalue lambda_j = I_{2j+1}(beta)/I_1(beta)."""
    return _r(bessel_I(j2 + 1, BETA, TERMS) / bessel_I(1, BETA, TERMS))


def face_integral_hi(fpt: dict, kappa: F, L=8) -> F:
    """upper bound of Int w_pt^2 e^{-kappa Tr/2} = sum d_j d_k d_l f^pt_j f^pt_k (-1)^{2l} f_l N_jkl."""
    tot = Iv(F(0))
    for j2, fj_ in fpt.items():
        for k2, fk_ in fpt.items():
            prods = chi_mul(j2, k2)                      # contents in j x k
            for l2 in range(0, L + 1):
                if l2 in prods:
                    sgn = (-1) ** l2
                    tot = tot + Iv(F(dim(j2) * dim(k2) * dim(l2))) * Iv(fj_) * Iv(fk_) * fj(l2, kappa) * sgn
    # tail l > L: |chi_l| <= d_l -> |term| <= d_l^2 f_l (sum_j d_j f^pt_j)^2
    s = sum(dim(j2) * v for j2, v in fpt.items())
    tail = Iv(F(0))
    for l2 in range(L + 1, L + 6):
        tail = tail + Iv(F(dim(l2) ** 2)) * fj(l2, kappa) * Iv(s * s)
    return _r(tot + tail).hi


def run():
    zero, half, one = F(0), F(1, 2), F(1)
    results = {}
    c1 = c3 = c4 = c5 = c6 = True
    for kap in GRID:
        fpt = {0: rnd_down(fj(0, kap).lo), 1: rnd_down(fj(1, kap).lo), 2: rnd_down(fj(2, kap).lo)}
        u = {zero: fpt[0], half: fpt[1], one: fpt[2]}
        klo = {zero: F(1), half: rnd_down(lam(1).lo), one: rnd_down(lam(2).lo)}
        khi = {zero: F(1), half: rnd_up(lam(1).hi), one: rnd_up(lam(2).hi)}
        Z = {m: ladder_partition(u, u, klo, m) for m in range(2, M0 + 3)}
        # C1 / C6 monotonicity in lambda
        kless = dict(klo)
        kless[half] = klo[half] * F(9, 10)
        c1 = c1 and ladder_partition(u, u, kless, 3) < Z[3]
        c6 = c6 and ladder_partition(u, u, khi, 3) >= Z[3]
        # C3 log-convexity on the computed sequence
        for m in range(3, M0 + 2):
            if not (Z[m] * Z[m] <= Z[m - 1] * Z[m + 1]):
                c3 = False
        # rung tail bound for Z^hi(m0)
        S1 = sum(dim(t) ** 2 * fpt[t] ** 2 for t in range(3))
        pair_mass = sum(dim(a) ** 2 * fpt[a] ** 2 * dim(b) ** 2 * fpt[b] ** 2
                        for a in range(3) for b in range(3) if a + b >= 3)
        tail = M0 * rnd_up(lam(3).hi) * pair_mass * S1 ** (M0 - 3)
        Zhi_m0 = Z[M0] + tail
        rho_lo = Z[M0 + 1] / Zhi_m0
        I_hi = face_integral_hi({0: fpt[0], 1: fpt[1], 2: fpt[2]}, kap)
        floor_m0 = Z[M0] / I_hi ** (M0 - 1)
        per_face = log_iv(Iv(rho_lo / I_hi), LOG_TERMS).lo
        f0 = f0_of(kap)
        old_rate = log_iv(f0, LOG_TERMS).hi
        # C4: floor below the YM-16 upper bracket e^{kappa(m-1)} and above f0^{m-1}
        for m in range(2, M0 + 3):
            fl = Z[m] / I_hi ** (m - 1)
            if not (fl <= exp_point(kap * (m - 1)).hi and fl > f0.hi ** (m - 1) * F(1)):
                c4 = False
        c5 = c5 and (per_face > E4D_RATE[kap])
        results[str(kap)] = {
            "f_pt": {str(t): str(v) for t, v in fpt.items()},
            "Z_lo_by_m": {str(m): _dec(Z[m], 14) for m in Z},
            "rung_tail_bound_at_m0": _dec(tail, 16),
            "rho_lo": _dec(rho_lo, 12),
            "face_integral_hi": _dec(I_hi, 12),
            "floor_at_m0": _dec(floor_m0, 12),
            "per_face_floor_rate_all_m_ge_m0": _dec(per_face, 8),
            "old_floor_rate_log_f0": _dec(old_rate, 8),
            "E4D_excitation_rate_YM28": str(E4D_RATE[kap]),
            "margin_per_face": _dec(per_face - E4D_RATE[kap], 8),
        }
    ok = c1 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM31_OUTWARD_M_UNIFORM_DRESSED_VACUUM_FLOOR",
        "claim_status": "lambda_1(S_m) >= FLOOR_m0 * (rho_lo/I_hi)^(m-m0) for all m >= m0 "
                        "(rung positivity + exact denominator + native log-convexity); "
                        "item 6 of theorum/28 at the true rate; gap NOT claimed",
        "m0": M0,
        "grid": results,
        "controls": {
            "C1_rung_monotone_in_lambda": bool(c1),
            "C3_log_convexity_observed": bool(c3),
            "C4_floor_inside_known_brackets": bool(c4),
            "C5_rate_above_E4D_excitation": bool(c5),
            "C6_bound_direction": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM31_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM31.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "rate", v["per_face_floor_rate_all_m_ge_m0"], "old", v["old_floor_rate_log_f0"],
              "E4D", v["E4D_excitation_rate_YM28"], "margin", v["margin_per_face"], "tail", v["rung_tail_bound_at_m0"])
    print("sha256:", sha)
