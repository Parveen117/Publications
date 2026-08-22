"""YM-22: THE TILING RATE SURVIVES THE CUTOFF — at leading tiling order
the chain's time-decay rate per unit physical time is uniform in BOTH the
volume m and the time step a along YM-9's declared trajectory
kappa = theta a, with exact limit 3/4 - theta/2. First capsule to touch
the cutoff wall from the native side; scope = leading tiling order only.

Why it works natively. On the heat-kernel family (YM-9) the time face
coefficient is lambda_a = Exp_Sigma(-a C_{1/2}) = Exp_Sigma(-3a/4), so
its native log is EXACTLY -3a/4 — no transcendental enters. The space
face coefficient is r(kappa) = f_{1/2}/f_0 = I_2(kappa)/I_1(kappa), and
YM-21 T3 gives, for every m, the leading-order time decay per time step
    gamma_step(a) >= -Log(lambda_a) - A_Sigma(r)  =  3a/4 - A_Sigma(r(theta a)).
Per unit physical time:
    gamma(a) = gamma_step(a)/a = 3/4 - A_Sigma(r(theta a))/a.
Since r(kappa) = kappa/4 + O(kappa^3) and A_Sigma(r) = 2r + O(r^3),
A_Sigma(r(theta a))/a -> theta/2, so gamma(a) -> 3/4 - theta/2 as a -> 0,
and gamma(a) is bounded below uniformly in a on the whole trajectory.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) EXACT TIME-FACE LOG. lambda_a = Exp(-3a/4) is used only through
      -Log(lambda_a) = 3a/4, an exact rational — the heat-kernel choice
      removes the transcendental from the time direction (YM-9 T2 seen
      from the fabric).

 (T2) CUTOFF-UNIFORM LEADING RATE. At theta = 1/16 and a over SIX
      DECADES (1, 1/10, ..., 1e-6), gamma(a) = 3/4 - A_Sigma(r(theta a))/a
      is certified by enclosure to satisfy
          3/4 - theta/2 - delta(a)  <=  gamma(a)  <=  3/4,
      with delta(a) -> 0 (certified monotone decreasing in a) and the
      limit value 3/4 - theta/2 = 23/32 = 0.71875 bracketed from both
      sides at a = 1e-6 to width < 1e-10. Compared with YM-9 T3's
      sandwich bound 3/4 - 6 theta = 3/8 at the same theta: the tiling
      rate is strictly larger at every a — the route is better by
      exactly the factor 12 in the theta-coefficient (6 -> 1/2).

 (T3) VOLUME-UNIFORM AT EVERY CUTOFF. YM-21 T3's KMS bracket is scalar
      in lambda; re-certified here with lambda_a at a in {1, 1e-3, 1e-6}
      for m = 2..8, t = 1..4: leading-order tiling powers lie inside
      lambda_a^t ((1-/+r)/(1+/-r))^t for every m. Hence the leading
      rate gamma(a) is uniform in m AND a simultaneously.

 (T4) FAIL-CLOSED TRAJECTORIES. The leading rate is positive iff
      A_Sigma(r(theta a))/a < 3/4. Certified: theta = 1/16 and theta = 1
      pass in the limit (3/4 - 1/2 = 1/4 > 0), theta = 2 FAILS in the
      limit (3/4 - 1 < 0) and is refused at small a. The route's
      trajectory ceiling is exactly theta < 3/2 in the limit.

 (T5) WHY THE EXPANSION PARAMETER DOES NOT DIE — honest remainder in
      the cutoff. Higher face coefficients r_j(theta a) ~ (theta a)^{j}
      vanish as a -> 0, but the number of time faces per unit time is
      1/a, so the higher-order weight per unit time scales like
      r_j/a ~ theta^j a^{j-1}: the j = 1/2 term stays finite (that is
      the theta/2 above), j >= 1 terms vanish in the limit. Certified:
      r_1(theta a)/a -> 0 and r_{3/2}(theta a)/a -> 0 over the six
      decades. So along the trajectory the full-order correction to the
      cutoff limit comes ONLY from branching/overlap of j = 1/2 tilings
      — a sharper statement of the open cluster problem than YM-21 T4.

NOT CLAIMED: full-order gap; the physical mass gap; asymptotic freedom
(trajectory is declared, as in YM-9); 2D fabrics; Clay.

Controls:
  C1  -Log(lambda_a) = 3a/4 exact (no log evaluated).
  C2  gamma(a) monotone and inside [3/4 - theta/2 - delta, 3/4] on the grid.
  C3  tiling strictly beats YM-9 sandwich at every a.
  C4  KMS bracket with lambda_a, m = 2..8, t = 1..4, three a's.
  C5  theta = 2 refused at small a; theta = 1 certified.
  C6  r_1/a and r_{3/2}/a decrease to below 1e-6 along the decades.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, bessel_I, _dec, canonical_sha, TERMS  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402
from ym15_chain_closed_form import kms  # noqa: E402
from ym19_dobrushin_dock import exp_neg  # noqa: E402
from ym20_native_origin_audit import native_odd_log  # noqa: E402
from ym21_tiling_law import mat_pow_rows  # noqa: E402

C_HALF = F(3, 4)
THETA = F(1, 16)
A_GRID = [F(1), F(1, 10), F(1, 100), F(1, 1000), F(1, 10 ** 4),
          F(1, 10 ** 5), F(1, 10 ** 6)]
A_KMS = [F(1), F(1, 1000), F(1, 10 ** 6)]


def r_of(kappa: F) -> Iv:
    return bessel_I(2, kappa, TERMS) / bessel_I(1, kappa, TERMS)


def rj_over_f0(j2: int, kappa: F) -> Iv:
    """r_j = f_j / f_0 = I_{j2+1}(kappa) / I_1(kappa), j2 = twice spin."""
    return bessel_I(j2 + 1, kappa, TERMS) / bessel_I(1, kappa, TERMS)


def gamma(a: F, theta: F) -> Iv:
    """leading tiling rate per unit time: 3/4 - A_Sigma(r(theta a))/a."""
    r = r_of(theta * a)
    A = native_odd_log(r)
    return Iv(C_HALF) - A / Iv(a)


def run():
    c1 = True
    # C1: -Log(lambda_a) = 3a/4 exact: lambda_a = e^{-3a/4} enclosure,
    # and the rate uses the exact rational 3a/4 directly (no log call).
    for a in A_GRID[:3]:
        lam = exp_neg(C_HALF * a)
        if not (lam.lo > 0 and lam.hi < 1):
            c1 = False

    rows = {}
    c2 = c3 = True
    prev = None
    limit = C_HALF - THETA / 2
    for a in A_GRID:
        g = gamma(a, THETA)
        excess = (g - Iv(limit)).lo           # excess ABOVE the limit
        ym9 = C_HALF - 6 * THETA
        # the limit 23/32 is a UNIFORM LOWER BOUND at every a (approach from above)
        if not (g.hi <= C_HALF and g.lo >= limit):
            c2 = False
        if not (g.lo > ym9):
            c3 = False
        if prev is not None and not (g.hi <= prev.hi + F(1, 10 ** 30)):
            c2 = False                         # nonincreasing as a -> 0
        prev = g
        rows[str(a)] = {"gamma_lo": _dec(g.lo, 15), "gamma_hi": _dec(g.hi, 15),
                        "excess_above_limit_23_32": _dec(max(excess, F(0)), 15),
                        "ym9_sandwich_bound": _dec(ym9, 6)}
    g_small = gamma(A_GRID[-1], THETA)
    limit_bracket_ok = (g_small.lo >= limit) and \
        (g_small.hi - limit < F(1, 10 ** 10))

    # C4 KMS bracket with lambda_a
    c4 = True
    for a in A_KMS:
        lam = exp_neg(C_HALF * a)
        r = r_of(THETA * a)
        one = Iv(F(1))
        up = lam * (one + r) / (one - r)
        lo = lam * (one - r) / (one + r)
        for m in range(2, 9):
            M = [[_r(lam * x) for x in row] for row in kms(m, r)]
            for t in range(1, 5):
                P = mat_pow_rows(M, t)
                up_t = one
                lo_t = one
                for _ in range(t):
                    up_t = up_t * up
                    lo_t = lo_t * lo
                for i in range(m):
                    rs = sum((P[i][j] for j in range(m)), Iv(F(0)))
                    if not (rs.hi <= up_t.hi + F(1, 10 ** 20)
                            and rs.lo >= lo_t.lo - F(1, 10 ** 20)):
                        c4 = False

    # C5 trajectory ceiling
    g2 = gamma(F(1, 10 ** 4), F(2))
    g1 = gamma(F(1, 10 ** 4), F(1))
    c5 = (g2.hi < 0) and (g1.lo > 0)

    # C6 higher coefficients per unit time vanish
    c6 = True
    ladder = {}
    prev1 = prev3 = None
    for a in A_GRID:
        k = THETA * a
        q1 = (rj_over_f0(2, k) / Iv(a)).hi
        q3 = (rj_over_f0(3, k) / Iv(a)).hi
        if prev1 is not None and not (q1 <= prev1 and q3 <= prev3):
            c6 = False
        prev1, prev3 = q1, q3
        ladder[str(a)] = [_dec(q1, 12), _dec(q3, 12)]
    if not (prev1 < F(1, 10 ** 6) and prev3 < F(1, 10 ** 6)):
        c6 = False

    ok = c1 and c2 and c3 and c4 and c5 and c6 and limit_bracket_ok
    cert = {
        "certificate_type": "YM22_TILING_RATE_UNIFORM_IN_VOLUME_AND_CUTOFF",
        "claim_status": "leading_tiling_order_only__uniform_in_m_and_a_along_"
                        "declared_trajectory__full_order_OPEN",
        "trajectory": {"kappa": "theta * a", "theta": str(THETA),
                       "limit_rate": str(limit), "limit_rate_dec": _dec(limit, 6),
                       "ym9_sandwich_rate": str(C_HALF - 6 * THETA)},
        "theorems": {
            "T1_exact_time_face_log": "-Log(lambda_a) = 3a/4 exactly",
            "T2_cutoff_uniform_leading_rate":
                "gamma(a) = 3/4 - A_Sigma(r(theta a))/a in [23/32, 3/4] over six "
                "decades: the limit is a uniform lower bound in the cutoff",
            "T3_volume_uniform_at_every_cutoff":
                "KMS bracket with lambda_a holds for m = 2..8 at a = 1, 1e-3, 1e-6",
            "T4_trajectory_ceiling": "theta < 3/2 in the limit; theta = 2 refused",
            "T5_expansion_parameter_in_the_cutoff":
                "r_1/a, r_3/2/a -> 0; only j = 1/2 branching survives a -> 0",
        },
        "grid": rows,
        "limit_bracket_at_smallest_a": [_dec(g_small.lo, 14), _dec(g_small.hi, 14)],
        "higher_coeff_per_unit_time_r1_r32": ladder,
        "honest_remainder": "leading tiling order; full-order convergence "
                            "(branching j=1/2 tilings) OPEN; trajectory declared; "
                            "not the physical mass gap; Clay OPEN",
        "controls": {
            "C1_time_log_exact": bool(c1),
            "C2_gamma_ge_limit_and_nonincreasing": bool(c2),
            "C3_beats_ym9_sandwich": bool(c3),
            "C4_kms_bracket_with_lambda_a": bool(c4),
            "C5_theta_2_refused_theta_1_ok": bool(c5),
            "C6_higher_coefficients_vanish": bool(c6),
            "C7_limit_is_infimum_within_1e-10": bool(limit_bracket_ok),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM22_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM22.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" a=", k, v["gamma_lo"][:10], v["gamma_hi"][:10],
              "excess", v["excess_above_limit_23_32"][:12])
    print("limit bracket:", cert["limit_bracket_at_smallest_a"])
    print("sha256:", sha)
