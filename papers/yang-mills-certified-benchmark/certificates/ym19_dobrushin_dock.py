"""YM-19: THE DOBRUSHIN DOCK — first certified VOLUME-UNIFORM gap of the
bounded-overlap chain S_m at strong coupling, via Dobrushin uniqueness
instead of the Kotecky-Preiss expansion; plus a STANDING CORRECTION to
YM-17 T3.

STANDING CORRECTION (YM-17 T3). YM-17 reduced the chain's volume
problem to the "vacuum-tracking inequality" lambda_1(T) vs
sigma_1(X)sigma_1(Y). That reduction is NOT sufficient: the ratio
lambda_1(T)/(sigma_1(X)sigma_1(Y)) is generically exponentially small in
m (the interleaved free energy is not the sum of the halves'), so even
perfect tracking leaves the singular-value Weyl route exponentially
loose. YM-18's measurement (per-bridge rate settles) is consistent with
this and is retained as calibration. The gap must be controlled through
CORRELATION DECAY, not through vacuum weights. This capsule does that.

SETTING. Chain S_m with heat-kernel free transfer (YM-9 family, time
step a) and bridges e^{(kappa/2)Tr(A_i A_{i+1}^{-1})}. The symmetric
transfer operator T = K_a^{1/2 (x) m} m_kappa K_a^{1/2 (x) m} is the
transfer matrix of the nearest-neighbour Gibbs field on the strip
{1..m} x Z with bond weights
    time bonds : w_t(x,y) = K_a(x y^{-1}),  K_a = sum_j d_j e^{-a C_j} chi_j
    space bonds: w_s(x,y) = exp[(kappa/2) Tr(x y^{-1})].
(K_a is the heat kernel normalised to mean 1; it is a positive class
function for every a > 0.)

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) DOBRUSHIN COEFFICIENT BOUND (proved here, elementary). For a
      bond weight w with delta := sup w / inf w, the conditional law at
      a site i given its neighbours changes, when ONE neighbour j is
      changed, by total variation at most
          C_ij  <=  1 - 1/delta^2 .
      Proof: p ∝ R w(.,x_j), q ∝ R w(.,x_j'); pointwise
      w(.,x_j)/w(.,x_j') in [1/delta, delta] and Z_q/Z_p >= 1/delta, so
      p/q >= 1/delta^2 and TV(p,q) = Int q (1 - p/q)_+ <= 1 - 1/delta^2.
      Exact deltas:  delta_s = e^{2 kappa}  (|Tr U| <= 2, attained);
      delta_t <= (1+S_a)/(1-S_a), S_a = sum_{j>=1/2} d_j^2 e^{-a C_j}
      (K_a(1) = 1 + S_a, K_a >= 1 - S_a since |chi_j| <= d_j), with the
      j-tail certified by geometric domination.

 (T2) DOBRUSHIN CONDITION ON THE STRIP, UNIFORM IN m. Every site has at
      most 2 time- and 2 space-neighbours, so
          alpha(a,kappa) := 2(1 - 1/delta_t^2) + 2(1 - 1/delta_s^2)
      bounds sup_i sum_j C_ij for EVERY m (boundary sites have fewer
      neighbours). alpha < 1 is certified on the grid
      (a, kappa) in {(6,1/16), (8,1/8), (12,1/8)} and REFUSED at
      (6,1/8), (4,1/16) and (10,1/4) (fail-closed). Note the space term
      alone forces kappa < (log 2)/4 = 0.173 for ANY a: the route's
      coupling ceiling is exact and small.

 (T3) UNIFORM GAP (pinned classical anchor, CIRC-1 style as YM-8).
      ANCHOR [Dobrushin 1968; Foellmer 1982 / Kuensch 1982; Georgii,
      Gibbs Measures and Phase Transitions, Ch. 8]: if alpha < 1 then
      the Gibbs measure on any volume is unique and, for bounded local
      f, g at graph distance t,  |Cov(f,g)| <= osc(f) osc(g)
      sum_{n>=t} alpha^n = osc(f)osc(g) alpha^t/(1-alpha).
      EXTRACTION (proved here): T is compact (K_a Hilbert-Schmidt),
      self-adjoint, positive, with simple Perron vacuum psi_0 > 0
      (YM-8 floors). For bounded O, Cov(O(0),O(t)) =
      sum_{k>=2} |<psi_k, O psi_0>|^2 (lambda_k/lambda_1)^t  >= 
      |<psi_2, O psi_0>|^2 (lambda_2/lambda_1)^t. Since psi_0 > 0 a.e.,
      {O psi_0 : O bounded} is dense, so some bounded O has
      <psi_2, O psi_0> != 0; letting t -> infinity,
          lambda_2/lambda_1  <=  alpha(a,kappa)     FOR EVERY m.
      Hence the chain's reduced gap per time step satisfies
          Delta(S_m; a, kappa) >= -log alpha(a,kappa)   for all m >= 1,
      certified:  >= 0.263 at (6,1/16),  >= 0.145 at (8,1/8),
      >= 0.23 at (12,1/8)   [values pinned below with enclosures].
      The invariant (gauge) sector contains psi_0, so the bound holds
      there too. The free B-sector (YM-16 T1) is compatible: its ratio
      e^{-3a/4} is far below alpha.

WHAT THIS CLOSES AND WHAT IT DOES NOT.
  * CLOSED: volume-uniformity of the interacting gap on the
    bounded-overlap chain — the object YM-12 named, YM-14 opened,
    YM-16 certified to finite m, YM-17 localised — now certified for
    ALL m at strong coupling on the heat-kernel family. YM-11 gate 2
    (volume/IR) interacting half: closed on the chain (bounded degree),
    the physically relevant sharing pattern, at strong coupling.
  * NOT CLOSED: continuum. alpha < 1 needs a LARGE time step a (coarse
    time lattice) and SMALL kappa; along YM-9's trajectory kappa = theta
    a with theta = 1/64 the grid point (8,1/8) lies ON the trajectory,
    but the condition degrades as a -> 0 (delta_t -> infinity: the heat
    kernel sharpens) — so uniformity in the CUTOFF is exactly what
    Dobrushin cannot give. AF trajectory, tightness, OS reconstruction,
    non-triviality, metric universality, 2D lattice, Clay: OPEN.
  * The classical anchor is CITED, not rederived (CIRC-1). Everything
    else in the verdict is exact arithmetic.

Controls:
  C1  delta_s exact: e^{2 kappa} two-route enclosure (YM-5 pattern).
  C2  S_a tail certified; K_a lower bound 1 - S_a is positive on the
      grid (else delta_t undefined -> refuse).
  C3  fail-closed cells present and refused.
  C4  kappa = 0 consistency: the Dobrushin bound -log alpha(a,0) is
      strictly BELOW the exact free gap a C_{1/2} = 3a/4 (bound never
      exceeds truth).
  C5  B-sector consistency: e^{-3a/4} < alpha on every certified cell.
  C6  monotonicity: alpha decreasing in a at fixed kappa, increasing in
      kappa at fixed a (certified on the grid).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, log_iv, _dec, canonical_sha, LOG_TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402

_EXP_CACHE = {}

GRID = [(F(6), F(1, 16)), (F(8), F(1, 8)), (F(12), F(1, 8))]
FAIL_CELLS = [(F(6), F(1, 8)), (F(4), F(1, 16)), (F(10), F(1, 4))]
J_CUT = 4            # twice-spin cutoff for S_a, tail certified
C_HALF = F(3, 4)


def casimir(t: int) -> F:
    """C_j = j(j+1) with t = 2j."""
    return F(t, 2) * (F(t, 2) + 1)


def exp_neg(x: F) -> Iv:
    """enclosure of e^{-x}, x >= 0 rational: e^{-x} = (e^{-x/n})^n with
    x/n <= 1 so the Taylor bracket converges (YM-13 compound-route)."""
    if x in _EXP_CACHE:
        return _EXP_CACHE[x]
    n = max(1, int(x) + 1)
    e = exp_point(x / n)
    base = _r(Iv(F(1)) / Iv(e.lo, e.hi))
    out = Iv(F(1))
    for _ in range(n):
        out = _r(out * base)
    _EXP_CACHE[x] = out
    return out


def S_a(a: F) -> Iv:
    """S_a = sum_{t>=1} (t+1)^2 e^{-a C_t}, tail by geometric domination:
    for t >= J, ratio of consecutive terms <= ((t+2)/(t+1))^2 e^{-a(t+1)/... }
    we use the crude but certified bound term_{t+1}/term_t <= 4 e^{-a(t+3/2)/2}
    <= 4 e^{-a} for t >= 1 (C_{t+1}-C_t = (2t+3)/4 >= 5/4)."""
    acc = Iv(F(0))
    for t in range(1, J_CUT + 1):
        acc = acc + Iv(F((t + 1) ** 2)) * exp_neg(a * casimir(t))
    # tail: first omitted term times 1/(1-r), r <= 4 e^{-a*5/4}
    first = Iv(F((J_CUT + 2) ** 2)) * exp_neg(a * casimir(J_CUT + 1))
    r = (Iv(F(4)) * exp_neg(a * F(5, 4))).hi
    if r >= 1:
        raise ValueError("tail ratio not < 1; increase a or J_CUT")
    tail_hi = first.hi / (1 - r)
    return Iv(acc.lo, acc.hi + tail_hi)


def delta_t(a: F):
    s = S_a(a)
    if not (s.hi < 1):
        return None
    return (Iv(F(1)) + s) / (Iv(F(1)) - s)


def delta_s(kappa: F) -> Iv:
    e = exp_point(2 * kappa)
    return Iv(e.lo, e.hi)


def coeff(delta: Iv) -> Iv:
    return Iv(F(1)) - Iv(F(1)) / (delta * delta)


def alpha(a: F, kappa: F):
    dt = delta_t(a)
    if dt is None:
        return None
    return Iv(F(2)) * coeff(dt) + Iv(F(2)) * coeff(delta_s(kappa))


def cell(a: F, kappa: F):
    al = alpha(a, kappa)
    if al is None:
        return {"status": "REFUSED_KA_LOWER_BOUND_NOT_POSITIVE"}
    row = {"alpha_lo": _dec(al.lo, 20), "alpha_hi": _dec(al.hi, 20)}
    if al.hi < 1:
        g = -log_iv(Iv(al.hi), LOG_TERMS)
        row["status"] = "CERTIFIED_UNIFORM_IN_m"
        row["gap_per_time_step_lower_bound"] = _dec(g.lo, 20)
        row["gap_per_unit_time_lower_bound"] = _dec(g.lo / a, 20)
    else:
        row["status"] = "REFUSED_ALPHA_NOT_BELOW_1"
    return row


def run():
    rows = {}
    c2 = c3 = c4 = c5 = c6 = True
    all_cert = True
    for a, k in GRID:
        r = cell(a, k)
        rows[f"a={a},kappa={k}"] = r
        if r["status"] != "CERTIFIED_UNIFORM_IN_m":
            all_cert = False
        # C5 B-sector: e^{-3a/4} < alpha
        al = alpha(a, k)
        if not (exp_neg(C_HALF * a).hi < al.lo):
            c5 = False
        # C4 kappa = 0 consistency
        al0 = alpha(a, F(0))
        g0 = -log_iv(Iv(al0.hi), LOG_TERMS)
        if not (g0.hi < C_HALF * a):
            c4 = False
    for a, k in FAIL_CELLS:
        r = cell(a, k)
        rows[f"a={a},kappa={k}"] = r
        if r["status"] == "CERTIFIED_UNIFORM_IN_m":
            c3 = False
    # C1 delta_s two-route
    e1, e2 = exp_point(F(1, 8)), exp_point(F(1, 4))  # e^{2k}, k=1/16 vs 1/8
    c1 = e1.lo * e1.lo <= e2.hi and e2.lo <= e1.hi * e1.hi
    # C2 K_a lower bound positive on grid
    c2 = all(S_a(a).hi < 1 for a, _ in GRID)
    # C6 monotonicity
    a6 = alpha(F(6), F(1, 16))
    a8 = alpha(F(8), F(1, 16))
    k1 = alpha(F(8), F(1, 16))
    k2 = alpha(F(8), F(1, 8))
    c6 = (a8.hi < a6.lo) and (k1.hi < k2.lo)

    ok = c1 and c2 and c3 and c4 and c5 and c6 and all_cert
    cert = {
        "certificate_type": "YM19_DOBRUSHIN_DOCK_VOLUME_UNIFORM_GAP",
        "claim_status": "first_certified_m_uniform_interacting_gap_on_the_"
                        "bounded_overlap_chain__strong_coupling_fixed_time_step",
        "standing_correction": {
            "target": "YM-17 T3",
            "content": "vacuum-tracking inequality is not sufficient for the "
                       "gap (interleaved free energy != sum of halves, "
                       "exponential mismatch generic); replaced by the "
                       "correlation-decay route; YM-18 measurement retained "
                       "as calibration only",
        },
        "theorems": {
            "T1_dobrushin_coefficient_bound":
                "C_ij <= 1 - 1/delta^2, delta = sup w / inf w (elementary "
                "proof in docstring); delta_s = e^{2kappa} exact, delta_t <= "
                "(1+S_a)/(1-S_a) with certified tail",
            "T2_dobrushin_condition_uniform_in_m":
                "alpha = 2(1-1/delta_t^2) + 2(1-1/delta_s^2) bounds the "
                "row sum on the strip for every m; alpha < 1 certified on "
                "the grid, refused on the fail cells",
            "T3_uniform_gap_via_pinned_anchor":
                "Dobrushin uniqueness => |Cov| <= osc osc alpha^t/(1-alpha) "
                "(cited); extraction lemma (proved): lambda_2/lambda_1 <= "
                "alpha for every m; gap >= -log alpha uniformly in volume",
        },
        "anchor": {
            "circ_status": "CIRC-1: cited, not rederived (as YM-8 Jentzsch)",
            "references": ["Dobrushin 1968 (uniqueness)",
                           "Foellmer 1982 / Kuensch 1982 (covariance decay)",
                           "Georgii, Gibbs Measures and Phase Transitions, "
                           "ch. 8"],
        },
        "grid": rows,
        "closes": "volume-uniformity of the interacting gap on the "
                  "bounded-overlap chain (YM-12 named object; YM-11 gate 2 "
                  "interacting half, bounded degree) — at strong coupling, "
                  "fixed time step",
        "honest_remainder": {
            "cutoff": "alpha -> infinity as a -> 0 (delta_t blows up): "
                      "NOT uniform in the cutoff; continuum OPEN",
            "coupling": "small kappa only; weak-coupling/AF regime OPEN",
            "scope": "toy chain (2 holonomies per block); 2D lattice, "
                     "tightness, OS, non-triviality, metric universality, "
                     "Clay OPEN",
        },
        "controls": {
            "C1_delta_s_two_route": bool(c1),
            "C2_Ka_lower_bound_positive": bool(c2),
            "C3_fail_closed_cells_refused": bool(c3),
            "C4_kappa0_bound_below_exact_free_gap": bool(c4),
            "C5_B_sector_consistent": bool(c5),
            "C6_monotone_in_a_and_kappa": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM19_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM19.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" ", k, v["status"], "alpha<=", v.get("alpha_hi", "")[:8],
              "gap/step>=", v.get("gap_per_time_step_lower_bound", "")[:8])
    print("sha256:", sha)
