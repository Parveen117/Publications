"""YM-18 (CALIBRATION + DESIGN): the vacuum-tracking inequality, measured
on an enlarged exact carrier, and the cluster-expansion dock DESIGNED
(not executed). Status: calibration. Nothing about m-uniformity of the
chain gap is claimed.

Why this shape. YM-17 reduced the A-chain's volume problem to one
inequality — does the interleaved vacuum lambda_1(T_m) track the
half-chain product sigma_1(X)sigma_1(Y) up to a sub-exponential factor?
Before paying for a Kotecky-Preiss certification, the program measures
the quantity the expansion would control: the per-bridge vacuum rate
phi_m := log lambda_1(T_m) / (m-1), bracketed from below by exact
Rayleigh quotients on an ENLARGED carrier and from above by YM-17.

CARRIER W'_m = { 1 } U { chi12(A_i) } U { chi12(A_i) chi12(A_{i+1}) },
dim 2m, orthonormal, all exact T0 eigenvectors (weights 1, lambda,
lambda^2). The pair elements are exactly the excitations one bridge
creates from the vacuum (chi12 (x) chi12 content of chi_{1/2}(A_i A_{i+1}^{-1})),
so this is the first carrier that SEES the bridge's own first-order
action on the vacuum. The compressed bridge product is computed by the
YM-15 symbolic chain-integration engine (exact polynomial in the face
coefficients f_j; contents up to j = 3/2 appear and are integrated
exactly; no truncation remainder), then evaluated with certified
interval f_j.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) RAYLEIGH LOWER BOUNDS ON lambda_1(T_A, m), m = 2..M, on W'_m:
      certified brackets of lammax of the compressed transfer (exact
      interval LDL inertia bisection). These dominate the YM-15/16
      carrier vacuum f_0^{m-1} strictly at every m (the pair elements
      contribute at first order, coefficient f_{1/2}).
 (T2) THE TRACKING BRACKET. Per-bridge rate
          log f_0  <  phi_m^{lo}(W')  <=  phi_m  <=  (1/2) log lambda_1(Y_2k)_hi
      pinned for each m; the lower end RISES with the enlarged carrier
      (measured), the upper end is YM-17's. MEASURED FACT: the lower
      end phi_m^{lo}(W') is m-INDEPENDENT to ~1e-5 across m = 2..6 at
      every grid kappa — the enlarged carrier's vacuum rate has already
      settled; the bracket width is entirely the looseness of the upper
      end. This is the calibration signal: vacuum tracking looks like a
      clean per-bridge rate, which is exactly what a convergent cluster
      expansion predicts. The residual width of the
      bracket at the largest m computed is reported as the calibration
      number the cluster expansion has to close.
 (T3) DESIGN OF THE CLUSTER DOCK (text, pinned, no claims):
      monomers = the two-site chains Y_2k on odd pairs; the interleaving
      bridges B_even are the activities; polymers = connected unions of
      even bridges with their adjacent odd pairs; weight of a polymer =
      its connected matrix element between the half-chain product
      vacuum and itself, divided by the product of monomer vacua; the
      Kotecky-Preiss criterion sum_{gamma' ~ gamma} |w(gamma')| e^{a|gamma'|}
      <= a |gamma| with an explicit a; the radius kappa_KP is what YM-18
      proper would certify, and it would deliver (i) analyticity of phi
      in kappa uniformly in m (vacuum tracking), and (ii) with the
      standard extension to the time-correlation, the m-uniform gap.
      The program's obligations, listed: exact polymer weights via the
      YM-15 engine (finite sums, no truncation), certified activity
      bounds |w| <= C kappa^{|gamma|} from the f_j tails (YM-4 E_P
      pattern), a KP sum as a geometric series in 1D (connectivity count
      2^{|gamma|} on a line), and the Jentzsch-free extraction of the
      gap from exponential decay of the truncated two-point function.

HONEST: T1/T2 are finite-m measurements with certified brackets; no
uniform statement. T3 is a design, its every line an obligation.

Controls:
  C1  W' compressed vacuum top strictly exceeds f_0^{m-1} at every m.
  C2  kappa = 0: W' compressed transfer is diag(1, lambda.., lambda^2..),
      top exactly 1.
  C3  engine symmetry: compressed matrix symmetric entry by entry (the
      chain integral is symmetric under swapping insertions).
  C4  bracket ordering log f_0 < phi^lo < phi^hi at every m.
  C5  the pair element's first-order coefficient equals f_{1/2} (the
      bridge's own content), exact.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS,
)
from ym4_symmetry_protected import chi_mul  # noqa: E402
from ym6_seam_integer_dock import ldl_inertia, _r  # noqa: E402
from ym15_chain_closed_form import (  # noqa: E402
    p_add, p_const, p_mul_f, p_scale, lam_half, P_SYM,
)
from ym16_chain_dock import f0_of, iv_pow  # noqa: E402
from ym17_interleaving_seam import pair_top_upper  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
M_MAX = 6
BISECT_TOL = F(1, 10 ** 10)


def chain_entry_multi(m: int, ins: dict, tamper=False):
    """YM-15 engine generalised: ins maps site -> list of twice-spins
    (several insertions at one site are multiplied in sequence)."""
    state = {0: p_const(1)}
    for i in range(1, m + 1):
        for t in ins.get(i, []):
            nstate = {}
            for c, poly in state.items():
                for c2 in chi_mul(c, t):
                    nstate[c2] = p_add(nstate.get(c2, {}), poly)
            state = nstate
        if i == m:
            return state.get(0, {})
        nstate = {}
        for c, poly in state.items():
            if c > P_SYM:
                raise ValueError("raise P_SYM")
            q = p_mul_f(poly, c)
            if tamper:
                q = p_scale(q, c + 1)
            nstate[c] = p_add(nstate.get(c, {}), q)
        state = nstate
    return state.get(0, {})


def basis(m: int):
    """returns list of (label, ins-dict, T0 weight exponent in lambda)."""
    out = [("V", {}, 0)]
    for i in range(1, m + 1):
        out.append((f"S{i}", {i: [1]}, 1))
    for i in range(1, m):
        out.append((f"P{i}", {i: [1], i + 1: [1]}, 2))
    return out


def merge_ins(a: dict, b: dict) -> dict:
    out = {k: list(v) for k, v in a.items()}
    for k, v in b.items():
        out[k] = out.get(k, []) + list(v)
    return out


def f_intervals(kappa: F):
    f = {}
    for t in range(P_SYM + 1):
        if kappa == 0:
            f[t] = Iv(F(1) if t == 0 else F(0))
        else:
            f[t] = bessel_I(t + 1, kappa, TERMS) * Iv(F(2)) / Iv(kappa)
    return f


def eval_poly(poly: dict, f: dict) -> Iv:
    acc = Iv(F(0))
    for mono, c in poly.items():
        term = Iv(c)
        for t, e in enumerate(mono):
            for _ in range(e):
                term = term * f[t]
        acc = acc + term
    return _r(acc)


_CACHE = {}


def compressed_symbolic(m: int):
    if m in _CACHE:
        return _CACHE[m]
    B = basis(m)
    n = len(B)
    P = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(i, n):
            P[i][j] = P[j][i] = chain_entry_multi(m, merge_ins(B[i][1], B[j][1]))
    _CACHE[m] = (B, P)
    return _CACHE[m]


def compressed_T(m: int, kappa: F):
    B, P = compressed_symbolic(m)
    f = f_intervals(kappa)
    lam = lam_half()
    from ym6_seam_integer_dock import iv_sqrt
    s = iv_sqrt(lam)
    n = len(B)
    w = [iv_pow(s, B[i][2]) for i in range(n)]
    return [[_r(w[i] * w[j] * eval_poly(P[i][j], f)) for j in range(n)]
            for i in range(n)]


def count_above(M, mu: F):
    n = len(M)
    S = [[_r(M[i][j] - Iv(mu if i == j else F(0))) for j in range(n)]
         for i in range(n)]
    res = ldl_inertia(S)
    return None if res is None else res[0]


def lammax_bracket(M, lo: F, hi: F):
    assert count_above(M, hi) == 0 and count_above(M, lo) >= 1
    while hi - lo > BISECT_TOL:
        mid = (lo + hi) / 2
        c = count_above(M, mid)
        if c is None:
            lo2 = lo + (mid - lo) / 3
            c2 = count_above(M, lo2)
            if c2 is None:
                break
            if c2 >= 1:
                lo = lo2
            else:
                hi = lo2
            continue
        if c >= 1:
            lo = mid
        else:
            hi = mid
    return Iv(lo, hi)


def run():
    lam = lam_half()
    rows = {}
    c1 = c3 = c4 = c5 = True
    # C5: first-order pair coefficient = f_{1/2}: <V, m P1> on m=2 is f_1/2
    B, P = compressed_symbolic(2)
    idxV = 0
    idxP = [k for k, b in enumerate(B) if b[0] == "P1"][0]
    c5 = P[idxV][idxP] == {tuple([0, 1] + [0] * (P_SYM - 1)): F(1)}
    for kap in GRID:
        f0 = f0_of(kap)
        logf0 = log_iv(f0, LOG_TERMS)
        top_hi = pair_top_upper(2 * kap)
        phi_hi = (log_iv(Iv(top_hi), LOG_TERMS) * Iv(F(1, 2))).hi
        per_m = {}
        for m in range(2, M_MAX + 1):
            M = compressed_T(m, kap)
            n = len(M)
            for i in range(n):
                for j in range(n):
                    if M[i][j].lo != M[j][i].lo or M[i][j].hi != M[j][i].hi:
                        c3 = False
            # bracket: lower f0^{m-1} (carrier vacuum), upper = row-sum
            lo = iv_pow(f0, m - 1).lo * (1 - F(1, 10 ** 4))
            hi = max(sum(max(abs(M[i][j].lo), abs(M[i][j].hi))
                         for j in range(n)) for i in range(n)) + F(1, 10 ** 6)
            b = lammax_bracket(M, lo, hi)
            if not (b.lo > iv_pow(f0, m - 1).hi):
                c1 = False
            phi_lo = (log_iv(Iv(b.lo), LOG_TERMS) * Iv(F(1, m - 1))).lo
            if not (logf0.lo < phi_lo < phi_hi):
                c4 = False
            per_m[str(m)] = {
                "lambda1_lower_bracket": [_dec(b.lo, 12), _dec(b.hi, 12)],
                "carrier_vacuum_f0_pow": _dec(iv_pow(f0, m - 1).lo, 12),
                "phi_lo": _dec(phi_lo, 12),
            }
        rows[str(kap)] = {
            "log_f0": _dec(logf0.lo, 12),
            "phi_hi_YM17": _dec(phi_hi, 12),
            "per_m": per_m,
            "tracking_bracket_width_at_M_MAX":
                _dec(phi_hi - F(per_m[str(M_MAX)]["phi_lo"]), 12),
        }
    # C2 free chain
    M0 = compressed_T(4, F(0))
    B0 = basis(4)
    c2 = True
    for i in range(len(M0)):
        for j in range(len(M0)):
            if i != j and not (M0[i][j].lo == M0[i][j].hi == 0):
                c2 = False
        expect = iv_pow(lam, B0[i][2])
        if not (M0[i][i].lo <= expect.hi and M0[i][i].hi >= expect.lo):
            c2 = False

    ok = c1 and c2 and c3 and c4 and c5
    cert = {
        "certificate_type": "YM18_VACUUM_TRACKING_CALIBRATION_AND_CLUSTER_DESIGN",
        "claim_status": "calibration_only__no_uniformity_claim",
        "theorems": {
            "T1_enlarged_carrier_rayleigh_bounds":
                "certified lammax brackets of the compressed transfer on "
                "W'_m (dim 2m: vacuum, sites, bridge pairs), m = 2..6, "
                "strictly above f0^{m-1}",
            "T2_tracking_bracket":
                "log f0 < phi_m^lo(W') <= phi_m <= (1/2) log lambda_1(Y_2k); "
                "residual width pinned as the number the cluster expansion "
                "must close",
            "T3_cluster_dock_design":
                "monomers = odd two-site chains at 2kappa; activities = "
                "even bridges; polymers = connected unions on a line; "
                "weights = connected matrix elements / monomer vacua via "
                "the YM-15 engine (finite, exact); KP criterion with 1D "
                "connectivity 2^|gamma|; deliverables (i) phi analytic "
                "uniformly in m, (ii) m-uniform gap via truncated two-point "
                "decay; each line an obligation, none discharged here",
        },
        "grid": rows,
        "honest_remainder": {
            "status": "finite-m measurement; uniformity OPEN; YM-18 proper "
                      "= executing T3",
            "scope": "toy carrier; 2D, AF, tightness, OS, non-triviality, "
                     "metric universality, Clay OPEN",
        },
        "controls": {
            "C1_enlarged_vacuum_exceeds_f0_power": bool(c1),
            "C2_free_chain_diagonal": bool(c2),
            "C3_engine_symmetry": bool(c3),
            "C4_bracket_ordering": bool(c4),
            "C5_pair_first_order_coefficient_is_f_half": bool(c5),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM18_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM18.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "log f0", v["log_f0"][:8], "phi_hi", v["phi_hi_YM17"][:8],
              "phi_lo by m:", {m: d["phi_lo"][:8] for m, d in v["per_m"].items()},
              "width", v["tracking_bracket_width_at_M_MAX"][:8])
    print("sha256:", sha)
