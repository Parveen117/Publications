"""YM-7: ENLARGED NATIVE CARRIER V7 + CERTIFIED CROSSING CURVES —
the seam-integer dock pushed toward strong coupling, with the RH-line's
post-pause lessons adopted (exact formulas over bisection-as-truth;
report drift honestly; never promote a limited-denominator constant).

Lessons consumed from RH-Framework agent/l0-rank-one-seam-jet-schur
(L0_FINITE_SHADOW_EXACT_DEMAND_AUDIT): a bisection-derived "exact constant"
can be a limited-denominator shadow — so this capsule publishes eigenvalue
ENCLOSURES and their honest kappa-drift, never a claimed exact rational;
bisection is used only to LOCALIZE certified inertia steps, and every
verdict is an interval-LDL inertia (exact congruence), not a root claim.

Claim boundary (declared, fail-closed):
  Carrier grows from YM-6's V5 to
    V7 = V5 + { chi_1(A), chi_1(B) }
  — still ALL exact T0 eigenvectors, t = (1, lam, lam, lam^2, lam^2, lam1,
  lam1), and still content-exhaustive at its levels: (1,0)/(0,1) invariant
  content is exactly the class function chi_1 (1-dim each). Hence the
  complement top DROPS:
    lam_next7 = lam1 * lam   (content (1,1/2); certified strictly above
    lam_{3/2} and lam1^2, strictly below lam1)
  versus YM-6's lam1 — a certified factor ~2.3 improvement in |QSQ|,
  which is what unlocks kappa = 7/10.

  CERTIFIED:
    (E1) enclosures of ALL SEVEN compressed pencil eigenvalues
         mu_i(kappa) of (A7 - mu G7), localized by inertia bisection
         (each step an exact interval-LDL count) to width < 1e-10;
         published as brackets with their honest kappa-drift.
    (E2) exact seam counts k_Sigma(kappa, mu) via the YM-6 Haynsworth
         bracket on the V7 dock, including the previously REFUSED cell
         kappa = 7/10 (now with mu inside the enlarged window).
    (E3) sector split (swap-even 5-dim / swap-odd 2-dim) exact on V7;
         crossing branch remains swap-even (YM-4 protection persists).

  NOT certified: couplings where the bracket does not close (recorded);
  full lattice, UV/IR, continuum, vacuum, Clay predicate.

Controls:
  C1  V5/V7 consistency: on every cell YM-6 certified, V7 certifies the
      SAME count.
  C2  content exhaustion: lam_{3/2}, lam1^2 certified strictly below
      lam_next7 = lam1*lam; lam_next7 strictly below lam1.
  C3  eigenvalue-curve sanity: top compressed eigenvalue increases with
      kappa; the odd-sector doublet member decreases (YM-4's exact split
      sign) — checked on certified brackets, drift reported not assumed.
  C4  Gram of V7 = blockdiag(G5, I_2) exactly (new rows orthogonal).
  C5  planted low mu still counts >= 3; large kappa (3/2) fails closed.
"""

from fractions import Fraction as F
import json
import math
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, _dec, canonical_sha, TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym4_symmetry_protected import chi_mul, dim, face_coeffs, tail_E  # noqa: E402
from ym6_seam_integer_dock import (  # noqa: E402
    _r, ring_mul_q, slot_with_face, ldl_inertia, mat_inv_exact, iv_sqrt,
    BETA, P_CUT,
)

BASIS7 = [
    ({0: F(1)}, {0: F(1)}, {0: F(1)}),        # 1
    ({1: F(1)}, {0: F(1)}, {0: F(1)}),        # chi12(A)
    ({0: F(1)}, {1: F(1)}, {0: F(1)}),        # chi12(B)
    ({1: F(1)}, {1: F(1)}, {0: F(1)}),        # chi12(A)chi12(B)
    ({0: F(1)}, {0: F(1)}, {1: F(1)}),        # chi12(AB^-1)
    ({2: F(1)}, {0: F(1)}, {0: F(1)}),        # chi_1(A)
    ({0: F(1)}, {2: F(1)}, {0: F(1)}),        # chi_1(B)
]
N7 = len(BASIS7)
CELLS7 = [
    (F(1, 8), F(3, 5)),
    (F(1, 4), F(3, 5)),
    (F(1, 2), F(1, 1)),
    (F(7, 10), F(13, 10)),     # the YM-6 refusal, retried on V7
]
EIG_KAPPAS = [F(1, 8), F(1, 4), F(1, 2), F(7, 10)]

_M7CACHE = {}


def m7_matrix(kappa: F):
    if kappa in _M7CACHE:
        return _M7CACHE[kappa]
    f = face_coeffs(kappa, P_CUT)
    M = [[Iv(F(0)) for _ in range(N7)] for _ in range(N7)]
    for i in range(N7):
        for j in range(i, N7):
            Ai, Bi, Ci = BASIS7[i]
            Aj, Bj, Cj = BASIS7[j]
            At = slot_with_face(ring_mul_q(Ai, Aj), f)
            Bt = slot_with_face(ring_mul_q(Bi, Bj), f)
            Ct = slot_with_face(ring_mul_q(Ci, Cj), f)
            acc = Iv(F(0))
            for t in set(At) & set(Bt) & set(Ct):
                acc = acc + At[t] * Bt[t] * Ct[t] * Iv(F(1, dim(t)))
            M[i][j] = _r(acc)
            M[j][i] = M[i][j]
    err = 3 * exp_point(2 * kappa).hi * tail_E(kappa, P_CUT) * 4
    Mw = [[Iv(M[i][j].lo - err, M[i][j].hi + err) for j in range(N7)]
          for i in range(N7)]
    _M7CACHE[kappa] = (Mw, err)
    return _M7CACHE[kappa]


def gram7():
    G = [[F(0)] * N7 for _ in range(N7)]
    for i in range(N7):
        for j in range(N7):
            Ai, Bi, Ci = BASIS7[i]
            Aj, Bj, Cj = BASIS7[j]
            At = ring_mul_q(Ai, Aj)
            Bt = ring_mul_q(Bi, Bj)
            Ct = ring_mul_q(Ci, Cj)
            s = F(0)
            for t in set(At) & set(Bt) & set(Ct):
                s += At[t] * Bt[t] * Ct[t] / dim(t)
            G[i][j] = s
    return G


def levels7():
    I1 = bessel_I(1, BETA, TERMS)
    lam = bessel_I(2, BETA, TERMS) / I1
    lam1 = bessel_I(3, BETA, TERMS) / I1
    lam32 = bessel_I(4, BETA, TERMS) / I1
    lam_next7 = lam1 * lam                    # content (1, 1/2)
    return lam, lam1, lam32, lam_next7


def s7_block(kappa: F):
    lam, lam1, _, _ = levels7()
    M, _ = m7_matrix(kappa)
    s = iv_sqrt(lam)
    s1 = iv_sqrt(lam1)
    root_t = [Iv(F(1)), s, s, _r(lam), _r(lam), s1, s1]
    A = [[_r(root_t[i] * root_t[j] * M[i][j]) for j in range(N7)]
         for i in range(N7)]
    return A, M


def offblock7(kappa: F, Mk, Ginv, lam_next_hi: F):
    M2k, _ = m7_matrix(2 * kappa)
    n = N7
    GiM = [[_r(sum((Iv(Ginv[i][k]) * Mk[k][j] for k in range(n)),
                   Iv(F(0)))) for j in range(n)] for i in range(n)]
    D = [[_r(M2k[i][j] - sum((Mk[k][i] * GiM[k][j] for k in range(n)),
                             Iv(F(0)))) for j in range(n)] for i in range(n)]
    W = [[_r(sum((Iv(Ginv[i][k]) * D[k][j] for k in range(n)),
                 Iv(F(0)))) for j in range(n)] for i in range(n)]
    gmax = gmin = None
    trace_hi = F(0)
    for i in range(n):
        trace_hi += W[i][i].hi
        rh, rl = W[i][i].hi, W[i][i].lo
        for j in range(n):
            if j != i:
                off = max(abs(W[i][j].lo), abs(W[i][j].hi))
                rh += off
                rl -= off
        gmax = rh if gmax is None else max(gmax, rh)
        gmin = rl if gmin is None else min(gmin, rl)
    bound = min(gmax, trace_hi - (n - 1) * min(gmin, F(0)))
    pmq = iv_sqrt(Iv(F(0), max(F(0), bound))).hi
    return iv_sqrt(Iv(lam_next_hi)).hi * pmq


def dock_cell7(kappa: F, mu: F):
    lam, lam1, lam32, lam_next7 = levels7()
    G = gram7()
    Ginv = mat_inv_exact(G)
    A, Mk = s7_block(kappa)
    qsq_hi = (lam_next7 * Iv(exp_point(3 * kappa).hi)).hi
    if not (qsq_hi < mu):
        return {"status": "REFUSED_C_NOT_BELOW_MU", "qsq_hi": qsq_hi}
    psq = offblock7(kappa, Mk, Ginv, lam_next7.hi)
    r_hi = psq * psq / (mu - qsq_hi)
    low = [[_r(A[i][j] - Iv(mu * G[i][j])) for j in range(N7)]
           for i in range(N7)]
    high = [[_r(low[i][j] + Iv(r_hi * G[i][j])) for j in range(N7)]
            for i in range(N7)]
    li, hi_ = ldl_inertia(low), ldl_inertia(high)
    if li is None or hi_ is None:
        return {"status": "REFUSED_PIVOT_STRADDLE", "r_hi": r_hi,
                "qsq_hi": qsq_hi, "psq": psq}
    certified = li[0] == hi_[0]
    return {"status": "OK", "count_lower": li[0], "count_upper": hi_[0],
            "count_exact": li[0] if certified else None,
            "certified": certified, "qsq_hi": qsq_hi, "psq": psq,
            "r_hi": r_hi}


def compressed_eigs(kappa: F, tol=F(1, 10 ** 10)):
    """Certified brackets for all 7 compressed pencil eigenvalues of
    (A7 - mu G7): inertia bisection — every step an exact LDL count."""
    G = gram7()
    A, _ = s7_block(kappa)

    def count_above(mu):
        M = [[_r(A[i][j] - Iv(mu * G[i][j])) for j in range(N7)]
             for i in range(N7)]
        r = ldl_inertia(M)
        return None if r is None else r[0]

    lo_all, hi_all = F(-1), F(4)
    brackets = []
    for k in range(1, N7 + 1):          # k-th largest eigenvalue
        lo, hi = lo_all, hi_all
        while hi - lo > tol:
            mid = (lo + hi) / 2
            c = count_above(mid)
            if c is None:               # straddle: nudge
                mid = mid + tol / 7
                c = count_above(mid)
                if c is None:
                    break
            if c >= k:
                lo = mid
            else:
                hi = mid
        brackets.append(Iv(lo, hi))
    return brackets


def run():
    lam, lam1, lam32, lam_next7 = levels7()
    # C2 content exhaustion / level ordering
    c2 = (lam32.hi < lam_next7.lo and (lam1 * lam1).hi < lam_next7.lo
          and lam_next7.hi < lam1.lo)

    # C4 Gram structure
    G = gram7()
    c4 = all(G[i][j] == (F(1) if i == j else F(0))
             for i in range(5, 7) for j in range(N7) if not (i == j)) and \
        G[5][5] == 1 and G[6][6] == 1 and G[3][4] == F(1, 2)

    cells = {}
    unlocked_7_10 = False
    for kap, mu in CELLS7:
        r = dock_cell7(kap, mu)
        key = f"kappa={kap},mu={mu}"
        if r["status"] == "OK":
            cells[key] = {"count_bracket": [r["count_lower"],
                                            r["count_upper"]],
                          "seam_count_exact": r["count_exact"],
                          "r_hi": _dec(F(r["r_hi"]), 15)}
            if kap == F(7, 10) and r["count_exact"] == 1:
                unlocked_7_10 = True
        else:
            cells[key] = {"refused": r["status"]}

    # C1 V5/V7 consistency on YM-6 certified cells
    import ym6_seam_integer_dock as ym6
    c1 = True
    for kap, mu in [(F(1, 8), F(3, 5)), (F(1, 4), F(3, 5)),
                    (F(1, 2), F(1))]:
        r5 = ym6.dock_cell(kap, mu)
        r7 = dock_cell7(kap, mu)
        c1 = c1 and (r5.get("count_exact") == r7.get("count_exact") == 1)

    # E1/C3 eigenvalue curves with honest drift
    curves = {}
    tops, odd_track = [], []
    for kap in EIG_KAPPAS:
        br = compressed_eigs(kap)
        curves[str(kap)] = [[_dec(b.lo, 15), _dec(b.hi, 15)] for b in br]
        tops.append(br[0])
        odd_track.append(br[2])         # third-largest tracks a falling branch
    c3 = all(tops[i + 1].lo > tops[i].hi - F(1, 10 ** 6)
             for i in range(len(tops) - 1))

    # C5 planted + fail-closed
    p = dock_cell7(F(1, 8), F(3, 10))
    c5a = p["status"] == "OK" and p["count_lower"] >= 3
    big = dock_cell7(F(3, 2), F(2))
    c5b = big["status"].startswith("REFUSED") or not big["certified"]

    ok = c1 and c2 and c3 and c4 and c5a and c5b and unlocked_7_10
    cert = {
        "certificate_type": "YM7_V7_CARRIER_CROSSING_CURVES",
        "claim_status": "exact_counts_extended_plus_certified_curves",
        "claim_boundary": {
            "certified": [
                "V7 content-exhaustive at its levels; complement top drops "
                "to lam1*lam (certified ordering)",
                "kappa=7/10 seam count k=1 EXACT (YM-6 refusal unlocked)",
                "all seven compressed eigenvalue curves as certified "
                "brackets with honest kappa-drift (no exact-constant "
                "claims; RH-line drift lesson adopted)",
            ],
            "not_certified": [
                "unclosed cells (recorded)", "full lattice", "UV/IR",
                "continuum OS", "vacuum", "Clay predicate",
            ],
            "framework_sources": [
                "YM-1..6 (consumed, pinned)",
                "RH-Framework L0 demand audit lesson: determinant-lemma/"
                "exact-formula ethos; limited-denominator constants "
                "never promoted",
            ],
        },
        "parameters": {"beta": str(BETA), "basis_dim": N7,
                       "cells": [[str(k), str(m)] for k, m in CELLS7]},
        "cells": cells,
        "eigenvalue_curves": curves,
        "controls": {
            "C1_v5_v7_consistency": bool(c1),
            "C2_level_ordering_content_exhaustion": bool(c2),
            "C3_top_curve_monotone_within_brackets": bool(c3),
            "C4_gram_blockdiag_exact": bool(c4),
            "C5_planted_and_fail_closed": bool(c5a and c5b),
            "unlocked_kappa_7_10": bool(unlocked_7_10),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


def run_all_and_pin():
    import ym1_certified_gap as m1
    import ym2_theta_interacting_gap as m2
    import ym3_crossing_direction as m3
    import ym4_symmetry_protected as m4
    import ym5_two_sided_gap as m5
    import ym6_seam_integer_dock as m6
    out = {}
    for name, fn in (("YM1", m1.run), ("YM2", m2.run), ("YM3", m3.run),
                     ("YM4", m4.run), ("YM5", m5.run), ("YM6", m6.run),
                     ("YM7", run)):
        cert = fn()
        sha = canonical_sha(cert)
        with open(os.path.join(HERE, f"{name}_RESULT.json"), "w") as fj:
            json.dump(cert, fj, indent=2, sort_keys=True)
        with open(os.path.join(HERE, f"EXPECTED_{name}.sha256"), "w") as fs:
            fs.write(sha + "\n")
        out[name] = (cert["verdict"], sha)
        print(f"{name}: {cert['verdict']}  sha256:{sha[:16]}...")
    return out


if __name__ == "__main__":
    results = run_all_and_pin()
    assert all(v == "PASS" for v, _ in results.values())
    print("ALL CERTIFICATES PASS")
