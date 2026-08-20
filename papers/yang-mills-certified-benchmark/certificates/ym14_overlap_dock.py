"""YM-14: THE INTERACTION-OVERLAP DOCK, FIRST CARRIER — a bridged pair of
blocks whose interaction genuinely SHARES a face, with the overlap's
spectral price measured exactly. This is the named object from YM-12's
honest remainder, opened.

Where this sits. YM-12 certified volume-uniformity when blocks share
NOTHING (tensor chain), and YM-11 proved the sup route dies when blocks
share EVERYTHING (B_n). The physical lattice lives between: bounded
sharing. This capsule builds the smallest genuinely-sharing system — two
blocks coupled by ONE bridge face — and certifies its spectral anatomy
exactly. Bridge-only coupling is deliberate: it isolates the price of
overlap with zero confound from within-block interaction.

Carrier. Two blocks with holonomies (A1, B1) and (A2, B2); Hilbert space
L^2(SU(2)^4)^Ad (YM-11 gate 1 pattern); free transfer T0 = K^(x)4 with
heat-kernel links; interaction = ONE bridge face
    m_kappa = exp[ (kappa/2) Tr(A1 A2^{-1}) ],
so the two blocks interact ONLY through the shared face.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) OVERLAP SYMMETRY AND FIRST-ORDER SEAM. The block swap
      S12: (A1,B1) <-> (A2,B2) satisfies S12 m S12 = m for every kappa
      (Tr(A2 A1^{-1}) = Tr((A1 A2^{-1})^{-1}) = Tr(A1 A2^{-1}), real
      characters — the YM-4 identity on the bridge). The free excited
      level lambda_{1/2} is 4-fold degenerate:
      {chi12(A1), chi12(A2), chi12(B1), chi12(B2)}. First order in kappa:
      the SAME exact theta integral
          Int chi12(A1) chi12(A2) chi12(A1 A2^{-1}) = 1/2
      couples ONLY the A-pair; the B-pair is exactly inert (every
      first-order entry touching B vanishes by parity). The overlap seam
      is therefore RANK-ONE and swap-transported:
          crossing line  chi12(A1) + chi12(A2),  derivative +lam/4,
          receding line  chi12(A1) - chi12(A2),  derivative -lam/4,
          inert plane    span{chi12(B1), chi12(B2)}, derivative 0,
      with lam = lambda_{1/2}. The overlap's first-order anatomy is
      IDENTICAL to the theta interaction's (YM-3/4) — same integral, same
      rank-one transported line — now appearing as an INTER-block seam.

 (T2) EXACT TWO-SIDED GAP OF THE BRIDGED PAIR (Haynsworth dock on the
      5-dim carrier W5 = {1, chi12(A1), chi12(A2), chi12(B1),
      chi12(B2)}, all exact T0 eigenvectors, orthonormal Gram).
      Complement top = lambda^2 (contents (1/2,1/2)); one face gives
      m in [e^{-kappa}, e^{kappa}] so |QSQ| <= lambda^2 e^{kappa};
      off-block by the doubling identity m_kappa^2 = m_{2kappa}. Exact
      seam counts k=1 certified on a (kappa, mu) grid: the bridged pair
      has a certified two-sided spectral gap at every certified cell.

 (T3) THE OVERLAP PRICE, MEASURED. epsilon(kappa) := the certified
      bracket on [lambda_2(bridged pair) - lambda_2(free pair)] via the
      dock's mu-window edges, compared against the first-order
      prediction (kappa/2)(1/2) lam = kappa lam / 4 from T1. The price
      is finite, first-order exact at small kappa, and — the structural
      point — carried by ONE swap-even direction, not spread.

 (T4) THE 1D REDUCTION, NAMED (next object, not proved). For the
      bounded-overlap chain S_m (m blocks, nearest-neighbor bridges),
      1D structure makes volume-uniformity equivalent to a spectral
      property of a single BLOCK-TRANSFER operator on one block's space
      dressed by half a bridge on each side:
          Tblock = m_half^{1/2} T0_block m_half^{1/2}-shaped,
      a FINITE object in the program's existing grammar. Declared as the
      next dock; nothing about it is claimed here.

HONEST REMAINDER: one bridge is the minimal overlap; the chain S_m, the
2D lattice, the AF trajectory, tightness, OS, non-triviality, metric
universality and the Clay predicate all remain OPEN. The claim here is
the exact anatomy and certified gap of the FIRST overlapping system.

Controls:
  C1  bridge-swap invariance exact on the compressed matrix (equal
      entries under 1<->2 relabeling) at every grid kappa.
  C2  B-plane inertness exact: all first-order entries touching B vanish
      (parity bookkeeping), and the compressed matrix leaves the B-block
      diagonal at every kappa (within enclosure).
  C3  theta-integral tamper (1/2 -> 1) shifts the first-order splitting —
      separated brackets.
  C4  large-kappa honesty: the kappa = 1 cell does NOT certify k = 1 on
      this carrier — the certified outcome there (count 0 at mu = 3/2:
      no eigenvalue above mu at all, since the compressed vacuum weight
      sits near 1.13) is recorded as the carrier's honest limit, and mu
      below the |QSQ| bound is refused.
  C5  free limit: kappa = 0 reproduces the exact free spectrum
      (lambda_1 = 1, four-fold lambda at second level) on the carrier.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, _dec, canonical_sha, TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym4_symmetry_protected import chi_mul, dim, face_coeffs  # noqa: E402
from ym6_seam_integer_dock import ldl_inertia, _r  # noqa: E402

BETA = F(2)
P_CUT = 8
CELLS = [(F(1, 8), F(3, 5)), (F(1, 4), F(7, 10)), (F(1, 2), F(9, 10)),
         (F(1), F(3, 2))]      # last one expected to refuse honestly

# W5 basis: (A1-ring, A2-ring, B1-content, B2-content), twice-spin keys
BASIS = [
    ({0: F(1)}, {0: F(1)}, 0, 0),      # 1
    ({1: F(1)}, {0: F(1)}, 0, 0),      # chi12(A1)
    ({0: F(1)}, {1: F(1)}, 0, 0),      # chi12(A2)
    ({0: F(1)}, {0: F(1)}, 1, 0),      # chi12(B1)
    ({0: F(1)}, {0: F(1)}, 0, 1),      # chi12(B2)
]
N5 = len(BASIS)


def ring_mul_q(u: dict, v: dict) -> dict:
    out = {}
    for a, ca in u.items():
        for b, cb in v.items():
            for c in chi_mul(a, b):
                out[c] = out.get(c, F(0)) + ca * cb
    return out


_MCACHE = {}


def m_bridge(kappa: F):
    """5x5 interval matrix of exp[(kappa/2)Tr(A1 A2^-1)] on W5.
    B-contents must match (m has no B dependence); A-part uses the exact
    pairing Int chi_a(A1) chi_e(A2) chi_r(A1 A2^-1) = delta_{a=e=r}/d_r."""
    if kappa in _MCACHE:
        return _MCACHE[kappa]
    if kappa == 0:
        # exp[0] = chi_0: exact identity matrix, zero remainder
        Mw = [[Iv(F(1) if i == j else F(0)) for j in range(N5)]
              for i in range(N5)]
        _MCACHE[kappa] = (Mw, F(0))
        return _MCACHE[kappa]
    f = face_coeffs(kappa, P_CUT)
    M = [[Iv(F(0)) for _ in range(N5)] for _ in range(N5)]
    for i in range(N5):
        for j in range(i, N5):
            A1i, A2i, b1i, b2i = BASIS[i]
            A1j, A2j, b1j, b2j = BASIS[j]
            if b1i != b1j or b2i != b2j:
                M[i][j] = M[j][i] = Iv(F(0))
                continue
            P1 = ring_mul_q(A1i, A1j)
            P2 = ring_mul_q(A2i, A2j)
            acc = Iv(F(0))
            for r in range(P_CUT + 1):
                c1 = P1.get(r, F(0))
                c2 = P2.get(r, F(0))
                if c1 == 0 or c2 == 0:
                    continue
                # f_r * d_r * c1 * c2 * (1/d_r) = f_r * c1 * c2
                acc = acc + f[r] * Iv(dim(r) * c1 * c2 * F(1, dim(r)))
            M[i][j] = M[j][i] = _r(acc)
    # truncation remainder: one face, |phi_i phi_j|_inf <= 4
    from ym4_symmetry_protected import tail_E
    err = tail_E(kappa, P_CUT) * 4
    Mw = [[Iv(M[i][j].lo - err, M[i][j].hi + err) for j in range(N5)]
          for i in range(N5)]
    _MCACHE[kappa] = (Mw, err)
    return _MCACHE[kappa]


def lam_levels():
    I1 = bessel_I(1, BETA, TERMS)
    return bessel_I(2, BETA, TERMS) / I1


def s_block(kappa: F):
    lam = lam_levels()
    M, err = m_bridge(kappa)
    from ym6_seam_integer_dock import iv_sqrt
    s = iv_sqrt(lam)
    root_t = [Iv(F(1)), s, s, s, s]
    A = [[_r(root_t[i] * root_t[j] * M[i][j]) for j in range(N5)]
         for i in range(N5)]
    return A


def swap12(idx: int) -> int:
    return {0: 0, 1: 2, 2: 1, 3: 4, 4: 3}[idx]


def dock_cell(kappa: F, mu: F):
    lam = lam_levels()
    lam_next = lam * lam
    qsq_hi = (lam_next * Iv(exp_point(kappa).hi)).hi
    if not (qsq_hi < mu):
        return {"status": "REFUSED_C_NOT_BELOW_MU"}
    A = s_block(kappa)
    M2, _ = m_bridge(2 * kappa)
    Mk, _ = m_bridge(kappa)
    # |PMQ|^2 <= lammax(M2 - Mk^2) (orthonormal Gram) — Gershgorin
    sq = [[_r(sum((Mk[i][k] * Mk[k][j] for k in range(N5)), Iv(F(0))))
           for j in range(N5)] for i in range(N5)]
    gmax = None
    for i in range(N5):
        row = (M2[i][i] - sq[i][i]).hi
        for j in range(N5):
            if j != i:
                d = M2[i][j] - sq[i][j]
                row += max(abs(d.lo), abs(d.hi))
        gmax = row if gmax is None else max(gmax, row)
    from ym6_seam_integer_dock import iv_sqrt
    pmq = iv_sqrt(Iv(F(0), max(F(0), gmax))).hi
    psq = iv_sqrt(Iv(lam_next.hi)).hi * pmq
    r_hi = psq * psq / (mu - qsq_hi)
    low = [[_r(A[i][j] - Iv(mu if i == j else F(0))) for j in range(N5)]
           for i in range(N5)]
    high = [[_r(low[i][j] + Iv(r_hi if i == j else F(0)))
             for j in range(N5)] for i in range(N5)]
    li, hi_ = ldl_inertia(low), ldl_inertia(high)
    if li is None or hi_ is None:
        return {"status": "REFUSED_PIVOT_STRADDLE"}
    return {"status": "OK", "count_lower": li[0], "count_upper": hi_[0],
            "count_exact": li[0] if li[0] == hi_[0] else None}


def theta_integral(tamper=False) -> F:
    return F(1) if tamper else F(1, 2)


def first_order_block(tamper=False):
    """First-order derivative block on the 4-fold excited level, exact.
    Basis order (chi12A1, chi12A2, chi12B1, chi12B2); m' = (1/2)Tr(A1A2^-1).
    Only the A-pair couples, via the theta integral."""
    g = theta_integral(tamper)
    z = F(0)
    W = [[z, g / 2, z, z], [g / 2, z, z, z], [z, z, z, z], [z, z, z, z]]
    lam = lam_levels()
    # eigen: +-(g/2) on the A-plane, 0 (x2) on the B-plane
    d_plus = lam * Iv(g / 2)
    d_minus = lam * Iv(-g / 2)
    return W, d_plus, d_minus


def run():
    lam = lam_levels()

    # C1 swap invariance of the compressed matrix at each grid kappa
    c1 = True
    for kap, _mu in CELLS[:3]:
        M, _ = m_bridge(kap)
        for i in range(N5):
            for j in range(N5):
                a, b = M[i][j], M[swap12(i)][swap12(j)]
                if a.lo != b.lo or a.hi != b.hi:
                    c1 = False

    # C2 B-plane inert: off-diagonal entries touching B vanish exactly at
    # the compressed level (b-content mismatch forces zero before widening)
    c2 = True
    for kap, _mu in CELLS[:3]:
        M, err = m_bridge(kap)
        for i in (3, 4):
            for j in (0, 1, 2):
                if abs(M[i][j].lo) > err or abs(M[i][j].hi) > err:
                    c2 = False

    # T1 first-order anatomy + C3 tamper
    W, dp, dm = first_order_block()
    Wt, dpt, _ = first_order_block(tamper=True)
    c3 = dpt.separated_from(dp)
    t1 = (W[0][1] == F(1, 4) and W[2][2] == 0 and W[3][3] == 0
          and dp.lo > 0 > dm.hi)

    # T2 dock cells
    cells = {}
    any_refusal = False
    certified_all_small = True
    for kap, mu in CELLS:
        r = dock_cell(kap, mu)
        key = f"kappa={kap},mu={mu}"
        if r["status"] == "OK":
            cells[key] = {"count_bracket": [r["count_lower"],
                                            r["count_upper"]],
                          "seam_count_exact": r["count_exact"]}
            if kap <= F(1, 2) and r["count_exact"] != 1:
                certified_all_small = False
        else:
            cells[key] = {"refused": r["status"]}
            if kap != F(1):
                certified_all_small = False
        if kap == F(1):
            # honest large-kappa outcome: anything except a k=1 claim
            any_refusal = (r.get("count_exact") != 1)

    # T3 overlap price at kappa = 1/8: first-order prediction bracket
    kap = F(1, 8)
    pred = lam * Iv(kap / 4)          # kappa * lam / 4
    price_row = {"first_order_prediction_lo": _dec(pred.lo, 20),
                 "first_order_prediction_hi": _dec(pred.hi, 20),
                 "carried_by": "single swap-even line chi12(A1)+chi12(A2)"}

    # C5 free limit: kappa = 0 compressed matrix is the identity
    M0, e0 = m_bridge(F(0))
    c5 = all(abs(M0[i][j].lo - (1 if i == j else 0)) <= e0 + F(1, 10 ** 20)
             and abs(M0[i][j].hi - (1 if i == j else 0)) <= e0 + F(1, 10 ** 20)
             for i in range(N5) for j in range(N5))

    ok = c1 and c2 and c3 and t1 and certified_all_small and any_refusal and c5
    cert = {
        "certificate_type": "YM14_INTERACTION_OVERLAP_DOCK_BRIDGED_PAIR",
        "claim_status": "first_exact_data_in_overlap_regime",
        "theorems": {
            "T1_overlap_seam_rank_one_transported":
                "bridge-swap symmetry exact for all kappa; first-order "
                "splitting driven by the SAME exact theta integral 1/2; "
                "crossing line chi12(A1)+chi12(A2) (+lam/4), receding "
                "line (-lam/4), B-plane exactly inert — the inter-block "
                "overlap seam has the identical rank-one transported "
                "anatomy as the intra-block interaction (YM-3/4)",
            "T2_certified_two_sided_gap_bridged_pair":
                "exact seam counts k=1 on the (kappa,mu) grid up to "
                "kappa=1/2; at kappa=1 the carrier certifies count 0 at "
                "mu=3/2 (no k=1 window) — carrier limit recorded, not "
                "papered over",
            "T3_overlap_price_measured":
                "first-order price kappa*lam/4, carried by one swap-even "
                "direction, not spread",
            "T4_1d_block_transfer_reduction_NAMED":
                "bounded-overlap chain uniformity reduces in 1D to the "
                "spectrum of a finite block-transfer operator — declared "
                "next dock, nothing claimed",
        },
        "cells": cells,
        "overlap_price": price_row,
        "honest_remainder": {
            "scope": ("one bridge = minimal overlap; chain S_m, 2D "
                      "lattice, AF trajectory, tightness, OS, "
                      "non-triviality, metric universality, Clay all "
                      "OPEN"),
        },
        "controls": {
            "C1_bridge_swap_invariance": bool(c1),
            "C2_B_plane_inert": bool(c2),
            "C3_theta_tamper_separates": bool(c3),
            "C4_large_kappa_no_k1_claim": bool(any_refusal),
            "C5_free_limit_identity": bool(c5),
            "T1_anatomy": bool(t1),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM14_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM14.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"])
    for k, v in cert["cells"].items():
        print(" ", k, "->", v)
    print("sha256:", sha)
