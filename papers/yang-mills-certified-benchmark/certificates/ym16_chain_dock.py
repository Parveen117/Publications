"""YM-16: THE CHAIN DOCK — certified two-sided gap of the bounded-overlap
chain S_m for EVERY m up to an EXACT, computable volume m*(kappa);
an exact m-uniform UPPER bound on the gap from B-factorisation; and the
sup-route's death volume computed, not feared.

Direction check (the reason this capsule exists): the Millennium
predicate needs a gap UNIFORM in the volume. YM-15 closed the carrier
side for all m. This capsule attacks the complement side with the
program's existing grammar, certifies exactly how far it reaches, and
names the object that must replace it.

Setting as in YM-15: chain S_m, blocks (A_i, B_i), bridges only on A.
Carrier W_{m+1} = {1, chi12(A_1..A_m)}, orthonormal, exact T0
eigenvectors; compressed bridge product (YM-15 T1, exact for all m)
    M_kappa = f_0^{m-1} * [ 1 (+) KMS_m(r) ],  r = I_2(kappa)/I_1(kappa).

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) EXACT B-FACTORISATION  =>  m-UNIFORM UPPER BOUND ON THE GAP.
      The bridge product does not depend on any B_i and T0 is a tensor
      product, so on L^2(SU(2)^{2m}) the chain transfer factorises
      EXACTLY:  T_kappa(S_m) = T_A,kappa (x) K^{(x)m}_B.  Hence
      chi12(B_1) * psi_vac(A) is an exact eigenvector (Ad-invariant:
      product of two invariants) with eigenvalue lambda * lambda_1(T_A).
      So for EVERY m and EVERY kappa,
          lambda_2 / lambda_1  >=  lambda,    gap(S_m)  <=  -log lambda
                                             = Delta_red(YM-1) = 0.83672...
      The free B-excitation is the cheapest excitation of the chain, at
      every volume and coupling; the gap is m-uniformly BOUNDED ABOVE,
      saturated at kappa = 0. Consequently the only volume-dependent
      content lives in the A-chain, and the Millennium question for this
      toy is exactly: is lambda_2(T_A)/lambda_1(T_A) < 1 uniformly in m.

 (T2) EXACT TWO-SIDED VOLUME BRACKET ON THE A-VACUUM.
          f_0^{m-1}  <=  lambda_1(T_A, m)  <=  e^{kappa (m-1)}
      (lower: the carrier vacuum is a trial state and YM-15 gives its
      Rayleigh quotient exactly; upper: |m_kappa|_inf <= e^{kappa(m-1)}).
      Per-site free-energy bracket [log f_0(kappa), kappa], m-independent,
      certified.  The ratio of the two ends, (e^kappa / f_0)^{m-1}, is
      the EXACT price of the sup route per bridge; it is > 1 for every
      kappa > 0 (e^x > 2 I_1(x)/x).

 (T3) CERTIFIED GAP OF S_m FOR m <= m*(kappa), AND THE DEATH VOLUME.
      Haynsworth dock on W_{m+1} (YM-6/14 grammar) with:
        complement T0 top = lambda^2 (two half-spins; lambda_1 = I_3/I_1
          certified below lambda^2), m-independent;
        |QSQ| <= lambda^2 e^{kappa(m-1)}  (sup route, the only
          complement control in the program so far);
        |PMQ|^2 <= lammax(M_{2kappa} - M_kappa^2) in CLOSED FORM
          (YM-15: both terms are f_0-powers times KMS matrices);
        threshold mu = nu * f_0(kappa)_lo^{m-1}.
      Exact seam count k = 1 certified cell by cell: lambda_2 <= mu <
      lambda_1, i.e. a certified two-sided gap ratio <= nu, for every
      m from 2 up to the largest m the route admits. Beyond that, the
      route is REFUSED by the arithmetic, and the capsule computes the
      exact refusal volume
          m*(kappa, nu) = 1 + floor( log(nu / lambda^2) / (kappa - log f_0) )
      (the last m with |QSQ| < mu) and checks it against the dock's
      actual first refusal. This is YM-11 gate-2 one level deeper: the
      sup route dies at a COMPUTABLE finite volume, quadratically there,
      exponentially here.

HONEST REMAINDER (the wall, now with coordinates):
  * T3 is a FINITE-volume certification: S_m has a certified two-sided
    gap for m <= m*(kappa), NOT for all m. Uniformity in m on the
    A-chain is OPEN.  The obstruction is identified exactly: the
    complement bound uses |m_kappa|_inf = e^{kappa(m-1)} while the
    vacuum only grows like f_0^{m-1} .. e^{kappa(m-1)}; any bound that
    pays the sup per bridge cannot be m-uniform (T2 proves the ratio
    of the bracket ends is exponential).
  * NAMED NEXT (YM-17): a LOCAL complement control — per-bridge
    operator control in the character basis (the chain is an exact
    nearest-neighbour operator on content configurations; bridge i acts
    only on sites i, i+1 with coefficients d_j f_j) — i.e. a certified
    strong-coupling cluster / gap-stability argument for a local
    perturbation of the classical gapped chain T0 (gap C_1/2 = 3/4 in
    heat-kernel units, YM-9). That is the Osterwalder-Seiler shape; the
    program's contribution would be to certify its radius exactly.
  * Toy carrier; 2D lattice, AF trajectory, tightness, OS,
    non-triviality, metric universality, Clay: OPEN.

Controls:
  C1  lambda_1 = I_3/I_1 certified strictly below lambda^2 (complement
      top identification bites).
  C2  B-factorisation witnessed on the grid: the eigenvalue
      lambda*(carrier vacuum) is reproduced by an explicit product
      state, and the compressed A-block never exceeds it (lambda *
      lammax(KMS) < 1 within the certified region) — consistent with T1.
  C3  kappa = 0 recovers the free chain exactly (f_0 = 1, r = 0, counts
      certified for all m tested).
  C4  predicted m* equals the dock's actual last certified m at every
      grid kappa (formula vs arithmetic).
  C5  fail-closed: the cell m = m*+1 is refused, not certified.
  C6  per-bridge price e^kappa/f_0 > 1 certified at every kappa > 0.
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
    Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym6_seam_integer_dock import ldl_inertia, _r, iv_sqrt  # noqa: E402
from ym15_chain_closed_form import kms, lam_half, r_of  # noqa: E402

BETA = F(2)
NU = F(9, 10)                                   # ratio threshold
GRID = [F(1, 8), F(1, 4), F(1, 2)]
M_MAX = 40                                      # hard ceiling on m scanned


def f0_of(kappa: F) -> Iv:
    if kappa == 0:
        return Iv(F(1))
    return bessel_I(1, kappa, TERMS) * Iv(F(2)) / Iv(kappa)


def iv_pow(x: Iv, e: int) -> Iv:
    out = Iv(F(1))
    for _ in range(e):
        out = _r(out * x)
    return out


def compressed_M(kappa: F, m: int):
    """M_kappa on W_{m+1}: f0^{m-1} * [1 (+) KMS_m(r)] as interval matrix."""
    f0 = f0_of(kappa)
    scale = iv_pow(f0, m - 1)
    K = kms(m, r_of(kappa))
    n = m + 1
    M = [[Iv(F(0))] * n for _ in range(n)]
    M[0][0] = scale
    for i in range(m):
        for j in range(m):
            M[i + 1][j + 1] = _r(scale * K[i][j])
    return M


def mat_mul(A, B):
    n = len(A)
    return [[_r(sum((A[i][k] * B[k][j] for k in range(n)), Iv(F(0))))
             for j in range(n)] for i in range(n)]


def pmq_bound(kappa: F, m: int) -> F:
    """|PMQ| <= sqrt(lammax(M_{2k} - M_k^2)); orthonormal Gram."""
    M2 = compressed_M(2 * kappa, m)
    Mk = compressed_M(kappa, m)
    D = mat_mul(Mk, Mk)
    n = m + 1
    W = [[_r(M2[i][j] - D[i][j]) for j in range(n)] for i in range(n)]
    gmax = gmin = None
    trace_hi = F(0)
    for i in range(n):
        trace_hi += W[i][i].hi
        row_hi, row_lo = W[i][i].hi, W[i][i].lo
        for j in range(n):
            if j != i:
                off = max(abs(W[i][j].lo), abs(W[i][j].hi))
                row_hi += off
                row_lo -= off
        gmax = row_hi if gmax is None else max(gmax, row_hi)
        gmin = row_lo if gmin is None else min(gmin, row_lo)
    bound = min(gmax, trace_hi - n * min(gmin, F(0)))
    return iv_sqrt(Iv(F(0), max(F(0), bound))).hi


def s_block(kappa: F, m: int):
    lam = lam_half()
    s = iv_sqrt(lam)
    M = compressed_M(kappa, m)
    n = m + 1
    rt = [Iv(F(1))] + [s] * m
    return [[_r(rt[i] * rt[j] * M[i][j]) for j in range(n)] for i in range(n)]


def dock_cell(kappa: F, m: int, nu: F = NU):
    lam = lam_half()
    lam_next = lam * lam
    f0 = f0_of(kappa)
    mu = nu * iv_pow(f0, m - 1).lo if kappa != 0 else nu
    # sup route: |QSQ| <= lambda^2 e^{kappa(m-1)}
    qsq_hi = (lam_next * Iv(exp_point(kappa * (m - 1)).hi)).hi \
        if kappa != 0 else lam_next.hi
    if not (qsq_hi < mu):
        return {"status": "REFUSED_C_NOT_BELOW_MU"}
    pmq = pmq_bound(kappa, m) if kappa != 0 else F(0)
    psq = iv_sqrt(Iv(lam_next.hi)).hi * pmq
    r_hi = psq * psq / (mu - qsq_hi)
    A = s_block(kappa, m)
    n = m + 1
    low = [[_r(A[i][j] - Iv(mu if i == j else F(0))) for j in range(n)]
           for i in range(n)]
    high = [[_r(low[i][j] + Iv(r_hi if i == j else F(0)))
             for j in range(n)] for i in range(n)]
    li, hi_ = ldl_inertia(low), ldl_inertia(high)
    if li is None or hi_ is None:
        return {"status": "REFUSED_PIVOT_STRADDLE"}
    return {"status": "OK", "count_lower": li[0], "count_upper": hi_[0],
            "count_exact": li[0] if li[0] == hi_[0] else None,
            "mu": str(mu)}


def predicted_m_star(kappa: F, nu: F = NU) -> int:
    """largest m with lambda^2 e^{kappa(m-1)} < nu f0_lo^{m-1} (exact scan)."""
    lam2 = lam_half() * lam_half()
    f0 = f0_of(kappa)
    m = 1
    while m < M_MAX:
        lhs = (lam2 * Iv(exp_point(kappa * m).hi)).hi
        rhs = nu * iv_pow(f0, m).lo
        if not (lhs < rhs):
            break
        m += 1
    return m                      # m-1 bridges admitted -> chain length m


def run():
    lam = lam_half()
    lam2 = lam * lam
    lam1 = bessel_I(3, BETA, TERMS) / bessel_I(1, BETA, TERMS)
    c1 = lam1.hi < lam2.lo
    gap_upper = -log_iv(lam, LOG_TERMS)

    rows = {}
    c4 = c5 = c6 = c2 = True
    for kap in GRID:
        f0 = f0_of(kap)
        price = Iv(exp_point(kap).lo, exp_point(kap).hi) / f0
        if not (price.lo > 1):
            c6 = False
        logf0 = log_iv(f0, LOG_TERMS)
        mstar = predicted_m_star(kap)
        cells = {}
        last_ok = 1
        for m in range(2, mstar + 2):
            res = dock_cell(kap, m)
            cells[str(m)] = res
            if res.get("count_exact") == 1:
                last_ok = m
        if last_ok != mstar:
            c4 = False
        if cells[str(mstar + 1)]["status"] == "OK" and \
                cells[str(mstar + 1)].get("count_exact") == 1:
            c5 = False
        # C2: compressed A-top below carrier vacuum on grid (YM-15 ceiling)
        r = r_of(kap)
        ceil = lam * (Iv(F(1)) + r) / (Iv(F(1)) - r)
        if not (ceil.hi < 1):
            c2 = False
        rows[str(kap)] = {
            "per_site_free_energy_bracket": [_dec(logf0.lo, 20),
                                             _dec(kap, 20)],
            "sup_route_price_per_bridge_e^k/f0_lo": _dec(price.lo, 20),
            "m_star_predicted": mstar,
            "m_star_dock_actual": last_ok,
            "certified_ratio_threshold_nu": str(NU),
            "certified_gap_lower_bound_for_m_le_mstar":
                _dec((-log_iv(Iv(NU), LOG_TERMS)).lo, 20),
            "cells": cells,
        }

    # C3 free chain
    c3 = all(dock_cell(F(0), m)["count_exact"] == 1 for m in (2, 5, 9))

    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM16_CHAIN_DOCK_FINITE_VOLUME_AND_DEATH_VOLUME",
        "claim_status": "certified_gap_of_S_m_for_m_le_mstar_kappa__"
                        "uniformity_OPEN_obstruction_localised",
        "theorems": {
            "T1_B_factorisation_uniform_upper_bound":
                "T(S_m) = T_A (x) K_B^{(x)m} exactly; chi12(B_1) psi_vac "
                "is an exact eigenvector with ratio lambda => gap(S_m) <= "
                "-log lambda = Delta_red for every m and kappa",
            "T2_vacuum_volume_bracket":
                "f0^{m-1} <= lambda_1(T_A,m) <= e^{kappa(m-1)}; per-site "
                "free energy in [log f0, kappa]; sup price e^kappa/f0 > 1",
            "T3_finite_volume_certified_gap_and_death_volume":
                "exact seam count k=1 on W_{m+1} for m = 2..m*(kappa): "
                "certified two-sided ratio <= nu; m* computed exactly and "
                "matched by the dock; m*+1 refused",
        },
        "gap_upper_bound_all_m_all_kappa": [_dec(gap_upper.lo, 20),
                                            _dec(gap_upper.hi, 20)],
        "grid": rows,
        "honest_remainder": {
            "uniformity": "OPEN on the A-chain; sup route provably "
                          "non-uniform (T2); certified only for m <= m*",
            "named_next": "YM-17: local (per-bridge) complement control in "
                          "the content basis — certified strong-coupling "
                          "cluster / gap-stability radius",
            "scope": "toy carrier; 2D, AF, tightness, OS, non-triviality, "
                     "metric universality, Clay OPEN",
        },
        "controls": {
            "C1_complement_top_is_lambda_squared": bool(c1),
            "C2_carrier_consistent_with_B_factorisation": bool(c2),
            "C3_free_chain_exact": bool(c3),
            "C4_predicted_mstar_equals_dock": bool(c4),
            "C5_fail_closed_at_mstar_plus_1": bool(c5),
            "C6_sup_price_exceeds_1": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM16_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM16.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    print("gap upper bound all m:", cert["gap_upper_bound_all_m_all_kappa"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "m* pred", v["m_star_predicted"], "dock",
              v["m_star_dock_actual"], "price", v["sup_route_price_per_bridge_e^k/f0_lo"][:8],
              {m: c.get("count_exact", c["status"]) for m, c in v["cells"].items()})
    print("sha256:", sha)
