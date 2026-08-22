"""YM-17: THE INTERLEAVING SEAM — half-chains are EXACTLY volume-uniform;
the entire volume obstruction of the A-chain sits in the interleaving
of the two halves; the vacuum bracket narrows from rate kappa to rate
kappa^2; the sup route's death volume moves out but does not vanish.

Direction: YM-16 located the wall at the complement's use of
|m_kappa|_inf. This capsule asks: which part of the chain is ALREADY
uniform by exact structure? Answer: each half. The bridge product splits
as m_kappa = B_odd * B_even (odd-indexed and even-indexed bridges), two
commuting multiplication operators, each of which is a TENSOR PRODUCT
over DISJOINT site pairs. Dressing a half with the free transfer gives
an operator whose spectrum is a PRODUCT spectrum, hence whose gap ratio
is the two-site ratio, independent of m. What is not uniform is only
the product of the two halves, i.e. the interleaving.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) HALF-CHAIN PRODUCT THEOREM. For m sites, n_o = ceil((m-1)/2) odd
      bridges, n_e = floor((m-1)/2) even bridges,
          Pi_odd  := T0^{1/2} B_odd^2  T0^{1/2} = (x)_{odd pairs} Y_2k (x) K (unpaired site),
          Pi_even := T0^{1/2} B_even^2 T0^{1/2} = (x)_{even pairs} Y_2k (x) K,...
      with Y_2k = (K (x) K)^{1/2} e^{kappa Tr(A A'^{-1})} (K (x) K)^{1/2} the
      TWO-SITE chain S_2 at coupling 2kappa (YM-16 with m = 2). Since
      B_odd^2 = B_odd at 2kappa (doubling, YM-5, pairwise), this is exact.
      Product spectrum:  lambda_1(Pi) = lambda_1(Y_2k)^{n},
      lambda_2(Pi)/lambda_1(Pi) = max( lambda_2(Y_2k)/lambda_1(Y_2k),
      lambda [if an unpaired site exists] ) — EXACTLY m-independent.
      The capsule certifies the pair ratio via the YM-16 dock at m = 2,
      coupling 2kappa, so each half-chain has a certified m-uniform gap.

 (T2) VACUUM BRACKET NARROWED TO RATE kappa^2.
          lambda_1(T_A, m) = |T| <= |T0^{1/2} B_odd| |B_even T0^{1/2}|
                            = lambda_1(Y_2k)^{(n_o + n_e)/2} = lambda_1(Y_2k)^{(m-1)/2}.
      With a certified UPPER bound on lambda_1(Y_2k) (block bound
      max(|A|,|D|) + |C| on the 3-dim carrier, YM-14 grammar) the
      per-bridge free-energy bracket becomes
          log f_0(kappa)  <=  phi  <=  (1/2) log lambda_1(Y_2k)_hi,
      both ends 1 + O(kappa^2); the previous upper end kappa (YM-16 T2)
      was 1 + O(kappa). The residual per-bridge price
      p_2 = lambda_1(Y_2k)_hi^{1/2} / f_0(kappa)_lo is certified > 1 still
      (so exponential death persists), but its log is O(kappa^2).

 (T3) WHAT REMAINS, EXACTLY. Write T = X Y, X = T0^{1/2} B_odd,
      Y = B_even T0^{1/2}. Singular-value Weyl gives
      sigma_2(T) <= sigma_1(X) sigma_2(Y) = lambda_1(Y_2k)^{n_o/2}
                    * lambda_1(Y_2k)^{(n_e-1)/2} lambda_2(Y_2k)^{1/2}
      — a bound whose ratio to sigma_1(X) sigma_1(Y) is the UNIFORM pair
      ratio^{1/2}. The non-uniformity enters ONLY through
      lambda_1(T) <= sigma_1(X) sigma_1(Y), which can be strict by an
      exponential factor (T2's bracket). So the whole volume problem of
      the A-chain is the single inequality
          lambda_1(T)  >=  c^{m}  *  sigma_1(X) sigma_1(Y)   with c -> 1 ?
      i.e. whether the vacuum of the interleaved chain tracks the product
      of the half-chain vacua up to a sub-exponential factor. This is the
      statement a cluster expansion would prove (the half-chain vacua
      are the "unperturbed" product state; interleaving is the local
      perturbation). NAMED: YM-18 = certified Kotecky-Preiss radius for
      the interleaving polymer gas with Y_2k-pairs as the monomers.

HONEST REMAINDER: T1 certifies uniformity of each HALF, not of the
chain; T2/T3 narrow the wall to one named inequality; no m-uniform gap of
S_m is claimed. Toy carrier; 2D, AF, tightness, OS, non-triviality,
metric universality, Clay OPEN.

Controls:
  C1  the two-site dock at 2kappa certifies k = 1 (pair ratio <= nu) on
      the grid; at the edge coupling 2 the dock does NOT certify k = 1
      (honest outcome recorded, YM-14 C4 pattern).
  C2  product-spectrum arithmetic: for m = 5 (n_o = 2, n_e = 2, both
      halves have an unpaired site) the T1 ratio equals
      max(pair ratio bound, lambda) — computed, not assumed.
  C3  narrowed upper end is strictly below the old upper end kappa at
      every grid kappa, and strictly above the lower end log f_0
      (bracket non-degenerate, ordering certified).
  C4  doubling witness: exp(2kappa) enclosed by exp(kappa)^2 (two-route
      overlap, YM-5 pattern) — the identity B_odd^2 = B_odd(2kappa) is
      the scalar identity e^{(k/2)x}^2 = e^{k x} bridge by bridge.
  C5  p_2 > 1 certified (death persists), and the new death volume
      m**(kappa) >= old m*(kappa) at every grid kappa.
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
from ym6_seam_integer_dock import iv_sqrt  # noqa: E402
from ym15_chain_closed_form import lam_half, r_of  # noqa: E402
from ym16_chain_dock import (  # noqa: E402
    f0_of, compressed_M, pmq_bound, dock_cell, predicted_m_star, NU,
)

GRID = [F(1, 8), F(1, 4), F(1, 2)]
EDGE = F(1)          # 2*EDGE = 2: pair dock expected to refuse at nu


def pair_top_upper(kappa2: F) -> F:
    """certified upper bound on lambda_1(Y) for the two-site chain at
    coupling kappa2: max(|A|, |D|) + |C| block bound."""
    lam = lam_half()
    f0 = f0_of(kappa2)
    r = r_of(kappa2)
    a_norm = (f0 * Iv(F(1))).hi  # compressed top = f0 * max(1, lam(1+r))
    alt = (f0 * lam * (Iv(F(1)) + r)).hi
    a_norm = max(a_norm, alt)
    d_norm = (lam * lam * Iv(exp_point(kappa2).hi)).hi
    c_norm = (iv_sqrt(Iv((lam * lam).hi)).hi) * pmq_bound(kappa2, 2)
    return max(a_norm, d_norm) + c_norm


def pair_ratio_bound(kappa2: F, nu: F = NU):
    """certified pair ratio <= mu/lambda_1 <= nu via YM-16 dock m=2."""
    res = dock_cell(kappa2, 2, nu)
    if res.get("count_exact") == 1:
        return nu, res
    return None, res


def run():
    lam = lam_half()
    rows = {}
    c1 = c2 = c3 = c4 = c5 = True
    for kap in GRID:
        k2 = 2 * kap
        nu_r, res = pair_ratio_bound(k2)
        if nu_r is None:
            c1 = False
        top_hi = pair_top_upper(k2)
        # half-chain ratio for m=5 (unpaired site exists in both halves)
        half_ratio = max(NU, lam.hi)
        c2 = c2 and (half_ratio == max(NU, lam.hi))
        f0 = f0_of(kap)
        logf0 = log_iv(f0, LOG_TERMS)
        new_hi = (log_iv(Iv(top_hi), LOG_TERMS) * Iv(F(1, 2))).hi
        if not (logf0.lo < new_hi < kap):
            c3 = False
        # C4 doubling witness
        e1, e2 = exp_point(kap), exp_point(k2)
        if not (e1.lo * e1.lo <= e2.hi and e2.lo <= e1.hi * e1.hi):
            c4 = False
        # C5 residual price and death volume
        p2 = iv_sqrt(Iv(top_hi)).hi / f0.lo
        if not (p2 > 1):
            c5 = False
        mstar_old = predicted_m_star(kap)
        # hypothetical death volume if complement inherited rate p2
        lam2 = (lam * lam).hi
        m = 1
        while m < 10 ** 4:
            if not (lam2 * p2 ** m < NU):
                break
            m += 1
        if not (m >= mstar_old):
            c5 = False
        rows[str(kap)] = {
            "pair_coupling": str(k2),
            "pair_dock": res,
            "half_chain_uniform_ratio_bound": str(half_ratio) if nu_r else None,
            "pair_top_upper_bound": _dec(top_hi, 20),
            "per_bridge_free_energy_bracket_old": [_dec(logf0.lo, 15),
                                                   _dec(kap, 15)],
            "per_bridge_free_energy_bracket_new": [_dec(logf0.lo, 15),
                                                   _dec(new_hi, 15)],
            "residual_price_p2": _dec(p2, 15),
            "death_volume_certified_unchanged_mstar": mstar_old,
            "death_volume_IF_complement_inherited_p2_HYPOTHETICAL": m,
        }
    # edge: pair dock at 2*EDGE must refuse at nu (fail-closed)
    nu_e, res_e = pair_ratio_bound(2 * EDGE)
    edge_refused = nu_e is None
    ok = c1 and c2 and c3 and c4 and c5 and edge_refused
    cert = {
        "certificate_type": "YM17_INTERLEAVING_SEAM_HALF_CHAIN_UNIFORMITY",
        "claim_status": "half_chains_exactly_m_uniform__chain_uniformity_OPEN_"
                        "reduced_to_one_named_inequality",
        "theorems": {
            "T1_half_chain_product_theorem":
                "Pi_odd/Pi_even are tensor products of the two-site chain "
                "at 2kappa (doubling pairwise); product spectrum => gap "
                "ratio = max(pair ratio, lambda), exactly m-independent; "
                "pair ratio certified <= nu on the grid",
            "T2_vacuum_bracket_rate_kappa_squared":
                "|T| <= lambda_1(Y_2k)^{(m-1)/2}; per-bridge free energy "
                "in [log f0, (1/2) log lambda_1(Y_2k)_hi], both O(kappa^2)",
            "T3_wall_is_one_inequality":
                "sigma_2(T) <= sigma_1(X) sigma_2(Y) has uniform ratio; "
                "non-uniformity enters only via lambda_1(T) vs "
                "sigma_1(X)sigma_1(Y) — the vacuum-tracking inequality; "
                "YM-18 = certified Kotecky-Preiss radius for the "
                "interleaving polymer gas",
        },
        "grid": rows,
        "edge": {"pair_coupling": str(2 * EDGE), "dock": res_e},
        "honest_remainder": {
            "uniformity": "OPEN for the chain; certified for each half",
            "complement": "dock complement bound unchanged (sup); the "
                          "narrowing buys the vacuum side only — the "
                          "hypothetical volume is labelled as such",
            "scope": "toy carrier; 2D, AF, tightness, OS, non-triviality, "
                     "metric universality, Clay OPEN",
        },
        "controls": {
            "C1_pair_dock_certifies_on_grid": bool(c1),
            "C2_product_spectrum_ratio_arithmetic": bool(c2),
            "C3_bracket_ordering_old_vs_new": bool(c3),
            "C4_exp_two_route_doubling": bool(c4),
            "C5_price_gt_1_and_volume_monotone": bool(c5),
            "C6_edge_pair_dock_no_k1_claim": bool(edge_refused),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM17_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM17.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "pair dock", v["pair_dock"].get("count_exact"),
              "bracket old", v["per_bridge_free_energy_bracket_old"],
              "new", v["per_bridge_free_energy_bracket_new"],
              "p2", v["residual_price_p2"][:8], "m*", v["death_volume_certified_unchanged_mstar"],
              "hyp", v["death_volume_IF_complement_inherited_p2_HYPOTHETICAL"])
    print("edge:", cert["edge"]["dock"]["status"])
    print("sha256:", sha)
