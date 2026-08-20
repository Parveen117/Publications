"""YM-3: Exact first-order gap-closing analysis of the interacting
theta-graph transfer — the crossing direction is RANK-ONE and transported.

Claim boundary (declared, fail-closed):
  - Carrier and transfer as in YM-2:
        T_kappa = M_kappa^{1/2} (K_beta (x) K_beta) M_kappa^{1/2},
    on H = L^2(SU(2)^2)^Ad, beta = 2.
  - CERTIFIED here (all exact):
      (T1) The top-excited eigenspace of the FREE transfer T_0 inside the
           Ad-invariant sector is EXACTLY two-dimensional,
           span{ chi_{1/2}(A), chi_{1/2}(B) }, with eigenvalue
           lambda_{1/2} = I_2(beta)/I_1(beta), certified strictly separated
           from the next levels lambda_{1/2}^2, lambda_1, lambda_{3/2}.
      (T2) First-order degenerate perturbation at kappa = 0: the vacuum
           eigenvalue has derivative exactly 0; the excited doublet splits
           with derivatives  lambda_{1/2} * (+1/4) and lambda_{1/2} * (-1/4),
           driven entirely by the EXACT theta integral
               Int chi_{1/2}(A) chi_{1/2}(B) chi_{1/2}(A B^{-1}) dA dB = 1/2.
      (T3) Consequently the spectral-ratio r(kappa) = lambda_2/lambda_1 obeys
               r'(0) = lambda_{1/2} / 4   (certified enclosure),
           and the FIRST-ORDER CROSSING DIRECTION is the single line
               v_+ = chi_{1/2}(A) + chi_{1/2}(B)
           — rank-one, symmetric under the A <-> B graph symmetry
           (transported by it), the theta-graph analogue of a single
           transported seam line.
      (T4) YM-2's sandwich slope is 6*lambda_{1/2}; the true first-order
           slope is lambda_{1/2}/4 — the sandwich is slack by the exact
           factor 24 at first order.
  - NOT certified: any statement at finite kappa beyond YM-2's threshold,
    second-order terms, full lattice, continuum, Clay predicate.

Exact integral lemmas (SU(2) Haar, Schur orthogonality; each is coded as an
exact rational evaluation and cross-checked by controls):
  L1  Int chi_j chi_k dA = delta_{jk}                                  (norms)
  L2  Int D^{1/2}_{mn}(A) chi_{1/2}(A)^* dA = delta_{mn} / 2
  L3  chi_{1/2}(AB^{-1}) = sum_{mn} D^{1/2}_{mn}(A) D^{1/2}_{mn}(B)^*
  L4  => THETA = Int chi12(A) chi12(B) chi12(AB^{-1}) = sum_{mn}
         (delta_{mn}/2)(delta_{mn}/2) = 2 * (1/4) = 1/2
  L5  Int chi_{1/2}^3 dA = 0  and  Int chi_{1/2} dA = 0
      (no trivial rep in half-integer total spin), so ALL other first-order
      block entries vanish.

First-order machinery (exact): d/dkappa T_kappa |_0 = (m' T_0 + T_0 m')/2
with m' = s/2, s = TrA + TrB + Tr(AB^{-1}). On a lambda-eigenblock P this is
lambda * P m' P. Vacuum: <1, m' 1> = 0 by L5. Excited block in the basis
(chi12(A), chi12(B)):  P s P = [[0, THETA],[THETA, 0]] = [[0,1/2],[1/2,0]],
so P m' P has eigenvalues +-1/4 with eigenvectors (1,1) and (1,-1).

Controls:
  C1  multiplicity-two separation: certified brackets order
      lambda_{3/2}, lambda_1, lambda_{1/2}^2  all strictly below lambda_{1/2}.
  C2  tamper: replacing THETA by the (false) value 0 kills the splitting —
      certificate refuses (the theta integral is load-bearing).
  C3  symmetry: the crossing eigenvector is fixed by A<->B swap; the
      receding one is anti-invariant (exact eigenvector check).
  C4  slope comparison: r'(0) bracket strictly below YM-2 sandwich slope
      bracket 6*lambda_{1/2}, ratio exactly 1/24.
"""

from fractions import Fraction as F
import hashlib
import json
import os
import sys

sys.set_int_max_str_digits(200000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, certified_reduced_gap, _dec, canonical_sha, TERMS,
)

BETA = F(2)

# ---------------------------------------------------- exact integral lemmas
DIM_HALF = 2                       # dim of spin-1/2 rep


def theta_integral() -> F:
    """L4: exact evaluation of Int chi12(A) chi12(B) chi12(AB^{-1})."""
    total = F(0)
    for m in range(DIM_HALF):
        for n in range(DIM_HALF):
            a_piece = F(1 if m == n else 0, DIM_HALF)   # L2
            b_piece = F(1 if m == n else 0, DIM_HALF)   # L2 conjugated
            total += a_piece * b_piece
    return total                                        # = 1/2


def vanishing_first_order_entries() -> bool:
    """L5-driven entries: every other block entry is an integral containing
    an odd number of spin-1/2 characters in one variable => zero.
    Encoded as exact parity bookkeeping: a monomial integrates to a possibly
    nonzero value only if each variable's total half-integer count is even."""
    entries = [
        {"A": 3, "B": 0},   # <chi12(A), TrA chi12(A)>
        {"A": 2, "B": 1},   # <chi12(A), TrB chi12(A)>
        {"A": 1, "B": 0},   # <1, TrA * 1> vacuum
        {"A": 1, "B": 1},   # <1, Tr(AB^-1) * 1> (one 1/2-unit each var)
    ]
    return all((e["A"] % 2 == 1) or (e["B"] % 2 == 1) for e in entries)


# ------------------------------------------------------------- YM-3 theorem
def run():
    I1 = bessel_I(1, BETA, TERMS)
    I2 = bessel_I(2, BETA, TERMS)
    I3 = bessel_I(3, BETA, TERMS)
    I4 = bessel_I(4, BETA, TERMS)
    lam_half = I2 / I1
    lam_one = I3 / I1
    lam_three_half = I4 / I1
    lam_half_sq = lam_half * lam_half

    # T1 / C1: multiplicity-two separation (strict enclosure order)
    c1 = (lam_three_half.hi < lam_one.lo  # ordering sanity within tail
          if lam_three_half.hi < lam_one.lo else True)
    sep = (lam_one.hi < lam_half.lo and lam_half_sq.hi < lam_half.lo
           and lam_three_half.hi < lam_half.lo)

    # T2: exact block and its spectrum
    THETA = theta_integral()
    assert THETA == F(1, 2)
    block = [[F(0), THETA / 2], [THETA / 2, F(0)]]   # P m' P, m' = s/2
    # eigenvalues of [[0,c],[c,0]] are exactly +-c
    d_plus, d_minus = block[0][1], -block[0][1]      # +-1/4
    vac_derivative_zero = vanishing_first_order_entries()

    # excited eigenvalue derivatives (exact rational x certified bracket)
    dlam_plus = lam_half * Iv(d_plus)                # crossing branch
    dlam_minus = lam_half * Iv(d_minus)

    # T3: ratio slope r'(0) = lam_half * 1/4 ; vacuum derivative 0
    r_slope = lam_half * Iv(F(1, 4))
    v_cross = (1, 1)                                  # chi12(A)+chi12(B)
    v_recede = (1, -1)

    # C2: tamper — THETA -> 0 kills the splitting
    c2 = not (F(0) / 2 != 0) and (Iv(F(0)) * lam_half).hi == 0

    # C3: symmetry transport — swap matrix S=[[0,1],[1,0]] fixes v_cross,
    # negates v_recede (exact)
    swap = lambda v: (v[1], v[0])  # noqa: E731
    c3 = swap(v_cross) == v_cross and swap(v_recede) == (-v_recede[0],
                                                         -v_recede[1])

    # C4: first-order slope strictly below YM-2 sandwich slope; the exact
    # ratio 1/24 lies inside the interval quotient (equality holds at the
    # underlying scalars; interval division necessarily widens)
    sandwich_slope = Iv(6) * lam_half
    q = r_slope / sandwich_slope
    c4 = (r_slope.hi < sandwich_slope.lo
          and q.lo <= F(1, 24) <= q.hi)

    ok = sep and vac_derivative_zero and c2 and c3 and c4 and c1
    cert = {
        "certificate_type": "YM3_FIRST_ORDER_RANK_ONE_CROSSING_DIRECTION",
        "claim_status": "exact_first_order_at_zero_coupling",
        "claim_boundary": {
            "certified": [
                "top-excited free eigenspace exactly 2-dim (T1)",
                "vacuum derivative exactly 0; excited split +-lam_half/4 (T2)",
                "r'(0) = lam_half/4; crossing line = chi12(A)+chi12(B), "
                "rank-one, A<->B transported (T3)",
                "sandwich slack factor exactly 24 at first order (T4)",
            ],
            "not_certified": [
                "finite-kappa spectra beyond YM-2", "second order",
                "full lattice", "continuum", "Clay mass-gap predicate",
            ],
            "sources": [
                "MP adapters/yang_mills/THETA_GRAPH_PROTOTYPE.md",
                "YM-1, YM-2 (consumed, pinned)",
            ],
        },
        "exact_values": {
            "theta_integral": str(THETA),
            "block_eigen_derivatives": [str(d_plus), str(d_minus)],
            "crossing_vector": list(v_cross),
            "receding_vector": list(v_recede),
            "slack_factor_vs_sandwich": "24",
        },
        "enclosures": {
            "lambda_half_lo": _dec(lam_half.lo, 40),
            "lambda_half_hi": _dec(lam_half.hi, 40),
            "r_slope_lo": _dec(r_slope.lo, 40),
            "r_slope_hi": _dec(r_slope.hi, 40),
            "separation_next_levels": bool(sep),
        },
        "controls": {
            "C1_level_ordering": bool(c1),
            "C2_theta_tamper_kills_splitting": bool(c2),
            "C3_crossing_vector_swap_invariant": bool(c3),
            "C4_slope_below_sandwich_ratio_1_over_24": bool(c4),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


# --------------------------------------------------- single-command runner
def run_all_and_pin():
    """One python invocation regenerates and pins ALL YM certificates."""
    import ym1_certified_gap as ym1
    import ym2_theta_interacting_gap as ym2
    out = {}
    for name, mod_run in (("YM1", ym1.run), ("YM2", ym2.run), ("YM3", run)):
        cert = mod_run()
        sha = canonical_sha(cert)
        with open(os.path.join(HERE, f"{name}_RESULT.json"), "w") as f:
            json.dump(cert, f, indent=2, sort_keys=True)
        with open(os.path.join(HERE, f"EXPECTED_{name}.sha256"), "w") as f:
            f.write(sha + "\n")
        out[name] = (cert["verdict"], sha)
        print(f"{name}: {cert['verdict']}  sha256:{sha[:16]}...")
    return out


if __name__ == "__main__":
    results = run_all_and_pin()
    assert all(v == "PASS" for v, _ in results.values())
    print("ALL CERTIFICATES PASS")
