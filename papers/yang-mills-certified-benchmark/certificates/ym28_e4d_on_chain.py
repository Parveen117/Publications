"""YM-28: T01-E4D ON THE CHAIN — the mean-multiplier law applied to the
content-1/2 memory channel that theorum/28 hypothesis 3 asks about.

Source: RH-Framework T01-E4D (PR #19, PROVED, NATIVE_DERIVED):
   I_rec(h f) = <N(h)> I_rec(f) + Cov_nu(N(h), N(f))          (E4D-A, exact)
   I_rec(h f) <= <N(h)> I_rec(f) + osc N(h) * D_rec(f)        (E4D-B)
and its OPEN item E4D-C: control of the deviation along an ITERATED
multiplier.  Here the iterated multiplier is the face product
W = prod_l w(H_l) of the chain (YM-F1), N(W) = prod w^2, so <N(W)> is the
chain vacuum at doubled coupling: <w^2> = f_0(2kappa) per face (exact, the
w^2 = sum d_j f_j(2kappa) chi_j expansion).  States f are content-1/2
insertions chi_half(A_i), i in S, the multi-insertion channel of YM-26 T3.
Everything below is an exact formal polynomial identity in the f_j(2kappa)
evaluated by the YM-15 chain-integration engine (extended to content 8).

 (T1) E4D-C IS EXACT WITH ZERO DEVIATION ON SINGLE INSERTIONS. For every
      m and every site i,
          I_rec(W chi_half(A_i)) = f_0(2kappa)^{m-1} = <N(W)> * I_rec(chi_half(A_i)),
      i.e. the E4D-A covariance vanishes identically along the whole
      iteration: N(chi_half) = 1 + chi_1 and the chi_1 branch dies at the
      open end (YM-15).  The mean law holds exactly face after face,
      with no deviation growth in m.  Certified m = 2..8 as a formal
      identity (no f_1 monomial survives).

 (T2) TWO INSERTIONS: THE DEVIATION IS A GEOMETRIC CORRELATOR. For sites
      i < j,
          I_rec(W chi_half(A_i) chi_half(A_j)) / <N(W)> = 1 + r_1(2kappa)^{j-i},
      r_1 = f_1/f_0 — the chi_1-chi_1 correlator along the tiling (YM-21
      T2 at content 1 and doubled coupling).  Deviation decays
      geometrically in the separation; certified exactly, m = 2..8.

 (T3) k INSERTIONS: THE DEVIATION IS PER-INSERTION, NOT PER-FACE. For
      every insertion set S with |S| <= 4 and m <= 7,
          I_rec(W chi_S) / (<N(W)> I_rec(chi_S))  <=  (1 + 2 r_1(2kappa))^{|S|-1}.
      The iterated multiplier's deviation is bounded by a factor per
      INSERTION; the number of faces m does not enter.  At kappa = 1/8
      the factor is 1 + 2 r_1(1/4) = 1.0052.  This is the E4D-C control
      for the content-1/2 multi-insertion channel on its basis states —
      uniform in m, growing only with the excitation count.

 (T4) E4D EXPLAINS THE MEASURED SLACK. The per-face price of the A-line
      memory channel under the mean law, relative to the vacuum floor
      lambda_1 >= f_0(kappa)^{m-1}, is
          phi_E4D = log( sqrt(f_0(2kappa)) / f_0(kappa) )   (per face),
      against the T01-E4C sup price log(e^kappa / f_0(kappa)).  Certified:
      kappa=1/8: phi_E4D = 0.001949 vs YM-18's MEASURED true rate 0.002177
      vs E4C 0.1230;  kappa=1/4: 0.007742 vs 0.008698 vs 0.2422;
      kappa=1/2: 0.03016 vs 0.03461 vs 0.4689.  The mean law accounts
      for 90 / 89 / 87 percent of the true per-face rate; the sup price exceeds the
      E4D price by 63 / 31 / 16 x (YM-27 measured 56.5 / 27.8 / 13.5
      against the true rate).  The slack of YM-27 is now a derived
      number, not a measured mystery.

 (T5) CONDITIONAL PROJECTION (clearly labelled, not a theorem). If the
      outward memory bound of the theorum/28 dock were taken at the E4D
      rate on the whole content-1/2 channel (with T3's per-insertion
      factor at |S| <= 2), the chain length admitted by the certificate
      would be m*_E4D >= 400 (cap) / 200 / 50, computed exactly below,
      against the sup-route m* = 13 / 7 / 4.  This requires E4D-C on
      SUPERPOSITIONS of basis states (a quadratic-form statement), which
      T1-T3 do not give; it is the precise remaining obligation.

NOT CLAIMED: m-uniform gap; E4D-C on superpositions / the full memory
channel; anything about weak coupling.

Controls:
  C1  single insertion: formal identity, no f_1 monomial, m = 2..8.
  C2  two insertions: ratio == 1 + r_1^{j-i} exactly.
  C3  |S| <= 4, m <= 7: ratio <= (1 + 2 r_1)^{|S|-1}; tamper (1 + r_1)^{|S|-1}
      is violated at |S| = 3 (so the bound is not slack by accident).
  C4  E4D per-face rate within 15% of YM-18's measured rate at all three
      kappa, and below it (a lower bound on the true rate).
  C5  E4C price / E4D price reproduces YM-27's slack factors within 20%.
  C6  engine tamper (keep d_c factor) breaks C1.
"""

from fractions import Fraction as F
import itertools
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS  # noqa: E402
from ym4_symmetry_protected import chi_mul  # noqa: E402
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402
from ym16_chain_dock import predicted_m_star, NU, f0_of  # noqa: E402
from ym15_chain_closed_form import lam_half  # noqa: E402

P = 8                                   # formal f_0 .. f_8 (content up to spin 4)
GRID = [F(1, 8), F(1, 4), F(1, 2)]
YM18_PHI = {F(1, 8): F("0.002177"), F(1, 4): F("0.008698"), F(1, 2): F("0.03461")}
YM27_SLACK = {F(1, 8): F("56.5"), F(1, 4): F("27.8"), F(1, 2): F("13.5")}


# ------------------------------------------------ formal ring over f_0..f_P
def p_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, F(0)) + v
        if out[k] == 0:
            del out[k]
    return out


def p_mul_f(a, t):
    out = {}
    for k, v in a.items():
        kk = list(k)
        kk[t] += 1
        out[tuple(kk)] = v
    return out


def chain_entry(m, ins, tamper=False):
    state = {0: {tuple([0] * (P + 1)): F(1)}}
    for i in range(1, m + 1):
        t = ins.get(i, 0)
        nstate = {}
        for c, poly in state.items():
            for c2 in chi_mul(c, t):
                nstate[c2] = p_add(nstate.get(c2, {}), poly)
        state = nstate
        if i == m:
            return state.get(0, {})
        nstate = {}
        for c, poly in state.items():
            if c > P:
                raise ValueError("raise P")
            q = p_mul_f(poly, c)
            if tamper:
                q = {k: v * (c + 1) for k, v in q.items()}
            nstate[c] = p_add(nstate.get(c, {}), q)
        state = nstate
    return state.get(0, {})


def evaluate(poly, k2: F) -> Iv:
    """evaluate at f_t = I_{t+1}(k2) enclosures."""
    vals = [bessel_I(t + 1, k2, TERMS) for t in range(P + 1)]
    tot = Iv(F(0))
    for mono, c in poly.items():
        v = Iv(c)
        for t, e in enumerate(mono):
            for _ in range(e):
                v = v * vals[t]
        tot = tot + v
    return _r(tot)


def n_state_integral(m, S, k2, tamper=False):
    """I_rec(W chi_S) = Int prod w^2 prod_{i in S} chi_half(A_i)^2
       = sum over subsets of S of the chi_1-insertion chain entries."""
    tot = Iv(F(0))
    for n in range(len(S) + 1):
        for sub in itertools.combinations(S, n):
            tot = tot + evaluate(chain_entry(m, {i: 2 for i in sub}, tamper), k2)
    return _r(tot)


def run():
    c1 = c2 = c3 = c6 = True
    tamper_bound_fails = False
    table = {}
    for kap in GRID:
        k2 = 2 * kap
        f0 = bessel_I(1, k2, TERMS)
        f1 = bessel_I(3, k2, TERMS)
        r1 = _r(f1 / f0)
        rows = {}
        for m in range(2, 9):
            vac_poly = chain_entry(m, {})
            vac = evaluate(vac_poly, k2)
            # T1: single insertion formal identity
            for i in range(1, m + 1):
                poly = p_add(chain_entry(m, {}), chain_entry(m, {i: 2}))
                if poly != vac_poly:
                    c1 = False
                if any(mono[2] > 0 for mono in chain_entry(m, {i: 2})) and chain_entry(m, {i: 2}) != {}:
                    c1 = False
            # T2: two insertions
            for i in range(1, m + 1):
                for j in range(i + 1, m + 1):
                    ratio = n_state_integral(m, (i, j), k2) / vac
                    want = Iv(F(1)) + r1 ** (j - i) if hasattr(r1, "__pow__") else None
                    # exact check at polynomial level:
                    pij = chain_entry(m, {i: 2, j: 2})
                    mono = tuple([m - 1 - (j - i)] + [0] + [j - i] + [0] * (P - 2))
                    if pij != {mono: F(1)}:
                        c2 = False
            # T3: |S| <= 4, m <= 7
            if m <= 7:
                for size in (2, 3, 4):
                    for S in itertools.combinations(range(1, m + 1), size):
                        ratio = n_state_integral(m, S, k2) / vac
                        bound = Iv(F(1)) + Iv(F(2)) * r1
                        b = Iv(F(1))
                        bt = Iv(F(1))
                        for _ in range(size - 1):
                            b = b * bound
                            bt = bt * (Iv(F(1)) + r1)
                        if not (ratio.hi <= b.lo):
                            c3 = False
                        if size == 3 and ratio.lo > bt.hi:
                            tamper_bound_fails = True
            rows[str(m)] = "ok"
        # C6 engine tamper (two-insertion entry changes)
        if chain_entry(3, {1: 2, 2: 2}, tamper=True) == chain_entry(3, {1: 2, 2: 2}):
            c6 = False
        # T4 per-face rates with the chain's normalised f_0 (YM-16 f0_of)
        f0k = f0_of(kap)
        f0_2k = f0_of(k2)
        sqrt_f0_2k = None
        lo, hi = F(0), f0_2k.hi
        for _ in range(80):
            mid = (lo + hi) / 2
            if mid * mid < f0_2k.lo:
                lo = mid
            else:
                hi = mid
        sqrt_f0_2k = Iv(lo, hi)
        phi_e4d = log_iv(_r(sqrt_f0_2k / f0k), LOG_TERMS)
        e = exp_point(kap)
        phi_e4c = log_iv(_r(Iv(e.lo, e.hi) / f0k), LOG_TERMS)
        table[str(kap)] = {
            "r_1_at_2kappa": _dec(r1.lo, 8),
            "per_insertion_factor_1+2r1": _dec((Iv(F(1)) + Iv(F(2)) * r1).hi, 8),
            "phi_E4D_per_face": _dec(phi_e4d.lo, 8),
            "phi_measured_YM18": str(YM18_PHI[kap]),
            "phi_E4C_sup_per_face": _dec(phi_e4c.lo, 8),
            "E4D_share_of_true_rate": _dec(phi_e4d.lo / YM18_PHI[kap], 4),
            "E4C_over_E4D": _dec(phi_e4c.lo / phi_e4d.hi, 3),
            "m_star_sup_route": predicted_m_star(kap),
        }
        # T5 conditional projection: admit m while
        # lambda^2 * (sqrt f0(2k))^{m-1} * (1+2r1) < nu * f0(k)^{m-1}
        lam2 = (lam_half() * lam_half()).hi
        mm = 1
        while mm < 400:
            lhs = lam2 * sqrt_f0_2k.hi ** mm * (1 + 2 * r1.hi)
            rhs = NU * f0k.lo ** mm
            if not (lhs < rhs):
                break
            mm += 1
        table[str(kap)]["m_star_E4D_projected_CONDITIONAL"] = (
            ">=400 (cap)" if mm >= 400 else mm)
    c4 = all(F(table[str(k)]["E4D_share_of_true_rate"]) > F(85, 100) and
             F(table[str(k)]["phi_E4D_per_face"]) < YM18_PHI[k] for k in GRID)
    c5 = all(abs(F(table[str(k)]["E4C_over_E4D"]) - YM27_SLACK[k]) / YM27_SLACK[k] < F(1, 5)
             for k in GRID)
    ok = c1 and c2 and c3 and tamper_bound_fails and c4 and c5 and c6
    cert = {
        "certificate_type": "YM28_T01_E4D_ON_THE_CHAIN",
        "source": "RH-Framework T01-E4D (PR #19) applied to the content-1/2 memory channel",
        "claim_status": "mean_law_exact_on_single_insertions__deviation_per_insertion_"
                        "on_basis_states__E4D_explains_measured_slack__superposition_"
                        "statement_OPEN",
        "grid": table,
        "controls": {
            "C1_single_insertion_formal_identity": bool(c1),
            "C2_two_insertion_ratio_exact": bool(c2),
            "C3_per_insertion_bound_S_le_4": bool(c3),
            "C3b_tighter_tamper_bound_violated": bool(tamper_bound_fails),
            "C4_E4D_rate_explains_measured_rate": bool(c4),
            "C5_reproduces_YM27_slack": bool(c5),
            "C6_engine_tamper_bites": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM28_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM28.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, v)
    print("sha256:", sha)
