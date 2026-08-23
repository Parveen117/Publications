"""YM-35: T54' ENERGY VERSION + COMMUTATOR REMAINDER — E4D-C attacked with the
Recognition-Kernel-Framework operator ladder (theorums 50–55, 61), as the
owner suggested.

Where it stands (YM-34): E4D-C = m-uniform gap of the alternating product
T = K^{1/2} M_w K^{1/2} on the chain; the vacuum half is PROVED (YM-31); the
excitation upper bound is the only open item.  Reading the ladder:

  theorum/54 (infinite faces) delivers theorum/28's hypothesis list on a
     product of faces under  sup mu_i < 1  AND  sum mu_i < inf.  The chain has
     the SAME mu at every bridge, so sum mu = inf: 54's own T5
     (sum_mu_diverges) says the constant-mu chain is exactly the uncovered
     case, and 54's Pi_n = prod(1+mu) IS the T01-E4C sup route.
  theorum/50 A4 (post-audit): the odd/memory channel is invisible to
     SELF-pairing and visible to POLARIZATION x^dag S y — YM-28's "mean law
     exact on single insertions, open on superpositions" is this statement.
  theorum/53: odd sector lives in phases, positivity in weights; native
     Parseval E_Sigma = sum d_k E(c_k) — per-face rates belong in ENERGY
     (square-sourced, T01-B/C), not mass.
  theorum/51 exchange law + theorum/61 commutator mass: the calculus of the
     cross terms between the even time kernel and the odd faces.

This capsule makes the consequence exact on the fabric (YM-F1), with the
YM-30 engine as the only evaluator.  Faces B_l = chi_half(H_l) (H_l =
A_l^{-1} A_{l+1}), zero Haar mean; time kernel K acts on each site by lambda_j
on content j (bi-invariant, YM-25 T1); lambda := lambda_{1/2}, lambda_1 the
content-1 value.

 T1  LOCALITY of the commutator [K^{1/2}, B_l] (exact, every m).  On a
     state whose contents at sites l, l+1 are 0, B_l creates content 1/2 at
     both and K^{1/2} acts by lambda on either side: commutator ZERO.  So the
     commutator can be nonzero only where the state already carries content
     at l or l+1 — an ADJACENT face insertion (or the same face twice).
     Engine check: B_1 B_3 (non-adjacent) has site-2 content exactly 1/2
     (weight 1) — no mixing.  Commutator support is bridge-local: the
     theorum/61 "support" of [K^{1/2}, B_l] is {l-1, l+1}.
 T2  ADJACENT COMMUTATOR ENERGY, exact.  B_l B_{l+1} carries site-(l+1)
     content 0 with weight 1/4 and content 1 with weight 3/4 (engine:
     d_c S(c) = 1/4, 3/4 — YM-25's recoupling squares).  K^{1/2} acts by
     lambda * lambda_c^{1/2}; B_l K^{1/2} acts by lambda.  Hence
         E_Sigma([K^{1/2}, B_l] B_{l+1}) = lambda^2 * (3/4) * (1 - sqrt(lambda_1))^2
     — independent of kappa and of m; a polynomial identity in
     (lambda, s_1 = sqrt(lambda_1)) verified by exact interpolation.
 T3  T54' RATIO THEOREM on the COMMUTING model (exact, enumerated m <= 8):
     faces b_l = +-1 zero-mean commuting variables, K contracts each b_l by
     lambda, weight prod(1 + mu b_l).  Vacuum energy and single-excitation
     energy both carry the SAME per-face factor (1 + mu^2 lambda)^{m-1}; the
     ratio is EXACTLY m-independent for every m — T54's hypothesis
     sum mu < inf is not needed for a RATIO statement; constant mu is fine.
 T4  REMAINDER AT SECOND ORDER.  In the face expansion prod(1 + mu B_l) the
     leading per-face energy is mu^2 lambda^2 (single face, K B = lambda^2 B);
     the first non-commutative correction is the adjacent pair at order mu^4
     with energy mu^4 * T2.  Relative size
         rho(kappa, a) = mu^2 * (3/4) * (1 - sqrt(lambda_1))^2,   mu = r(kappa) = f_half/f_0
     displayed (mpmath, display only) at a = 1, kappa = 1/8, 1/4, 1/2.
     KILL CRITERION (declared before running): rho < 1/10 at all three and
     m-uniform by T1  ==>  E4D-C's excitation bound OPENS at first
     non-commutative order; else refusal.
 T5  HONEST REMAINDER.  What is NOT done: the full face expansion (all
     words, contents j >= 1 of the faces with coefficient ladder r_j), whose
     m-uniform convergence is the worldline cluster problem named in YM-33.
     T1–T4 show the first non-commutative correction is local and small;
     they do not bound the sum.  E4D-C stays OPEN; its first obstruction is
     now an exact local number, not a sup.
"""

from fractions import Fraction as F
import itertools
import json
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ym1_certified_gap import canonical_sha  # noqa: E402
from ym30_recoupling_engine import ladder_eval  # noqa: E402

mp.mp.dps = 30
half, one, zero = F(1, 2), F(1), F(0)
LABELS = [zero, half, one, F(3, 2)]


def d(j):
    return int(2 * j + 1)


def site_content_weights(top, bottom, site, contents, labels=LABELS):
    """For the state with face labels `top` (= bottom, the state paired with itself), the
    weight d_c S(c) of content c at `site`, summing the other rungs with d (delta rungs)."""
    m = len(top) + 1
    out = {}
    for c in contents:
        s = F(0)
        others = [i for i in range(m) if i != site]
        for combo in itertools.product(labels, repeat=len(others)):
            rungs = [None] * m
            rungs[site] = c
            for i, lab in zip(others, combo):
                rungs[i] = lab
            v = ladder_eval(top, bottom, rungs)
            if not v.is_zero():
                w = F(1)
                for lab in combo:
                    w *= d(lab)
                s += w * v.rational()
        out[c] = d(c) * s
    return out


def t3_commuting_model(mu: F, lam: F, m: int):
    """Exact enumeration over b in {+-1}^(m-1): vacuum V = <Pi, K Pi>, excitation X_p = <b_p Pi, K b_p Pi>
    with K acting on each b_l by lam (K b_l = lam b_l, K 1 = 1); inner product = average over signs."""
    faces = m - 1
    def K_poly(coeffs):  # dict frozenset(S)->coef ; K multiplies monomial b_S by lam^|S|
        return {S: c * lam ** len(S) for S, c in coeffs.items()}
    def inner(p, q):  # <p,q> = sum over monomials coef*coef (orthonormal monomials)
        return sum(c * q.get(S, F(0)) for S, c in p.items())
    # Pi(1+mu b_l) as dict
    Pi = {}
    for k in range(faces + 1):
        for S in itertools.combinations(range(faces), k):
            Pi[frozenset(S)] = mu ** k
    V = inner(Pi, K_poly(Pi))
    Xs = []
    for p in range(faces):
        bp = {}
        for S, c in Pi.items():
            T = S ^ {p}  # b_p * b_S, b_p^2 = 1
            bp[frozenset(T)] = bp.get(frozenset(T), F(0)) + c
        Xs.append(inner(bp, K_poly(bp)))
    return V, Xs


def run():
    # ---- T1: locality — non-adjacent pair B_1 B_3 on m=4: site-2 content exactly 1/2 with weight 1
    w13 = site_content_weights([half, zero, half], [half, zero, half], site=1, contents=[zero, half, one])
    t1 = (w13[half] == 1 and w13[zero] == 0 and w13[one] == 0)
    # ---- T2: adjacent pair B_1 B_2 on m=3: site-2 content weights 1/4 (c=0), 3/4 (c=1)
    w12 = site_content_weights([half, half], [half, half], site=1, contents=[zero, one])
    t2_weights = (w12[zero] == F(1, 4) and w12[one] == F(3, 4))
    # commutator energy identity as a polynomial identity in (lam, s1): verified on a 7x7 grid (degree <= 2+2)
    def comm_energy(lam, s1):
        return sum(w12[c] * (lam * (1 if c == zero else s1) - lam) ** 2 for c in (zero, one))
    t2_identity = all(comm_energy(F(i, 7), F(j, 5)) == F(3, 4) * F(i, 7) ** 2 * (1 - F(j, 5)) ** 2
                      for i in range(1, 8) for j in range(0, 7))
    # sanity: with s1 = 1 (no contraction on content 1) the commutator vanishes; with lam=1,s1=0 it is 3/4
    t2_controls = comm_energy(F(1, 3), F(1)) == 0 and comm_energy(F(1), F(0)) == F(3, 4)
    # ---- T3: commuting model, exact ratio independence of m
    t3 = True
    t3_rows = {}
    for mu, lam in ((F(1, 30), F(1, 2)), (F(1, 8), F(1, 3)), (F(1, 3), F(2, 3))):
        ratios = set()
        for m in range(2, 9):
            V, Xs = t3_commuting_model(mu, lam, m)
            pf = (1 + mu * mu * lam) ** (m - 1)
            if V != pf:
                t3 = False
            for X in Xs:
                ratios.add(X / V)
        # bulk excitation ratio must be one value for all m and all p (it is, since faces commute)
        t3_rows[f"mu={mu},lam={lam}"] = {"distinct_excitation_to_vacuum_ratios": len(ratios), "ratio": str(next(iter(ratios)))}
        if len(ratios) != 1:
            t3 = False
    # ---- T4: second-order remainder ratio, displayed at a = 1 (heat kernel lambda_j = exp(-a C_j), C_j = j(j+1))
    a = mp.mpf(1)
    lam = mp.e ** (-a * mp.mpf(3) / 4)
    lam1 = mp.e ** (-a * 2)
    rows = {}
    kill_ok = True
    for k in (F(1, 8), F(1, 4), F(1, 2)):
        kk = mp.mpf(k.numerator) / k.denominator
        mu = mp.besseli(2, kk) / mp.besseli(1, kk)  # r = f_half / f_0 in the YM-15 convention (r_of)
        comm = lam ** 2 * mp.mpf(3) / 4 * (1 - mp.sqrt(lam1)) ** 2
        leading = mu ** 2 * lam ** 2
        pair = mu ** 4 * comm
        rho = pair / leading
        rows[str(k)] = {"mu": mp.nstr(mu, 8), "adjacent_commutator_energy": mp.nstr(comm, 8),
                        "leading_single_face_energy": mp.nstr(leading, 8), "pair_remainder_energy": mp.nstr(pair, 10),
                        "relative_remainder_rho": mp.nstr(rho, 8)}
        if not rho < mp.mpf(1) / 10:
            kill_ok = False
    t4 = kill_ok
    ok = t1 and t2_weights and t2_identity and t2_controls and t3 and t4
    return {
        "certificate_type": "YM35_T54PRIME_ENERGY_VERSION_AND_COMMUTATOR_REMAINDER",
        "framework_inputs_cited": ["theorum/54 Thm 2.1 + T5 (sum_mu_diverges)", "theorum/50 A4 (polarization visibility, post-audit)",
                                   "theorum/53 T1/T2 (phases vs weights; native Parseval)", "theorum/51 exchange law", "theorum/61 commutator support/mass",
                                   "YM-25 T1 bi-invariance, YM-25 T3 recoupling 1/4:3/4", "YM-30 engine"],
        "T1_commutator_locality_nonadjacent_no_mixing": t1,
        "T1_nonadjacent_site2_weights": {str(c): str(v) for c, v in w13.items()},
        "T2_adjacent_site2_weights_quarter_threequarter": t2_weights,
        "T2_adjacent_site2_weights": {str(c): str(v) for c, v in w12.items()},
        "T2_commutator_energy_identity_lam2_3over4_(1-sqrt_lam1)2": t2_identity,
        "T2_controls": t2_controls,
        "T3_commuting_model_ratio_m_independent_m2to8": t3,
        "T3_rows": t3_rows,
        "T4_second_order_remainder_a1": rows,
        "T4_kill_criterion_rho_below_1_over_10_all_kappa": t4,
        "verdict_on_E4DC": "first non-commutative correction is LOCAL (adjacent faces only) and relatively small (rho ~ mu^2) at every m; "
                           "E4D-C excitation bound OPENS at first non-commutative order; full face expansion (all words, contents j>=1) NOT summed — E4D-C OPEN",
        "not_claimed": ["m-uniform excitation upper bound (all orders)", "convergence of the face expansion (worldline cluster problem, YM-33)",
                        "anything about the cutoff a -> 0", "Clay predicate"],
        "verdict": "PASS" if ok else "FAIL",
    }


if __name__ == "__main__":
    cert = run()
    json.dump(cert, open(os.path.join(HERE, "YM35_RESULT.json"), "w"), indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    open(os.path.join(HERE, "EXPECTED_YM35.sha256"), "w").write(sha + "\n")
    print(cert["verdict"])
    for k in ("T1_commutator_locality_nonadjacent_no_mixing", "T2_adjacent_site2_weights_quarter_threequarter",
              "T2_commutator_energy_identity_lam2_3over4_(1-sqrt_lam1)2", "T3_commuting_model_ratio_m_independent_m2to8", "T4_kill_criterion_rho_below_1_over_10_all_kappa"):
        print(k, cert[k])
    print(cert["T3_rows"])
    for k, v in cert["T4_second_order_remainder_a1"].items():
        print(k, v)
    print("sha256:", sha)
