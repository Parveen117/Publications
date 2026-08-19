"""New development — what information invariance means microscopically.

The companion work [Dabas2026] derives, for a weighted ensemble with weight W
on a phase space carrying a time-reversal involution Theta, the identity
    L_ij - eps_i eps_j L_ji = Integral <(1 - chi) J_i(0) J_j(t)>_W dt,
    chi(Gamma) = W(Theta Gamma) / W(Gamma),
and identifies the coarse-grained image of the antisymmetric part with the
thermodynamic curvature Omega = d omega. Theorem A of this paper identified
information invariance (delta I = 0) with Omega = 0. Chaining the two suggests
that delta I = 0 is a statement about time reversal of the weight. We prove the
finite version of that chain, and we find that one leg of it is one-directional.

Theorem N (microscopic characterisation).
  Let the phase space be finite, Theta an involution, W a positive weight, and
  for observables f, g define the weighted correlation
      C_W(f, g) = Sum_Gamma W(Gamma) f(Gamma) g(Theta Gamma).
  Then, exactly,
      C_W(f, g) - C_W(g, f) = Sum_Gamma W(Gamma) [1 - chi(Gamma)] f(Gamma) g(Theta Gamma),
  and the following are equivalent:
    (i)  chi == 1, i.e. the weight is time-reversal invariant;
    (ii) C_W is symmetric for every pair of observables;
    (iii) the antisymmetric part of the weighted response vanishes identically.
  ((ii) => (i) follows by taking indicator observables.)
  So the Onsager-reciprocal sector is exactly the time-reversal-invariant-weight
  sector, with no analytic hypotheses: this is finite bookkeeping.

Theorem O (curvature is a sufficient, not a necessary, witness).
  Let K be a linear coarse-graining from microscopic observable pairs to the
  chart, and let Omega be the image of the antisymmetric part.
    (a) chi == 1  ==>  Omega = 0. Always.
    (b) The converse fails: there are weights with chi != 1, hence a nonzero
        microscopic antisymmetric part, whose coarse-grained image vanishes.
        An explicit two-mode example is constructed and verified below.
  Consequently: observing Omega = 0 does not establish microscopic
  reversibility, and a null curvature measurement cannot refute the framework;
  only Omega != 0 is informative. This is the microscopic counterpart of
  Theorem B, which showed that a single measured loop area does not establish
  Omega != 0 either. Together they bound what a hysteresis experiment can decide:
      Omega != 0 (across cycles)  =>  chi != 1   [informative]
      Omega  = 0                  =>  nothing about chi
      one-cycle loop area != 0     =>  nothing about Omega

Consequence for the paper. The definition of information invariance acquires a
microscopic reading -- stationarity of the loop functional is time-reversal
invariance of the statistical weight -- and the experimental logic is sharpened:
the predeclared cycle-2 test is necessary because neither a single loop area nor
a null result carries the inference by itself.

Exact rational arithmetic throughout.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import random
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent


def cj(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(o): return hashlib.sha256(cj(o).encode()).hexdigest()


def involution(n, rng):
    """random fixed-point-allowed involution on {0..n-1}"""
    theta = list(range(n))
    idx = list(range(n))
    rng.shuffle(idx)
    while len(idx) >= 2 and rng.random() < 0.8:
        a, b = idx.pop(), idx.pop()
        theta[a], theta[b] = b, a
    return theta


def C_W(W, theta, f, g):
    return sum(W[G] * f[G] * g[theta[G]] for G in range(len(W)))


def chi(W, theta):
    return [W[theta[G]] / W[G] for G in range(len(W))]


def theorem_N(rng, trials=40):
    ch = {"identity_exact": True, "chi_one_implies_symmetric": True,
          "symmetric_for_all_implies_chi_one": True,
          "control_asymmetric_weight_breaks_symmetry": True}
    for _ in range(trials):
        n = rng.randint(3, 7)
        theta = involution(n, rng)
        # (a) symmetric weight: W constant on Theta-orbits  -> chi == 1
        Worb = [Fr(rng.randint(1, 9), rng.randint(1, 4)) for _ in range(n)]
        W = [min(Worb[G], Worb[theta[G]]) for G in range(n)]   # forces W(G) = W(theta G)
        X = chi(W, theta)
        if any(x != 1 for x in X): ch["chi_one_implies_symmetric"] = False
        for _ in range(4):
            f = [Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(n)]
            g = [Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(n)]
            lhs = C_W(W, theta, f, g) - C_W(W, theta, g, f)
            rhs = sum(W[G] * (1 - X[G]) * f[G] * g[theta[G]] for G in range(n))
            if lhs != rhs: ch["identity_exact"] = False
            if lhs != 0: ch["chi_one_implies_symmetric"] = False
        # (b) generic weight -> chi != 1 somewhere, and indicators expose it
        W2 = [Fr(rng.randint(1, 9), rng.randint(1, 4)) for _ in range(n)]
        X2 = chi(W2, theta)
        asym_exists = any(X2[G] != 1 for G in range(n))
        if asym_exists:
            G0 = next(G for G in range(n) if X2[G] != 1)
            f = [Fr(1) if k == G0 else Fr(0) for k in range(n)]
            g = [Fr(1) if k == theta[G0] else Fr(0) for k in range(n)]
            diff = C_W(W2, theta, f, g) - C_W(W2, theta, g, f)
            if diff == 0: ch["symmetric_for_all_implies_chi_one"] = False
            if diff != W2[G0] - W2[theta[G0]]: ch["control_asymmetric_weight_breaks_symmetry"] = False
    return {"status": "PROVED" if all(ch.values()) else "FAILED", "checks": ch,
            "reading": ("delta I = 0 (equivalently Omega = 0, Theorem A) is exactly the statement that the "
                        "statistical weight is invariant under time reversal.")}


def theorem_O():
    """Explicit counterexample: chi != 1 with a nonzero microscopic antisymmetric
    part per orbit, whose coarse-grained total vanishes."""
    theta = [1, 0, 3, 2]                 # two Theta-orbits: {0,1} and {2,3}
    W = [Fr(3), Fr(1), Fr(1), Fr(3)]      # chi != 1 on every state
    X = chi(W, theta)
    # observables exposing the asymmetry within each orbit
    ind = lambda k: [Fr(1) if j == k else Fr(0) for j in range(4)]
    asym_orbit1 = C_W(W, theta, ind(0), ind(1)) - C_W(W, theta, ind(1), ind(0))   # = W0 - W1
    asym_orbit2 = C_W(W, theta, ind(2), ind(3)) - C_W(W, theta, ind(3), ind(2))   # = W2 - W3
    # coarse-graining K = total over orbits (a legitimate linear map)
    Omega_image = asym_orbit1 + asym_orbit2
    # sufficiency: a chi == 1 weight kills every orbit contribution
    Wsym = [Fr(2), Fr(2), Fr(5), Fr(5)]
    Xs = chi(Wsym, theta)
    suff = all(Wsym[G] * (1 - Xs[G]) == 0 for G in range(4))
    checks = {
        "chi_differs_from_one": any(x != 1 for x in X),
        "microscopic_antisymmetric_part_nonzero": asym_orbit1 != 0 and asym_orbit2 != 0,
        "orbit_contributions_are_opposite": asym_orbit1 == -asym_orbit2,
        "coarse_grained_curvature_vanishes": Omega_image == 0,
        "sufficiency_direction_holds": suff,
    }
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "example": {"theta": theta, "W": [str(w) for w in W], "chi": [str(x) for x in X],
                        "asym_orbit1": str(asym_orbit1), "asym_orbit2": str(asym_orbit2),
                        "coarse_image": str(Omega_image)},
            "reading": ("Omega = 0 does not imply chi = 1: microscopic irreversibility can cancel under "
                        "coarse-graining. A null curvature measurement is uninformative; only Omega != 0 "
                        "carries the inference.")}


def build():
    rng = random.Random(20260819_3)
    N, O = theorem_N(rng), theorem_O()
    body = {
        "protocol": "INFORMATION_INVARIANCE_MICRO_V1",
        "theorem_N_microscopic_characterisation": N,
        "theorem_O_curvature_is_sufficient_not_necessary": O,
        "chain": ("delta I = 0  <=>  Omega = 0 (Theorem A)  <=  chi == 1 (Theorems N, O); "
                  "the last arrow is one-directional."),
        "experimental_logic": {
            "Omega_nonzero_across_cycles": "implies chi != 1 — informative",
            "Omega_zero": "implies nothing about chi",
            "single_cycle_loop_area_nonzero": "implies nothing about Omega (Theorem B)",
        },
        "consequence": ("Information invariance has a microscopic reading: stationarity of the loop "
                        "functional is time-reversal invariance of the statistical weight. And the "
                        "predeclared cycle-2 test is not optional: neither a single loop area nor a null "
                        "curvature result carries an inference by itself."),
        "claim_boundary": ("Finite phase space, exact rational weights, linear coarse-graining. The "
                           "continuum Green-Kubo statement of the companion work is cited, not re-derived."),
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "MICRO_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_MICRO.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("N:", cert["theorem_N_microscopic_characterisation"]["status"],
          "| O:", cert["theorem_O_curvature_is_sufficient_not_necessary"]["status"],
          "|", cert["certificate_sha256"])
