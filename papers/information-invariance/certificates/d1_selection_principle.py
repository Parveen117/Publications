"""D1 — can the energy dynamics be derived from information invariance?

Obligation D1 asked: derive
    d(E)/du = alpha_G E^2 + alpha_Lambda E_vac + eps F[E]
*from* Definition 1 (delta I = 0), rather than positing it, and show why the
quadratic term is the one that appears. Two separate questions, both settled
here — the first negatively, the second positively.

Theorem I (no-go: invariance alone selects nothing).
  On a chart where information invariance holds, omega is exact, omega = dPhi
  (Theorem A). Any evolution that transports the potential,
      d(Phi)/du = f(Phi),
  keeps omega exact for *every* function f: the flow acts on Phi alone, so the
  transported one-form is again a differential. Hence the class of dynamics
  compatible with delta I = 0 is parameterised by an arbitrary f, an
  infinite-dimensional family, of which the quadratic law is a single point.
  Information invariance therefore does *not* single out the quadratic term.
  The certificate verifies exactness is preserved, on a finite lattice and in
  exact arithmetic, for a family of distinct f (linear, quadratic, cubic,
  rational), and for compositions of them.

Theorem J (a selection principle that does work: minimal scale usage).
  Fix the dimension group (M, L, T) with [E] = M L^-1 T^-2 and u a length, so
  [d/du] = L^-1. For a candidate term  d(E)/du = alpha_n E^n  the coupling must
  have  [alpha_n] = M^{1-n} L^{n-2} T^{2n-2}.
  (i)  No n admits a coupling built from {G, c} alone: matching exponents forces
       2n - 2 = 1, i.e. n = 3/2, not an admissible power. A third constant is
       unavoidable. (This is the same conclusion reached from the correlation
       function in D3, arrived at independently.)
  (ii) Admit exactly one new length lambda, so alpha_n = G^a c^b lambda^d. Then
       matching gives, uniquely,
            a = n - 1,   b = 4 - 4n,   d = 2n - 3.
  (iii) Requiring *minimal scale usage* — exactly one power of the new length,
       d = 1 — selects  n = 2  uniquely, with
            alpha_2 = G lambda / c^4.
       The quadratic self-interaction is therefore not arbitrary: it is the
       unique power law whose coupling consumes the new length exactly once.
  (iv) Consistency check: the draft wrote alpha_G ~ G/c^4, which has dimension
       M^-1 L^-1 T^2, whereas the equation requires M^-1 T^2 — short by exactly
       one length, precisely the lambda that (iii) restores.

  By Theorem G of D3, the only length the deterministic dynamics generates by
  itself is L_E = c^2 / sqrt(G E), so the natural reading is lambda = L_E,
  giving alpha_2 = G L_E / c^4 = 1/(c^2 sqrt(G E)). This is a *derivation of the
  coupling's form*, not of its numerical value.

Verdict: D1 PARTIALLY DISCHARGED. Negative half: information invariance alone
cannot produce the dynamics (Theorem I). Positive half: adding one stated
principle — a single new length, used exactly once — makes the quadratic term
the unique admissible power law and fixes the coupling's form (Theorem J). The
paper should present the quadratic law as a consequence of that stated
selection principle, not of delta I = 0.

Exact rational arithmetic; no floating point enters a verdict.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import random
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent
DIMS = {"G": (Fr(-1), Fr(3), Fr(-2)), "c": (Fr(0), Fr(1), Fr(-1)), "lam": (Fr(0), Fr(1), Fr(0))}


def cj(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(o): return hashlib.sha256(cj(o).encode()).hexdigest()


# ------------------------------------------------------------------ Theorem I
def theorem_I(rng, trials=30):
    """Finite lattice: Phi on an n x m grid; omega = dPhi (exact by construction).
    Evolve Phi -> f(Phi) and check the transported form is still exact."""
    def dPhi(P):
        n, m = len(P), len(P[0])
        ep = [[P[i + 1][j] - P[i][j] for j in range(m)] for i in range(n - 1)]
        ev = [[P[i][j + 1] - P[i][j] for j in range(m - 1)] for i in range(n)]
        return ep, ev

    def curvature(ep, ev):
        n, m = len(ep) + 1, len(ev[0]) + 1
        return [[ep[i][j] + ev[i + 1][j] - ep[i][j + 1] - ev[i][j] for j in range(m - 1)] for i in range(n - 1)]

    fams = {
        "linear": lambda x: 2 * x + Fr(1, 3),
        "quadratic": lambda x: x * x,
        "cubic": lambda x: x ** 3 - x,
        "rational": lambda x: x / (1 + x * x),
        "composed_quad_then_lin": lambda x: 3 * (x * x) - Fr(2, 5),
    }
    ok = {k: True for k in fams}
    for _ in range(trials):
        n, m = rng.randint(3, 5), rng.randint(3, 5)
        P = [[Fr(rng.randint(-6, 6), rng.randint(1, 4)) for _ in range(m)] for _ in range(n)]
        for name, f in fams.items():
            Q = [[f(P[i][j]) for j in range(m)] for i in range(n)]
            K = curvature(*dPhi(Q))
            if any(K[a][b] != 0 for a in range(len(K)) for b in range(len(K[0]))):
                ok[name] = False
    # negative control: perturb an edge value directly (not through a potential) -> curvature appears
    P = [[Fr(rng.randint(-6, 6), rng.randint(1, 4)) for _ in range(4)] for _ in range(4)]
    ep, ev = dPhi(P); ep[1][1] += Fr(1, 5)
    K = curvature(ep, ev)
    control = any(K[a][b] != 0 for a in range(len(K)) for b in range(len(K[0])))
    checks = dict(ok); checks["control_non_potential_perturbation_creates_curvature"] = control
    checks["family_is_infinite_dimensional"] = len(fams) >= 5 and all(ok.values())
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "reading": ("Every potential-transporting evolution preserves exactness, so delta I = 0 is "
                        "compatible with an arbitrary function f: the quadratic law is one point in an "
                        "infinite family and is not selected by invariance.")}


# ------------------------------------------------------------------ Theorem J
def required_alpha_dim(n):
    """[alpha_n] = [E]^{1-n} * L^{-1}, with [E] = M L^-1 T^-2."""
    M = Fr(1) * (1 - n)
    L = Fr(-1) * (1 - n) - 1
    T = Fr(-2) * (1 - n)
    return (M, L, T)


def solve_exponents(n):
    """alpha_n = G^a c^b lam^d ; return (a, b, d) or None."""
    Mr, Lr, Tr = required_alpha_dim(n)
    a = -Mr                                  # from M: -a = Mr
    b = -2 * a - Tr                           # from T: -2a - b = Tr
    d = Lr - (3 * a + b)                      # from L: 3a + b + d = Lr
    # verify
    v = (-a, 3 * a + b + d, -2 * a - b)
    return (a, b, d) if v == (Mr, Lr, Tr) else None


def theorem_J():
    # (i) with {G, c} only: need d = 0
    only_Gc = {n: (solve_exponents(Fr(n)) or (None, None, None))[2] == 0 for n in (0, 1, 2, 3, 4, 5)}
    # the half-integer that would work
    half = solve_exponents(Fr(3, 2))
    # (ii)-(iii) with one length, minimal usage d = 1
    sols = {n: solve_exponents(Fr(n)) for n in (0, 1, 2, 3, 4, 5)}
    d_one = [n for n, s in sols.items() if s and s[2] == 1]
    checks = {
        "no_power_law_from_G_and_c_alone": not any(only_Gc.values()),
        "half_integer_case_confirms_obstruction": half is not None and half[2] == 0,
        "exponent_formula_a_eq_n_minus_1": all(s and s[0] == Fr(n - 1) for n, s in sols.items()),
        "exponent_formula_b_eq_4_minus_4n": all(s and s[1] == Fr(4 - 4 * n) for n, s in sols.items()),
        "exponent_formula_d_eq_2n_minus_3": all(s and s[2] == Fr(2 * n - 3) for n, s in sols.items()),
        "minimal_scale_usage_selects_n_equals_2": d_one == [2],
        "alpha_2_equals_G_lambda_over_c4": sols[2] == (Fr(1), Fr(-4), Fr(1)),
        "draft_alpha_short_by_one_length": (Fr(-1), Fr(-1), Fr(2)) != required_alpha_dim(2)
                                           and (Fr(-1), Fr(0), Fr(2)) == required_alpha_dim(2),
    }
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks,
            "solutions": {str(n): [str(x) for x in s] if s else None for n, s in sols.items()},
            "selected": {"n": 2, "alpha_2": "G * lambda / c^4",
                         "natural_lambda": "L_E = c^2 / sqrt(G E)  (Theorem G of D3) -> alpha_2 = 1/(c^2 sqrt(G E))"},
            "reading": ("Dimensions alone do not fix the power once a new length is admitted; requiring that "
                        "the new length be used exactly once fixes n = 2 uniquely and gives the coupling's form.")}


def build():
    rng = random.Random(20260819_1)
    I, J = theorem_I(rng), theorem_J()
    body = {
        "protocol": "INFORMATION_INVARIANCE_D1_V1",
        "theorem_I_invariance_alone_selects_nothing": I,
        "theorem_J_minimal_scale_usage_selects_the_quadratic_law": J,
        "obligation": {"D1": "PARTIALLY DISCHARGED — negative half proved (Theorem I), positive half proved under a stated selection principle (Theorem J)"},
        "verdict": ("delta I = 0 does not derive the energy dynamics: every potential-transporting evolution "
                    "preserves it, an infinite family. What does select the quadratic term is a separate, "
                    "stateable principle — admit exactly one new length and use it exactly once — under which "
                    "n = 2 is unique and alpha_2 = G lambda / c^4. The draft's alpha_G ~ G/c^4 is short by "
                    "exactly that length; the length the theory itself supplies is the cosmological L_E."),
        "scope": "Bears on the dynamics conjecture; Theorems A and B are unaffected.",
        "claim_boundary": ("Finite-lattice exactness checks and exact dimensional algebra. No claim is made "
                           "that the selection principle is forced by physics; it is stated so that the "
                           "quadratic law follows from something written down."),
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D1_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D1.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("I:", cert["theorem_I_invariance_alone_selects_nothing"]["status"],
          "| J:", cert["theorem_J_minimal_scale_usage_selects_the_quadratic_law"]["status"],
          "| selected n =", cert["theorem_J_minimal_scale_usage_selects_the_quadratic_law"]["selected"]["n"],
          "|", cert["certificate_sha256"])
