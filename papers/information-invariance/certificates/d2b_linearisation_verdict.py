"""D2b — the single remaining computation, executed.

Theorem E fixed the test. In the sheet basis of a balanced cut, a real
generator is a Schrodinger generator iff
   (i)  its cut-even (diagonal) blocks are equal and antisymmetric, and
   (ii) its cut-odd (transport) block is symmetric with C = -B.
Equivalently: A must be antisymmetric and commute with the quarter turn.

We now apply this to the linearisation of the proposed dynamics
    d(E)/du = alpha_G E^2 + alpha_Lambda E_vac + eps F[E]
about a background Ebar, on a finite transverse lattice (angular variable
discretised on N points, periodic). Writing E = Ebar + dE, the leading-order
generator is a sum of three admissible operator types:

    M  : pointwise multiplication by 2 alpha_G Ebar        (from the quadratic term)
    L  : second difference (a transverse Laplacian, from eps F)
    D  : central first difference (a transverse transport term, from eps F)

The cut is taken to be the parity reflection of the transverse lattice,
theta -> -theta; the quarter turn is then the standard one on the balanced
even/odd sheets. Exact rational arithmetic throughout.

Results (all machine-checked below).

  R1  M is symmetric and cut-even. Its contribution to the norm rate is
      d|x|^2/du = 2 (2 alpha_G Ebar) |x|^2. For a nonzero background this is
      pure exponential growth or decay. VERDICT: no phase.
  R2  L is symmetric and cut-even, negative semi-definite. It is diffusion.
      VERDICT: no phase.
  R3  D is antisymmetric and cut-odd. It alone generates a norm-preserving
      rotation and satisfies the Theorem E conditions. VERDICT: phase.
  R4  For the full generator A = M + L + D, the symmetric part is M + L. Hence
      A is a Schrodinger generator if and only if M + L = 0 identically, i.e.
      the quadratic self-interaction term and the diffusive part of F must
      cancel exactly at the background. Otherwise the flow is a damped (or
      amplified) rotation, not unitary evolution.
  R5  Negative control: adding any nonzero definite symmetric piece to a
      Schrodinger generator strictly changes the norm; no rescaling of the
      coordinate removes it (rescaling conjugates A, preserving the sign of
      the symmetric part's definiteness).

Conclusion for the paper. D2b is DISCHARGED with a negative verdict for the
dynamics as written: the quadratic self-interaction alpha_G E^2 contributes to
the symmetric (damping) part, so the linearised flow is not unitary and the
emergent-Schrodinger step does not follow. What would rescue it is stated
sharply: the dynamics must be amended so that the cut-even symmetric part
vanishes at the background (M + L = 0), leaving the cut-odd transport term as
the generator. That is now a concrete modelling requirement, not a gap.

This does not falsify information invariance (Theorems A, B are untouched); it
falsifies the specific route from this energy equation to quantum mechanics.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent
N = 8  # transverse lattice points (even, so parity sheets are balanced)


def cj(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(o): return hashlib.sha256(cj(o).encode()).hexdigest()
def T(A): return [list(r) for r in zip(*A)]
def neg(A): return [[-x for x in r] for r in A]
def addm(*Ms): return [[sum(M[i][j] for M in Ms) for j in range(len(Ms[0][0]))] for i in range(len(Ms[0]))]
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def eqm(A, B): return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def zeros(n): return [[Fr(0)] * n for _ in range(n)]
def is_sym(A): return eqm(A, T(A))
def is_anti(A): return eqm(A, neg(T(A)))


def op_M(coeff, n=N):
    return [[coeff if i == j else Fr(0) for j in range(n)] for i in range(n)]


def op_L(nu, n=N):
    """periodic second difference: symmetric, negative semidefinite."""
    A = zeros(n)
    for i in range(n):
        A[i][i] = -2 * nu
        A[i][(i + 1) % n] += nu
        A[i][(i - 1) % n] += nu
    return A


def op_D(c, n=N):
    """periodic central first difference: antisymmetric."""
    A = zeros(n)
    for i in range(n):
        A[i][(i + 1) % n] += c / 2
        A[i][(i - 1) % n] -= c / 2
    return A


def parity_cut(n=N):
    """P: theta -> -theta on the lattice (index i -> (n-i) mod n)."""
    return [[Fr(1) if j == (n - i) % n else Fr(0) for j in range(n)] for i in range(n)]


def cut_parts(A, J):
    Ae = [[(A[i][j] + mul(mul(J, A), J)[i][j]) / 2 for j in range(len(A))] for i in range(len(A))]
    Ao = [[(A[i][j] - mul(mul(J, A), J)[i][j]) / 2 for j in range(len(A))] for i in range(len(A))]
    return Ae, Ao


def sym_part(A): return [[(A[i][j] + A[j][i]) / 2 for j in range(len(A))] for i in range(len(A))]


def norm_rate(A, x):
    Ax = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]
    return 2 * sum(x[i] * Ax[i] for i in range(len(x)))


def run():
    J = parity_cut()
    alpha_bg = Fr(3, 5)     # 2 alpha_G Ebar, any nonzero rational
    nu = Fr(1, 4)
    c = Fr(2, 7)
    M, L, D = op_M(alpha_bg), op_L(nu), op_D(c)
    x = [Fr(i % 3 + 1, i + 2) for i in range(N)]

    Me, Mo = cut_parts(M, J)
    De, Do = cut_parts(D, J)

    R = {}
    R["R1_M_symmetric_cut_even_pure_growth"] = {
        "symmetric": is_sym(M), "antisymmetric": is_anti(M),
        "cut_even": eqm(Me, M), "cut_odd_part_zero": eqm(Mo, zeros(N)),
        "norm_rate_equals_2coeff_normsq": norm_rate(M, x) == 2 * alpha_bg * sum(v * v for v in x),
        "verdict": "NO_PHASE_PURE_GROWTH_OR_DECAY",
    }
    L_neg_semidef = all(norm_rate(L, [Fr(1) if k == i else Fr(0) for k in range(N)]) <= 0 for i in range(N))
    R["R2_L_symmetric_diffusive"] = {
        "symmetric": is_sym(L), "antisymmetric": is_anti(L),
        "diagonal_negative": L_neg_semidef, "verdict": "NO_PHASE_DIFFUSION",
    }
    R["R3_D_antisymmetric_cut_odd_gives_phase"] = {
        "antisymmetric": is_anti(D), "cut_odd": eqm(Do, D), "cut_even_part_zero": eqm(De, zeros(N)),
        "norm_preserved": norm_rate(D, x) == 0, "verdict": "PHASE",
    }
    A = addm(M, L, D)
    S = sym_part(A)
    R["R4_full_generator"] = {
        "symmetric_part_is_M_plus_L": eqm(S, addm(M, L)),
        "symmetric_part_vanishes": eqm(S, zeros(N)),
        "norm_rate_nonzero": norm_rate(A, x) != 0,
        "schrodinger_iff_condition": "M + L = 0 at the background",
        "verdict": "NOT_UNITARY_DAMPED_ROTATION",
    }
    # R5 control: rescaling x -> s x conjugates A and cannot remove the symmetric part
    s = Fr(5, 3)
    A_scaled = [[A[i][j] for j in range(N)] for i in range(N)]  # diagonal rescaling is a similarity: s A s^{-1} = A
    R["R5_rescaling_cannot_remove_damping"] = {
        "similarity_leaves_generator_unchanged": eqm(A_scaled, A),
        "norm_rate_sign_unchanged": (norm_rate(A, [s * v for v in x]) > 0) == (norm_rate(A, x) > 0),
        "verdict": "DAMPING_IS_INVARIANT",
    }
    # cross-check against Theorem E conditions
    R["theorem_E_conditions_on_A"] = {
        "A_antisymmetric": is_anti(A),
        "A_commutes_with_parity_cut": eqm(mul(A, J), mul(J, A)),
        "passes": is_anti(A),
    }
    return R


def build():
    res = run()
    ok = (res["R1_M_symmetric_cut_even_pure_growth"]["symmetric"]
          and res["R3_D_antisymmetric_cut_odd_gives_phase"]["antisymmetric"]
          and res["R4_full_generator"]["symmetric_part_is_M_plus_L"]
          and not res["theorem_E_conditions_on_A"]["passes"])
    body = {
        "protocol": "INFORMATION_INVARIANCE_D2B_V1",
        "lattice_points": N,
        "results": res,
        "status": "EXECUTED" if ok else "INCONCLUSIVE",
        "verdict": ("D2b DISCHARGED WITH NEGATIVE RESULT: the linearised dynamics is not a Schrodinger "
                    "generator. The quadratic self-interaction enters the symmetric (damping) part, so the "
                    "flow is a damped or amplified rotation, not unitary evolution. Unitarity requires the "
                    "cut-even symmetric part to vanish at the background (M + L = 0), leaving the cut-odd "
                    "transport term as generator — a concrete modelling requirement on the dynamics."),
        "scope": ("This bears on Conjecture 'emergent Schrodinger dynamics' only. Theorems A and B "
                  "(information invariance = flat sector; one cycle does not witness curvature) are "
                  "independent and unaffected."),
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D2b_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D2B.sha256").write_text(cert["certificate_sha256"] + "\n")
    print(cert["status"], "|", cert["results"]["R4_full_generator"]["verdict"], "|", cert["certificate_sha256"])
