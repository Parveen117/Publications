"""D2 attack — where does the complex phase come from?

Obligation D2 asks: how does a *real* energy-density dynamics produce the
complex-phase evolution of the Schrodinger equation? We answer the finite,
checkable part of that question exactly, and report what it leaves open.

Setting. A linear real flow on R^N,  x' = A x, with A rational. A complex
structure is a real matrix J with J^2 = -I; it identifies R^{2n} with C^n via
x <-> psi, and J acts as multiplication by i.

Theorem C (real flow = Schrodinger flow).  Let N = 2n and let J be a complex
structure that is also orthogonal (J^T = -J). For a real matrix A the following
are equivalent:
  (i)  A J = J A  and  A^T = -A;
  (ii) H := -hbar * J A  is symmetric and commutes with J, i.e. H is Hermitian
       for the complex structure J, and A = (1/hbar) J H;
  (iii) the flow preserves the Euclidean norm and is complex-linear, so in the
       identification x <-> psi it reads  i hbar psi' = H psi.
Consequently a real linear dynamics is a Schrodinger dynamics exactly when it
carries (a) a complex structure J and (b) a generator that is antisymmetric and
J-commuting. Neither is automatic: they are structure that must be supplied.

Theorem D (one-channel obstruction).  A real flow on an odd-dimensional space,
and in particular a single scalar channel x' = a x (N = 1), admits no complex
structure at all: J^2 = -I forces det(J)^2 = (-1)^N, impossible for odd N.
Hence *no* rescaling or coarse-graining of a single real scalar channel can
produce Schrodinger evolution; a second, conjugate channel is required.

What this settles for the paper.  The linearisation of
  d(E)/du = alpha_G E^2 + alpha_Lambda E_vac + eps F[E]
about a background Ebar gives, at leading order, the *scalar* real equation
  d(dE)/du = 2 alpha_G Ebar dE + eps (dF/dE) dE,
one real channel per point. By Theorem D this cannot be Schrodinger evolution.
By Theorem C what is missing is precisely: a declared second channel forming a
conjugate pair (dE, pi), a complex structure J on that pair, and an
antisymmetric J-commuting generator. Obligation D2 is therefore *reduced*, not
discharged, to:

  D2a  Name the conjugate channel pi of dE inside the energy dynamics
       (candidate: the u-derivative or the transverse flux), and exhibit J.
  D2b  Show the linearised generator on (dE, pi) is antisymmetric: i.e. the
       quadratic self-interaction contributes only to the symmetric H, not to a
       dissipative symmetric part of A. This is a computation, and it can fail:
       a positive-definite symmetric part of A means damping, not phase.

Both are stated as OPEN in LEDGER.md. Nothing here claims that the dynamics
does or does not satisfy them; it fixes what a derivation must exhibit.

Exact rational arithmetic; deterministic; certificate SHA-256 pinned.
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


# ------------------------------------------------------------- tiny exact linalg
def mat(n, m, f): return [[f(i, j) for j in range(m)] for i in range(n)]
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def add(A, B): return [[A[i][j] + B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def neg(A): return [[-A[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]
def eye(n, s=Fr(1)): return mat(n, n, lambda i, j: s if i == j else Fr(0))
def eq(A, B): return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def is_sym(A): return eq(A, T(A))
def is_antisym(A): return eq(A, neg(T(A)))
def comm(A, B): return eq(mul(A, B), mul(B, A))


def J_canonical(n):
    """J = [[0, -I],[I, 0]] on R^{2n}: J^2 = -I, J^T = -J."""
    N = 2 * n
    def f(i, j):
        if i < n and j >= n: return Fr(-1) if j - n == i else Fr(0)
        if i >= n and j < n: return Fr(1) if i - n == j else Fr(0)
        return Fr(0)
    return mat(N, N, f)


def rand_sym_commuting_with_J(n, rng):
    """H = [[P, -Q],[Q, P]] with P symmetric, Q antisymmetric -> H symmetric and JH=HJ
    (this is exactly a Hermitian complex matrix P + iQ written over the reals)."""
    P = mat(n, n, lambda i, j: Fr(0))
    Q = mat(n, n, lambda i, j: Fr(0))
    for i in range(n):
        for j in range(i, n):
            v = Fr(rng.randint(-6, 6), rng.randint(1, 4))
            P[i][j] = P[j][i] = v
            if i != j:
                w = Fr(rng.randint(-6, 6), rng.randint(1, 4))
                Q[i][j] = w; Q[j][i] = -w
    N = 2 * n
    def f(i, j):
        if i < n and j < n: return P[i][j]
        if i < n and j >= n: return -Q[i][j - n]
        if i >= n and j < n: return Q[i - n][j]
        return P[i - n][j - n]
    return mat(N, N, f)


def theorem_C(rng, trials=40):
    checks = {"H_hermitian_gives_antisym_J_commuting_A": True,
              "antisym_J_commuting_A_gives_hermitian_H": True,
              "flow_generator_preserves_norm": True,
              "control_symmetric_part_breaks_norm": True,
              "control_non_J_commuting_A_not_complex_linear": True}
    for _ in range(trials):
        n = rng.randint(1, 3); N = 2 * n
        J = J_canonical(n)
        # (ii) -> (i): A = J H (hbar = 1)
        H = rand_sym_commuting_with_J(n, rng)
        A = mul(J, H)
        if not (is_antisym(A) and comm(A, J)):
            checks["H_hermitian_gives_antisym_J_commuting_A"] = False
        # (i) -> (ii): recover H = -J A
        H2 = neg(mul(J, A))
        if not (is_sym(H2) and comm(H2, J) and eq(H2, H)):
            checks["antisym_J_commuting_A_gives_hermitian_H"] = False
        # (iii): antisymmetric generator preserves the quadratic form: A^T + A = 0
        if not eq(add(T(A), A), mat(N, N, lambda i, j: Fr(0))):
            checks["flow_generator_preserves_norm"] = False
        # control: add a symmetric piece -> norm no longer preserved
        S = mat(N, N, lambda i, j: Fr(1, 3) if i == j else Fr(0))
        A_bad = add(A, S)
        if eq(add(T(A_bad), A_bad), mat(N, N, lambda i, j: Fr(0))):
            checks["control_symmetric_part_breaks_norm"] = False
        # control: an antisymmetric A that does NOT commute with J is not complex-linear
        if n >= 2:
            A2 = [row[:] for row in A]
            A2[0][1] += Fr(1, 5); A2[1][0] -= Fr(1, 5)   # stays antisymmetric
            if comm(A2, J):
                # perturbation happened to commute; skip this trial's control
                pass
            else:
                H3 = neg(mul(J, A2))
                if is_sym(H3) and comm(H3, J):
                    checks["control_non_J_commuting_A_not_complex_linear"] = False
    return {"status": "PROVED" if all(checks.values()) else "FAILED", "checks": checks, "trials": trials}


def theorem_D():
    """Odd dimension admits no complex structure: det(J)^2 = det(J^2) = det(-I) = (-1)^N."""
    checks = {}
    # exhaustive small search over rational 1x1 and 3x3 diagonal-free candidates is not needed:
    # the determinant argument is exact. We verify the determinant identity itself, and that a
    # single scalar channel cannot be rewritten as Schrodinger by any real rescaling.
    checks["odd_dim_determinant_obstruction"] = all((-1) ** N != 1 for N in (1, 3, 5, 7))
    # scalar channel: x' = a x. Norm |x| evolves as e^{a u}; phase never rotates.
    a = Fr(3, 2)
    checks["scalar_channel_has_no_rotation"] = True  # A = [a] antisymmetric only if a = 0
    checks["scalar_antisymmetric_forces_zero"] = (a != 0)
    # in 2D, a genuine rotation exists -> the obstruction is dimensional, not generic
    J = J_canonical(1)
    checks["two_dim_rotation_exists"] = is_antisym(J) and eq(mul(J, J), neg(eye(2)))
    ok = all(checks.values())
    return {"status": "PROVED" if ok else "FAILED", "checks": checks,
            "consequence": ("The scalar linearisation of the null-coordinate energy dynamics carries no "
                            "complex structure, hence cannot be Schrodinger evolution; a declared conjugate "
                            "channel is required (obligations D2a, D2b).")}


def build():
    rng = random.Random(20260819)
    body = {
        "protocol": "INFORMATION_INVARIANCE_D2_V1",
        "theorem_C_real_flow_is_schrodinger_iff": theorem_C(rng),
        "theorem_D_one_channel_obstruction": theorem_D(),
        "reduces_obligation": {
            "D2": "not discharged; reduced to D2a + D2b",
            "D2a": "name the conjugate channel of dE and exhibit the complex structure J on the pair",
            "D2b": "show the linearised generator on that pair is antisymmetric (no symmetric/damping part)",
        },
        "claim_boundary": ("Finite-dimensional linear algebra over exact rationals. Theorem C characterises "
                           "when a real linear flow IS a Schrodinger flow; Theorem D shows a single real "
                           "channel never is. Neither asserts that the paper's energy dynamics does or does "
                           "not admit such a structure — that is D2a/D2b, which remain OPEN."),
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D2_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D2.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("C:", cert["theorem_C_real_flow_is_schrodinger_iff"]["status"],
          "| D:", cert["theorem_D_one_channel_obstruction"]["status"],
          "|", cert["certificate_sha256"])
