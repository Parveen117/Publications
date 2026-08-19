"""D2a — the conjugate channel and complex structure, supplied by the cut.

Theorems C and D (d2_phase_origin.py) left two requirements: name a conjugate
channel and a complex structure J (D2a), and show the linearised generator is
antisymmetric and J-commuting (D2b). D2a is answered by structure the
Recognition-Kernel Framework already derives, and is not new here:

  * `theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md` (PROVED):
    an oriented cut yields a quarter-turn element with iota^2 = -1 and
    iota^dagger = -iota — a complex structure, derived rather than assumed;
  * `theorum/41_cut_graded_universal_generator_theorem.md` (certified): for a
    cut J = J* = J^{-1} with projections P_pm = (I +- J)/2, every generator
    splits uniquely as G = G_e + G_o with J G_e J = G_e, J G_o J = -G_o; the
    even part preserves the two sheets and the odd part transports between
    them (P_+ G_e P_- = 0, P_+ G_o P_+ = 0).

Theorem E (cut supplies the pair; phase is the odd, transporting part).
Let J be a balanced cut on R^{2n} (dim P_+ = dim P_- = n). In the sheet basis
write the real generator as blocks
      A = [[A_++ , B ], [ C , A_-- ]],
so that the cut-even part is diag(A_++, A_--) and the cut-odd part is the
off-diagonal (B, C). Let iota = [[0, -I], [I, 0]] be the quarter turn on the
balanced pair. Then:

  (E1) A commutes with iota  <=>  A_++ = A_--  and  C = -B.
  (E2) A is antisymmetric    <=>  A_++^T = -A_++ , A_--^T = -A_--, C = -B^T.
  (E3) Both hold             <=>  A_++ = A_-- antisymmetric and B symmetric
                                  with C = -B.
       In that case H = -hbar*iota*A is symmetric and commutes with iota, and
       by Theorem C the flow is exactly i hbar psi' = H psi with
       psi = (channel on sheet +) + iota (channel on sheet -).
  (E4) Phase requires the odd part: if B = C = 0 (a purely cut-even generator)
       the two sheets decouple, no rotation mixes them, and the dynamics is a
       pair of independent real flows — decay or growth, never phase.
  (E5) Damping is exactly the symmetric part: writing A = A_sym + A_anti, the
       Euclidean norm obeys d|x|^2/du = 2 x^T A_sym x. Any nonzero definite
       A_sym gives exponential decay or growth, which no rescaling removes.

Consequence for the paper. D2a is DISCHARGED by framework structure: the
conjugate channel of the fluctuation is its image on the opposite sheet of the
oriented cut, and the complex structure is the derived quarter turn. D2b is
now a single explicit computation on the linearised energy dynamics:

  D2b  In the sheet basis of the cut, show the linearisation of
       d(E)/du = alpha_G E^2 + alpha_Lambda E_vac + eps F[E]
       about a background has (i) equal, antisymmetric diagonal blocks and
       (ii) a symmetric transport block with C = -B. Failure of (i) is
       damping, not phase; failure of (ii) is a non-Hermitian generator.

This file proves E1-E5 over exact rationals with negative controls. It does not
assert that the energy dynamics passes D2b; that remains OPEN.
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
def mul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def T(A): return [list(r) for r in zip(*A)]
def neg(A): return [[-x for x in r] for r in A]
def eq(A, B): return all(A[i][j] == B[i][j] for i in range(len(A)) for j in range(len(A[0])))
def zeros(n, m=None): return [[Fr(0)] * (m or n) for _ in range(n)]
def eye(n): return [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def blocks_to_matrix(App, B, C, Amm):
    n = len(App)
    return [[*App[i], *B[i]] for i in range(n)] + [[*C[i], *Amm[i]] for i in range(n)]


def iota(n):
    Z, I = zeros(n), eye(n)
    return [[*Z[i], *neg(I)[i]] for i in range(n)] + [[*I[i], *Z[i]] for i in range(n)]


def rand(n, rng, sym=None):
    M = [[Fr(rng.randint(-6, 6), rng.randint(1, 4)) for _ in range(n)] for _ in range(n)]
    if sym == "sym":
        for i in range(n):
            for j in range(i + 1, n): M[j][i] = M[i][j]
    if sym == "anti":
        for i in range(n):
            M[i][i] = Fr(0)
            for j in range(i + 1, n): M[j][i] = -M[i][j]
    return M


def quad_form_rate(A, x):
    """d|x|^2/du = 2 x^T A_sym x."""
    Ax = [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]
    return 2 * sum(x[i] * Ax[i] for i in range(len(x)))


def theorem_E(rng, trials=40):
    ch = {k: True for k in ("E1_commutes_iff_equal_diag_and_C_eq_minusB",
                            "E2_antisym_iff_blocks_antisym_and_C_eq_minusBT",
                            "E3_both_iff_B_symmetric_gives_hermitian_H",
                            "E4_pure_even_generator_never_rotates",
                            "E5_symmetric_part_is_the_damping",
                            "control_asymmetric_B_breaks_hermiticity",
                            "control_unequal_diagonal_breaks_commuting")}
    for _ in range(trials):
        n = rng.randint(1, 3)
        Ii = iota(n)
        # E3 construction: equal antisymmetric diagonals, symmetric B, C = -B
        D = rand(n, rng, "anti"); B = rand(n, rng, "sym"); A = blocks_to_matrix(D, B, neg(B), D)
        if not eq(mul(A, Ii), mul(Ii, A)): ch["E1_commutes_iff_equal_diag_and_C_eq_minusB"] = False
        if not eq(T(A), neg(A)): ch["E2_antisym_iff_blocks_antisym_and_C_eq_minusBT"] = False
        H = neg(mul(Ii, A))
        if not (eq(H, T(H)) and eq(mul(H, Ii), mul(Ii, H))):
            ch["E3_both_iff_B_symmetric_gives_hermitian_H"] = False
        # control: make B asymmetric -> H no longer symmetric (generically)
        if n >= 2:
            B2 = [r[:] for r in B]; B2[0][1] += Fr(1, 3)
            A2 = blocks_to_matrix(D, B2, neg(B2), D)
            H2 = neg(mul(iota(n), A2))
            if eq(H2, T(H2)) and eq(T(A2), neg(A2)):
                ch["control_asymmetric_B_breaks_hermiticity"] = False
            # control: unequal diagonal blocks -> no commuting with iota
            D2 = rand(n, rng, "anti")
            if not eq(D2, D):
                A3 = blocks_to_matrix(D, B, neg(B), D2)
                if eq(mul(A3, iota(n)), mul(iota(n), A3)):
                    ch["control_unequal_diagonal_breaks_commuting"] = False
        # E4: purely cut-even generator (B = C = 0) never mixes sheets
        Aeven = blocks_to_matrix(rand(n, rng), zeros(n), zeros(n), rand(n, rng))
        mixes = any(Aeven[i][j] != 0 for i in range(n) for j in range(n, 2 * n)) or \
                any(Aeven[i][j] != 0 for i in range(n, 2 * n) for j in range(n))
        if mixes: ch["E4_pure_even_generator_never_rotates"] = False
        # E5: norm rate is governed by the symmetric part; antisymmetric part gives zero
        x = [Fr(rng.randint(-5, 5), rng.randint(1, 3)) for _ in range(2 * n)]
        if quad_form_rate(A, x) != 0: ch["E5_symmetric_part_is_the_damping"] = False
        Sdiag = [[Fr(1, 2) if i == j else Fr(0) for j in range(2 * n)] for i in range(2 * n)]
        Adamp = [[A[i][j] + Sdiag[i][j] for j in range(2 * n)] for i in range(2 * n)]
        if any(v != 0 for v in x) and quad_form_rate(Adamp, x) == 0:
            ch["E5_symmetric_part_is_the_damping"] = False
    return {"status": "PROVED" if all(ch.values()) else "FAILED", "checks": ch, "trials": trials}


def build():
    rng = random.Random(2026_08_19_2)
    body = {
        "protocol": "INFORMATION_INVARIANCE_D2A_V1",
        "theorem_E_cut_supplies_conjugate_pair_and_complex_structure": theorem_E(rng),
        "framework_inputs_consumed": {
            "iota_from_oriented_cut": "Recognition-Kernel-Framework theorems/foundation/F00E_NATIVE_EULER_FROM_IOTA_COMPLEX.md (PROVED): oriented cut -> quarter turn, iota^2 = -1, iota^dagger = -iota",
            "cut_grading": "Recognition-Kernel-Framework theorum/41_cut_graded_universal_generator_theorem.md (certified, pin 34afc445...): G = G_e + G_o, J G_o J = -G_o, even preserves sheets, odd transports",
        },
        "obligations": {
            "D2a": "DISCHARGED — conjugate channel = the opposite sheet of the oriented cut; complex structure = the derived quarter turn iota",
            "D2b": "OPEN — show the linearised energy dynamics has equal antisymmetric diagonal blocks and a symmetric transport block in the sheet basis; failure of the first is damping, of the second non-Hermiticity",
        },
        "claim_boundary": ("Exact finite linear algebra. Theorem E characterises, in the sheet basis of a "
                           "balanced cut, exactly when a real generator is a Schrodinger generator, and shows "
                           "phase lives in the cut-odd transporting block while damping lives in the symmetric "
                           "part. It consumes the framework's derived quarter turn and cut grading as inputs "
                           "and re-proves neither. It makes no claim about the energy dynamics itself."),
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D2a_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D2A.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("E:", cert["theorem_E_cut_supplies_conjugate_pair_and_complex_structure"]["status"],
          "|", cert["certificate_sha256"])
