"""D3 — where does the regulator length come from?

Obligation D3 asked: the claimed correlation
    <E(x)E(x')> ~ 1/(r^2 + L_P^2),  L_P = sqrt(hbar G / c^3),
inserts a length that contains hbar, inside a dynamics declared deterministic
and non-quantised. Derive it, or admit it was inserted.

This is a dimensional question and it is decidable exactly. Work in the (M, L, T)
dimension group; a quantity is an integer/rational exponent vector, products add
vectors, so "can a length be built from a set of constants" is a linear algebra
question over the rationals.

Constants available to the dynamics as written:
    c      = (0, 1, -1)
    G      = (-1, 3, -2)
    E      = (1, -1, -2)      energy density (the field, incl. E_vac)
    Lambda = (0, -2, 0)       cosmological constant, an inverse area
and, only if quantum theory is presupposed,
    hbar   = (1, 2, -1).

Theorem F (no Planck length without hbar).
  (i)  No product of powers of {c, G} is a length: G^a c^b = M^{-a} L^{3a+b} T^{-2a-b}
       and L^1 M^0 T^0 forces a = 0, then b = 1 and -b = 0, a contradiction.
  (ii) hbar is dimensionally independent of {c, G}: it is not in their span.
       Hence L_P = sqrt(hbar G / c^3) cannot be produced from {c, G} alone;
       it requires hbar as an input. Inserting L_P therefore presupposes the
       quantum constant the framework claims to derive.
  (iii) The vectors {c, G, hbar} are linearly independent and span the whole
       (M, L, T) group, so *with* hbar every dimension is expressible: the
       obstruction in (i)-(ii) is exactly the absence of hbar, nothing else.

Theorem G (the theory's own regulator is cosmological, not Planckian).
  The dynamics does contain a length that needs no hbar:
      L_E = c^2 / sqrt(G E)     from {c, G, E},
  and equivalently L_Lambda = Lambda^{-1/2} from the cosmological term. Both are
  exact solutions of the exponent system and are produced below. Numerically, for
  a vacuum energy density of order the observed one, L_E is of order the Hubble
  scale, i.e. ~10^61 Planck lengths. Hence:
      * the finiteness the model can honestly claim is regulated at a
        cosmological scale, not at L_P;
      * substituting L_P into the correlation function is an insertion, not a
        consequence, and the ultraviolet-finiteness claim does not follow.

Consequence for the paper: D3 is DISCHARGED with a negative result. The
correlation-function and emergent-Einstein conjecture must either (a) declare
hbar as an input constant --- which forfeits "deterministic, non-quantised" as
the derivation of quantum structure --- or (b) use its own scale L_E and drop
the claim of Planck-scale regularisation. Neither is a small edit; both are
stated here so a reader is not left with an unsupported UV claim.

Exact rational linear algebra; deterministic; certificate pinned.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
from fractions import Fraction as Fr

HERE = pathlib.Path(__file__).resolve().parent

# dimension vectors in (M, L, T)
DIM = {
    "c":      (Fr(0), Fr(1), Fr(-1)),
    "G":      (Fr(-1), Fr(3), Fr(-2)),
    "E":      (Fr(1), Fr(-1), Fr(-2)),     # energy density
    "Lambda": (Fr(0), Fr(-2), Fr(0)),
    "hbar":   (Fr(1), Fr(2), Fr(-1)),
    "length": (Fr(0), Fr(1), Fr(0)),
    "mass":   (Fr(1), Fr(0), Fr(0)),
}


def cj(o): return json.dumps(o, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
def sha(o): return hashlib.sha256(cj(o).encode()).hexdigest()


def solve(cols, target):
    """Solve sum_i x_i * cols[i] = target over Q by Gaussian elimination.
    Returns (solvable, one_solution_or_None, rank)."""
    n = len(cols)
    A = [[cols[i][r] for i in range(n)] + [target[r]] for r in range(3)]
    piv, where = 0, [-1] * n
    for col in range(n):
        sel = next((r for r in range(piv, 3) if A[r][col] != 0), None)
        if sel is None:
            continue
        A[piv], A[sel] = A[sel], A[piv]
        f = A[piv][col]
        A[piv] = [v / f for v in A[piv]]
        for r in range(3):
            if r != piv and A[r][col] != 0:
                g = A[r][col]
                A[r] = [A[r][k] - g * A[piv][k] for k in range(n + 1)]
        where[col] = piv
        piv += 1
    # inconsistent row?
    for r in range(3):
        if all(A[r][k] == 0 for k in range(n)) and A[r][n] != 0:
            return False, None, piv
    x = [Fr(0)] * n
    for col in range(n):
        if where[col] != -1:
            x[col] = A[where[col]][n]
    return True, x, piv


def check_product(names, exps, target):
    v = [Fr(0)] * 3
    for nm, e in zip(names, exps):
        d = DIM[nm]
        v = [v[k] + e * d[k] for k in range(3)]
    return tuple(v) == tuple(target)


def theorem_F():
    ok_i, sol_i, _ = solve([DIM["c"], DIM["G"]], DIM["length"])
    ok_h, sol_h, _ = solve([DIM["c"], DIM["G"]], DIM["hbar"])
    ok_full, sol_full, rank = solve([DIM["c"], DIM["G"], DIM["hbar"]], DIM["length"])
    lp_ok = ok_full and check_product(["c", "G", "hbar"], sol_full, DIM["length"])
    # verify the textbook Planck length exponents are the solution: hbar^1/2 G^1/2 c^-3/2
    lp_textbook = check_product(["hbar", "G", "c"], [Fr(1, 2), Fr(1, 2), Fr(-3, 2)], DIM["length"])
    return {
        "status": "PROVED" if (not ok_i and not ok_h and lp_ok and lp_textbook and rank == 3) else "FAILED",
        "checks": {
            "no_length_from_c_and_G_alone": not ok_i,
            "hbar_not_in_span_of_c_G": not ok_h,
            "length_solvable_once_hbar_admitted": lp_ok,
            "planck_exponents_verified": lp_textbook,
            "c_G_hbar_span_all_dimensions": rank == 3,
        },
        "planck_solution_exponents": {"c": str(sol_full[0]), "G": str(sol_full[1]), "hbar": str(sol_full[2])} if sol_full else None,
        "reading": "L_P requires hbar as an input; it is not derivable from the deterministic constants.",
    }


def theorem_G():
    ok_E, sol_E, _ = solve([DIM["c"], DIM["G"], DIM["E"]], DIM["length"])
    LE_ok = ok_E and check_product(["c", "G", "E"], sol_E, DIM["length"])
    # canonical form c^2 (G E)^(-1/2)
    canonical = check_product(["c", "G", "E"], [Fr(2), Fr(-1, 2), Fr(-1, 2)], DIM["length"])
    ok_L, sol_L, _ = solve([DIM["Lambda"]], DIM["length"])
    lam_ok = ok_L and check_product(["Lambda"], sol_L, DIM["length"])
    # numeric magnitude: rho_vac ~ 6e-10 J/m^3, G=6.674e-11, c=3e8  -> L_E in metres, and in Planck lengths
    import math
    G, c, rho = 6.674e-11, 2.998e8, 6.0e-10
    L_E = c ** 2 / math.sqrt(G * rho)
    L_P = math.sqrt(1.0546e-34 * G / c ** 3)
    ratio = L_E / L_P
    return {
        "status": "PROVED" if (LE_ok and canonical and lam_ok) else "FAILED",
        "checks": {"length_from_c_G_E_exists": LE_ok,
                   "canonical_form_c2_over_sqrt_GE": canonical,
                   "lambda_gives_a_length": lam_ok},
        "own_regulator": "L_E = c^2 / sqrt(G * E)  (equivalently Lambda^{-1/2})",
        "numeric_illustration": {"L_E_metres": f"{L_E:.3e}", "L_P_metres": f"{L_P:.3e}",
                                 "L_E_over_L_P": f"{ratio:.3e}",
                                 "note": "order Hubble scale; ~10^61 Planck lengths"},
        "reading": ("The only length the deterministic dynamics generates by itself is cosmological. "
                    "Regularisation at L_P is an insertion, and the ultraviolet-finiteness claim does not follow."),
    }


def build():
    F, G_ = theorem_F(), theorem_G()
    body = {
        "protocol": "INFORMATION_INVARIANCE_D3_V1",
        "theorem_F_no_planck_length_without_hbar": F,
        "theorem_G_own_regulator_is_cosmological": G_,
        "obligation": {"D3": "DISCHARGED WITH NEGATIVE RESULT"},
        "verdict": ("The correlation function 1/(r^2 + L_P^2) cannot be obtained from the stated "
                    "deterministic constants: hbar is dimensionally independent of c and G, so L_P is an "
                    "input, not a consequence. The scale the dynamics does generate, c^2/sqrt(G E), is "
                    "cosmological. The paper must therefore either admit hbar as a primitive — forfeiting "
                    "the claim to derive quantum structure — or drop Planck-scale regularisation."),
        "scope": "Bears on the emergent-gravity conjecture only; Theorems A and B are unaffected.",
        "claim_boundary": "Exact rational linear algebra in the (M, L, T) dimension group; no physics beyond dimensional analysis is used or claimed.",
    }
    body["certificate_sha256"] = sha(body)
    return body


if __name__ == "__main__":
    cert = build()
    (HERE / "D3_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    (HERE / "EXPECTED_D3.sha256").write_text(cert["certificate_sha256"] + "\n")
    print("F:", cert["theorem_F_no_planck_length_without_hbar"]["status"],
          "| G:", cert["theorem_G_own_regulator_is_cosmological"]["status"],
          "|", cert["theorem_G_own_regulator_is_cosmological"]["numeric_illustration"]["L_E_over_L_P"],
          "|", cert["certificate_sha256"])
