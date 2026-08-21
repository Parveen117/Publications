"""CFE-U: uniqueness of the obstruction — the flagship's last half.

The Cut-First Equivalence theorem had one remaining open obligation:

    (U) UNIQUENESS. iint Omega is the UNIQUE obstruction to the
        chi -> 1 limit: no second, independent obstruction exists.

This capsule discharges (U) on the DECLARED RESPONSE ALGEBRA — the same
algebra whose classical face CFE-2 certified at rank 5 — by the same
exact-rank machinery that closed (S). Everything is exact rational
linear algebra; no floats, no limits, no transcendental evaluation.

THE RESPONSE COMPLEX

    Lambda^0  --d0-->  Lambda^1  --d1-->  Lambda^2

    Lambda^0 = quadratic potentials      span{1, p, v, p^2, pv, v^2}   (6)
    Lambda^1 = affine-coefficient forms  (a+bp+cv)dp + (e+fp+gv)dv     (6)
    Lambda^2 = reachable curvatures      Q * (dp ^ dv)                 (1)

The CFE witness family omega(chi) = (3p+5v)dp + (2v+5*chi*p)dv lives in
Lambda^1 for EVERY chi, and its classical (chi=1) face generates exactly
the rank-5 classical response space of CFE-2 (= rank d0).

BLOCKS

  T1  THE COMPLEX AND ITS EXACT RANKS. rank d0 = 5 (kernel = the
      constants, dim 1), rank d1 = 1, and d1 o d0 = 0. All as exact
      integer-matrix ranks over Q.

  T2  CLOSED = EXACT (H^1 = 0). ker d1 = im d0, certified two ways:
      the rank equality dim ker d1 = 5 = rank d0, AND constructively —
      the explicit potential Phi = a p + b p^2/2 + c p v + e v +
      g v^2/2 integrates EVERY closed form, verified at the coefficient
      level. Consequence: a response with zero curvature leaves NO
      residual obstruction — it IS memoryless (a gradient). There is no
      hidden obstruction at level one.

  T3  THE UNIQUENESS THEOREM (H^2(response) rank 1). The space of
      obstruction functionals — linear functionals on Lambda^1 that
      vanish on every memoryless response — is the annihilator of
      ker d1, and its dimension is EXACTLY 1 (exact nullspace
      computation). The loop-residue functional R(omega) =
      oint_{dD} omega on the CFE-1 rectangle lies in it (Stokes,
      verified on every kernel basis vector) and is nonzero, hence
      SPANS: every obstruction to the memoryless limit is c * R for a
      scalar c. Moreover R = Area * d1 as linear maps on the whole
      basis: the residue and the curvature flux are the same rank-1
      map. iint Omega is THE obstruction — there is no second one.

  T4  THE DIAL EXCITES THE GENERATOR, AND NOTHING ELSE EXISTS TO
      EXCITE. The chi-derivative of the witness family is
      d omega/d chi = 5 p dv, whose curvature is 5 != 0: the memory
      dial drives the one obstruction dimension. The residue of
      omega(chi) is EXACTLY 80*(chi - 1), reproducing CFE-1's
      certified ladder -32, -16, 0, +16, +32 at chi = 3/5 ... 7/5 —
      consumed directly from the pinned CFE-1 module, not recomputed
      from a copy. UNIQUENESS IN ACTION: two deformations with equal
      residue differ by a closed form, and T2 integrates that
      difference to an explicit potential — equal residue IS equal
      obstruction class.

  T5  THE DISCRETE SIDE, RANK 1 CELL BY CELL. On the exact 4x4 cell
      grid over the CFE-1 rectangle, the (cells x basis) residue
      matrix — every entry an exact midpoint-rule circulation, which
      is EXACT for affine integrands — has rank EXACTLY 1, every cell
      row equal to (cell area) * d1, and kernel exactly ker d1. The
      per-cell obstruction is one number repeated, not a family of
      independent obstructions.

  T6  THE BOUNDARY IS REAL (CONTROL). Enlarging the algebra to
      quadratic-coefficient one-forms (dim 12) raises rank d1 to 3:
      in the larger algebra there are three independent curvature
      directions. Uniqueness is therefore a THEOREM ABOUT THE DECLARED
      RESPONSE ALGEBRA — the one whose classical face CFE-1/CFE-2
      certified — not a universal law, and the certificate proves
      where the boundary sits instead of hiding it.

WHAT THIS CLOSES. With CFE-1 (memoryless limit + residue), CFE-2
(surjectivity, rank 5 = 5, cokernel 0) and CFE-U (uniqueness,
obstruction space rank 1), all declared obligations of the Cut-First
Equivalence theorem are certified on the witness EOS and its declared
response algebra. Continuum, general-manifold and general-EOS
uniqueness are NOT claimed; the EOS is the witness, not the theorem.
Nothing here touches RH / K0 / L0 / YM continuum gates.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfe1 = _load("cfe1_for_u", "cfe1_cut_first_equivalence.py")

# ----------------------------------------------------------------------
# exact rank machinery over Q
# ----------------------------------------------------------------------


def rank(M):
    """Exact rank of a rational matrix by Gaussian elimination."""
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0]) if A else 0
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        A[r] = [x / A[r][c] for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                A[i] = [x - A[i][c] * y for x, y in zip(A[i], A[r])]
        r += 1
    return r


def nullspace(M):
    """Exact basis of the right nullspace of a rational matrix."""
    A = [row[:] for row in M]
    rows, cols = len(A), len(A[0]) if A else 0
    pivots = []
    r = 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        A[r] = [x / A[r][c] for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                A[i] = [x - A[i][c] * y for x, y in zip(A[i], A[r])]
        pivots.append(c)
        r += 1
    free = [c for c in range(cols) if c not in pivots]
    basis = []
    for fc in free:
        vec = [Fr(0)] * cols
        vec[fc] = Fr(1)
        for i, pc in enumerate(pivots):
            vec[pc] = -A[i][fc]
        basis.append(vec)
    return basis


def mv(M, x):
    return [sum(M[i][j] * x[j] for j in range(len(x)))
            for i in range(len(M))]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


# ----------------------------------------------------------------------
# the response complex on the declared algebra
# ----------------------------------------------------------------------
#
# Lambda^0 basis: [1, p, v, p^2, pv, v^2]
# Lambda^1 coordinates: (a, b, c, e, f, g) for (a+bp+cv)dp + (e+fp+gv)dv
# Lambda^2: multiples of dp ^ dv

D0 = [
    # columns: 1, p, v, p^2, pv, v^2 ; rows: a, b, c, e, f, g
    [Fr(0), Fr(1), Fr(0), Fr(0), Fr(0), Fr(0)],   # a  <- d/dp
    [Fr(0), Fr(0), Fr(0), Fr(2), Fr(0), Fr(0)],   # b
    [Fr(0), Fr(0), Fr(0), Fr(0), Fr(1), Fr(0)],   # c
    [Fr(0), Fr(0), Fr(1), Fr(0), Fr(0), Fr(0)],   # e  <- d/dv
    [Fr(0), Fr(0), Fr(0), Fr(0), Fr(1), Fr(0)],   # f
    [Fr(0), Fr(0), Fr(0), Fr(0), Fr(0), Fr(2)],   # g
]

D1 = [[Fr(0), Fr(0), Fr(-1), Fr(0), Fr(1), Fr(0)]]   # f - c

BASIS1 = ("dp", "p dp", "v dp", "dv", "p dv", "v dv")

P_LO, P_HI, V_LO, V_HI = cfe1.LOOP_RECT
AREA = (P_HI - P_LO) * (V_HI - V_LO)


def form_coeffs(vec):
    a, b, c, e, f, g = vec

    def f_dp(p, v):
        return a + b * p + c * v

    def f_dv(p, v):
        return e + f * p + g * v

    return f_dp, f_dv


def circulation_of(vec, loop):
    """Exact midpoint circulation of an affine one-form (midpoint rule is
    EXACT for affine integrands)."""
    f_dp, f_dv = form_coeffs(vec)
    total = Fr(0)
    n = len(loop)
    for i in range(n):
        p0, v0 = loop[i]
        p1, v1 = loop[(i + 1) % n]
        pm, vm = Fr(p0 + p1, 2), Fr(v0 + v1, 2)
        total += f_dp(pm, vm) * (p1 - p0) + f_dv(pm, vm) * (v1 - v0)
    return total


def rect_loop(p_lo, p_hi, v_lo, v_hi):
    return [(p_lo, v_lo), (p_hi, v_lo), (p_hi, v_hi), (p_lo, v_hi)]


def witness_vec(chi):
    """The CFE-1 witness omega(chi) as a Lambda^1 coordinate vector,
    read off the PINNED cfe1 coefficient functions (not re-declared)."""
    a = cfe1.lambda_p(Fr(0), Fr(0), chi)
    b = cfe1.lambda_p(Fr(1), Fr(0), chi) - a
    c = cfe1.lambda_p(Fr(0), Fr(1), chi) - a
    e = cfe1.lambda_v(Fr(0), Fr(0), chi)
    f = cfe1.lambda_v(Fr(1), Fr(0), chi) - e
    g = cfe1.lambda_v(Fr(0), Fr(1), chi) - e
    # affine check against a mixed point
    assert cfe1.lambda_p(Fr(2), Fr(3), chi) == a + 2 * b + 3 * c
    assert cfe1.lambda_v(Fr(2), Fr(3), chi) == e + 2 * f + 3 * g
    return [a, b, c, e, f, g]


# ----------------------------------------------------------------------
# T1  the complex and its exact ranks
# ----------------------------------------------------------------------


def certify_T1():
    r0, r1 = rank(D0), rank(D1)
    assert r0 == 5 and r1 == 1
    ker0 = nullspace(D0)
    assert len(ker0) == 1 and ker0[0][0] != 0     # constants only
    comp = [mv(D1, [D0[i][j] for i in range(6)]) for j in range(6)]
    assert all(x == [Fr(0)] for x in comp)        # d1 o d0 = 0
    return {
        "statement": (
            "The response complex Lambda^0 -> Lambda^1 -> Lambda^2 on "
            "the declared algebra has exact ranks: rank d0 = 5 (kernel "
            "= the constants), rank d1 = 1, d1 o d0 = 0. The rank-5 "
            "image of d0 IS the classical response space certified by "
            "CFE-2"),
        "rank_d0": r0,
        "rank_d1": r1,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  closed = exact: H^1 = 0, constructively
# ----------------------------------------------------------------------


def potential_of_closed(vec):
    a, b, c, e, f, g = vec
    assert f == c                                 # closed
    # Phi = a p + b p^2/2 + c p v + e v + g v^2/2 in Lambda^0 coords
    return [Fr(0), a, e, b / 2, c, g / 2]


def certify_T2():
    ker1 = nullspace(D1)
    assert len(ker1) == 5                        # dim ker d1 = rank d0
    for k in ker1:
        phi = potential_of_closed(k)
        assert mv_transpose_d0(phi) == k         # d0(Phi) = the form
    # a mixed rational closed form, integrated explicitly
    w = [Fr(3, 2), Fr(-1), Fr(7, 3), Fr(0), Fr(7, 3), Fr(5)]
    assert mv(D1, w) == [Fr(0)]
    phi = potential_of_closed(w)
    assert mv_transpose_d0(phi) == w
    return {
        "statement": (
            "H^1(response) = 0: ker d1 = im d0, by the rank equality "
            "5 = 5 AND constructively — the explicit potential Phi = "
            "a p + b p^2/2 + c p v + e v + g v^2/2 integrates every "
            "closed form at the coefficient level. Zero curvature "
            "leaves NO residual obstruction: the response is a "
            "gradient, i.e. memoryless"),
        "dim_ker_d1": len(ker1),
        "verdict": "PASS",
    }


def mv_transpose_d0(phi):
    """d0 applied to a Lambda^0 vector: D0 is stored rows=Lambda^1
    coords, columns=Lambda^0 basis, so this is D0 * phi."""
    return mv(D0, phi)


# ----------------------------------------------------------------------
# T3  the uniqueness theorem: obstruction space has rank exactly 1
# ----------------------------------------------------------------------


def certify_T3():
    ker1 = nullspace(D1)
    # annihilator of ker d1 inside (Lambda^1)*: functionals L with
    # L(k) = 0 for every kernel vector — nullspace of the kernel matrix
    ann = nullspace([k for k in ker1])
    assert len(ann) == 1                          # THE uniqueness rank

    # the residue functional on the CFE-1 rectangle
    loop = rect_loop(P_LO, P_HI, V_LO, V_HI)
    R = []
    for j in range(6):
        basis = [Fr(0)] * 6
        basis[j] = Fr(1)
        R.append(circulation_of(basis, loop))
    # R annihilates every memoryless response (Stokes on the kernel)
    for k in ker1:
        assert dot(R, k) == 0
    assert any(x != 0 for x in R)                 # nonzero
    # hence R spans the 1-dimensional annihilator: check proportionality
    gen = ann[0]
    ratio = None
    for rj, gj in zip(R, gen):
        if gj != 0:
            ratio = rj / gj
            break
    assert ratio is not None and ratio != 0
    assert all(rj == ratio * gj for rj, gj in zip(R, gen))

    # and R = Area * d1 as linear maps: residue IS curvature flux
    assert R == [AREA * D1[0][j] for j in range(6)]
    assert R == [Fr(0), Fr(0), -AREA, Fr(0), AREA, Fr(0)]

    return {
        "statement": (
            "THE UNIQUENESS THEOREM: the space of obstruction "
            "functionals — linear functionals vanishing on every "
            "memoryless response — is the annihilator of ker d1 and "
            "has dimension EXACTLY 1. The loop-residue functional "
            "R(omega) = oint_(dD) omega lies in it (Stokes verified on "
            "the kernel basis), is nonzero, and therefore SPANS: every "
            "obstruction to the chi -> 1 limit is a scalar multiple of "
            "the residue. Moreover R = Area * d1 exactly on the whole "
            "basis: the residue and the curvature flux iint Omega are "
            "the same rank-1 map. There is no second, independent "
            "obstruction"),
        "annihilator_dimension": len(ann),
        "residue_functional": [str(x) for x in R],
        "area": str(AREA),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  the dial excites the generator; equal residue = equal class
# ----------------------------------------------------------------------

CHIS = (Fr(3, 5), Fr(4, 5), Fr(1), Fr(6, 5), Fr(7, 5))


def certify_T4():
    # the chi-derivative of the witness family is 5 p dv
    w1 = witness_vec(Fr(2))
    w0 = witness_vec(Fr(1))
    deriv = [x - y for x, y in zip(w1, w0)]       # linear in chi: exact
    assert deriv == [Fr(0)] * 4 + [cfe1.BETA, Fr(0)]
    assert mv(D1, deriv) == [cfe1.BETA]           # excites the generator

    # the residue ladder, consumed from the PINNED cfe1 circulation
    loop = rect_loop(P_LO, P_HI, V_LO, V_HI)
    ladder = []
    for chi in CHIS:
        mine = circulation_of(witness_vec(chi), loop)
        pinned = cfe1.circulation(loop, chi)
        assert mine == pinned                     # two routes agree
        assert mine == cfe1.BETA * (chi - 1) * AREA
        ladder.append((chi, mine))
    assert [x for _, x in ladder] == [Fr(-32), Fr(-16), Fr(0),
                                      Fr(16), Fr(32)]

    # uniqueness in action: equal residue => difference is closed =>
    # integrate it to an explicit potential
    u = witness_vec(Fr(6, 5))
    w = [x + y for x, y in zip(
        witness_vec(Fr(6, 5)),
        # add a closed perturbation: gradient of p*v + 3*p
        mv(D0, [Fr(0), Fr(3), Fr(0), Fr(0), Fr(1), Fr(0)]))]
    assert circulation_of(u, loop) == circulation_of(w, loop)
    diff = [x - y for x, y in zip(w, u)]
    assert mv(D1, diff) == [Fr(0)]
    phi = potential_of_closed(diff)
    assert mv(D0, phi) == diff                    # explicit witness

    return {
        "statement": (
            "The memory dial's deformation direction is d omega/d chi "
            "= 5 p dv with curvature 5 != 0: the dial drives the one "
            "obstruction dimension. The residue of omega(chi) equals "
            "BETA*(chi-1)*Area = 80*(chi-1) exactly, reproducing "
            "CFE-1's certified ladder -32, -16, 0, +16, +32 — computed "
            "here AND consumed from the pinned CFE-1 circulation, two "
            "routes agreeing as rationals. Uniqueness in action: two "
            "responses with equal residue differ by a closed form, "
            "integrated to an explicit potential — equal residue IS "
            "equal obstruction class"),
        "ladder": [[str(c), str(x)] for c, x in ladder],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  discrete side: the residue matrix has rank 1, cell by cell
# ----------------------------------------------------------------------


def certify_T5():
    cells = []
    for i in range(4):
        for j in range(4):
            p_lo = P_LO + i * (P_HI - P_LO) / 4
            p_hi = P_LO + (i + 1) * (P_HI - P_LO) / 4
            v_lo = V_LO + j * (V_HI - V_LO) / 4
            v_hi = V_LO + (j + 1) * (V_HI - V_LO) / 4
            cells.append((p_lo, p_hi, v_lo, v_hi))
    cell_area = (P_HI - P_LO) / 4 * (V_HI - V_LO) / 4

    M = []
    for cell in cells:
        loop = rect_loop(*cell)
        row = []
        for j in range(6):
            basis = [Fr(0)] * 6
            basis[j] = Fr(1)
            row.append(circulation_of(basis, loop))
        M.append(row)

    assert rank(M) == 1
    expected_row = [cell_area * D1[0][j] for j in range(6)]
    assert all(row == expected_row for row in M)  # one obstruction,
    ker_M = nullspace(M)                          # repeated per cell
    ker1 = nullspace(D1)
    assert len(ker_M) == len(ker1) == 5
    assert rank(ker_M + ker1) == 5                # same 5-dim kernel

    # discrete Stokes for the witness at several chis, every cell
    for chi in (Fr(3, 5), Fr(1), Fr(7, 5)):
        wv = witness_vec(chi)
        for cell, row in zip(cells, M):
            assert dot(row, wv) == cfe1.BETA * (chi - 1) * cell_area

    return {
        "statement": (
            "On the exact 4x4 cell grid the (cells x basis) residue "
            "matrix — every entry an exact midpoint circulation — has "
            "rank EXACTLY 1, every cell row equal to (cell area) * d1, "
            "and kernel exactly ker d1 (the same 5-dimensional "
            "memoryless space). The per-cell obstruction is ONE number "
            "repeated, not a family of independent obstructions; "
            "discrete Stokes holds for the witness at every cell and "
            "every sampled chi"),
        "cells": len(cells),
        "residue_matrix_rank": 1,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  the boundary control: a larger algebra has rank 3
# ----------------------------------------------------------------------


def certify_T6():
    # quadratic-coefficient one-forms: coefficients in span{1,p,v,p^2,
    # pv,v^2} for each of dp, dv -> dim 12. Curvature density lives in
    # span{1, p, v}. d1_big columns = basis forms, rows = (1, p, v).
    def d_of(dp_coeff, dv_coeff):
        # returns density coords (const, p, v) of d(dp_c dp + dv_c dv)
        # density = d(dv_c)/dp - d(dp_c)/dv, computed on monomials
        # derivative tables on basis [1,p,v,p2,pv,v2] -> (const,p,v)
        ddp = {0: (0, 0, 0), 1: (1, 0, 0), 2: (0, 0, 0),
               3: (0, 2, 0), 4: (0, 0, 1), 5: (0, 0, 0)}
        ddv = {0: (0, 0, 0), 1: (0, 0, 0), 2: (1, 0, 0),
               3: (0, 0, 0), 4: (0, 1, 0), 5: (0, 0, 2)}
        out = [Fr(0), Fr(0), Fr(0)]
        for k in range(6):
            for t in range(3):
                out[t] += Fr(ddp[k][t]) * dv_coeff[k]
                out[t] -= Fr(ddv[k][t]) * dp_coeff[k]
        return out

    cols = []
    for part in range(2):
        for k in range(6):
            dp_c = [Fr(0)] * 6
            dv_c = [Fr(0)] * 6
            (dp_c if part == 0 else dv_c)[k] = Fr(1)
            cols.append(d_of(dp_c, dv_c))
    D1_big = [[cols[j][t] for j in range(12)] for t in range(3)]
    r_big = rank(D1_big)
    assert r_big == 3                             # three directions now

    # the declared algebra sits inside and still has rank 1
    sub_cols = [0, 1, 2, 6, 7, 8]                 # affine coefficients
    D1_sub = [[D1_big[t][j] for j in sub_cols] for t in range(3)]
    assert rank(D1_sub) == 1

    return {
        "statement": (
            "CONTROL: in the enlarged algebra of quadratic-coefficient "
            "one-forms the curvature map has rank 3 — three independent "
            "obstruction directions exist there — while the declared "
            "affine response algebra embedded inside it still has rank "
            "1. Uniqueness is therefore a THEOREM ABOUT THE DECLARED "
            "RESPONSE ALGEBRA (the one whose classical face CFE-1 and "
            "CFE-2 certified), not a universal law, and this "
            "certificate proves where the boundary sits"),
        "rank_enlarged": r_big,
        "rank_declared_inside": 1,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "CFE-U: uniqueness of the obstruction — the "
                   "flagship's last half",
        "consumes": {
            "CFE-1": ("witness family lambda_p, lambda_v, BETA, "
                      "LOOP_RECT and circulation consumed from the "
                      "pinned module, not re-declared"),
            "CFE-2": ("the rank-5 classical response space = im d0 of "
                      "this capsule's complex"),
        },
        "T1_response_complex_exact_ranks": certify_T1(),
        "T2_closed_equals_exact": certify_T2(),
        "T3_uniqueness_obstruction_rank_1": certify_T3(),
        "T4_dial_excites_generator_ladder": certify_T4(),
        "T5_discrete_residue_matrix_rank_1": certify_T5(),
        "T6_boundary_control_rank_3": certify_T6(),
        "flagship_status": {
            "part_1_memoryless_limit": "CERTIFIED (CFE-1 T1, T4)",
            "part_2_residue": "CERTIFIED (CFE-1 T2, T3)",
            "part_S_surjectivity": "CERTIFIED (CFE-2)",
            "part_U_uniqueness": "CERTIFIED (this capsule)",
            "scope": (
                "all declared obligations of the Cut-First Equivalence "
                "theorem are now certified ON THE WITNESS EOS AND ITS "
                "DECLARED RESPONSE ALGEBRA"),
        },
        "claim_boundary": {
            "generality": (
                "NOT CLAIMED: continuum, general-manifold and "
                "general-EOS uniqueness. The declared algebra is the "
                "affine-coefficient response space; T6 certifies that "
                "larger algebras carry more obstruction directions"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": (
                "NONE — (U) was the flagship's single remaining open "
                "obligation since CFE-2"),
            "machinery": (
                "the same exact-rank machinery that closed (S): "
                "rational Gaussian elimination, exact nullspaces, "
                "midpoint circulations exact for affine integrands"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "CFEU_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    print("CFE-U certificate written:", out)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
