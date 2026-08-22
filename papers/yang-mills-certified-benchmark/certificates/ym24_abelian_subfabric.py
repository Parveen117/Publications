"""YM-24: THE ABELIAN SUBFABRIC AND THE NON-ABELIAN RESIDUE — the exact
structure of the commutator sector that YM-23 named as the only possible
carrier of a weak-coupling gap. Everything here is exact algebra in the
EMK block over C_Sigma (YM-F1); nothing is a gap claim.

 (T1) THE ABELIAN SUBFABRIC IS F00-E'S CIRCULAR EULER ORBIT. Fix a unit
      odd direction u (u^2 = -1 in the quaternion algebra). Faces of the
      form H = c + s u with c^2 + s^2 = 1 commute, and compose by the
      NATIVE ADDITION LAW of Exp_Sigma(iota theta) (F00-E Thm 3.1 /
      Thm 5.2): (c1 + s1 u)(c2 + s2 u) = (c1 c2 - s1 s2) + (s1 c2 + c1 s2) u.
      On such a fabric Recognition-Stokes is ADDITIVE in the Euler
      coordinate — the boundary holonomy's (c, s) is the angle-sum of
      the faces — certified exactly on rational circle points
      (Pythagorean, no trig) for m = 2..8. The abelian subfabric of
      the chain is exactly the compact-rotor chain, and every such
      subfabric is a copy of the same native circle.

 (T2) NON-ABELIAN RESIDUE = GRAM DEFECT (native Lagrange identity).
      For odd p, q:  |p x q|^2 = |p|^2 |q|^2 - <p, q>^2   exactly, so the
      commutator residue of two faces vanishes IFF the faces lie in one
      abelian subfabric (p parallel q), and is otherwise strictly
      positive. The "non-healable" content of a pair of faces is the
      Gram determinant of their odd coordinates — the same object as
      EMK-1's rotational determinant channel, here between faces.

 (T3) EXACT SECOND-ORDER STOKES. For faces H_i = 1 + p_i (unnormalised
      chart), the product's odd part to second order is
          odd_2 = sum_i p_i + sum_{i<j} p_i x p_j,
      and its even part is 1 - sum_{i<j} <p_i, p_j>; the exact remainder
      is of degree >= 3 in the face coordinates (certified by exact
      polynomial extraction in a scaling parameter, m = 2..6, random
      rational directions). The second-order non-abelian term is the
      SUM OF PAIRWISE CROSS PRODUCTS OF ALL FACE PAIRS — not only
      neighbours: once faces are multiplied along the chain, every pair
      (i<j) contributes, ordered. Stokes additivity fails at second
      order by exactly this bivector.

 (T4) THE SOFT MODE IS ABELIAN; THE NON-ABELIAN ENERGY IS A POSITIVE
      FORM THAT VANISHES ON IT. YM-23's softest spatial mode v (linear
      profile) has all face coordinates PARALLEL (p_i = v_i u), hence
      lies in one abelian subfabric (T1), hence has non-abelian residue
      exactly zero (T2) — the gapless direction of the abelianized
      fabric is not an artefact of linearisation: it is an exact
      abelian subfabric of the full chain. The non-abelian energy
      E_na(p) = |sum_{i<j} p_i x p_j|^2 is a positive semidefinite
      quartic; certified: E_na = 0 on every parallel configuration,
      E_na > 0 on generic ones, and E_na is invariant under a common
      rotation of all faces (exact, rational rotations).

WHAT THIS MEANS, NATIVELY (no claim beyond the algebra): any
weak-coupling gap of the chain must come from the dynamics that
connects DIFFERENT abelian subfabrics — transitions that change the
direction u — because inside any single subfabric the chain is the
exactly additive rotor chain with a gapless soft mode. The named next
object (YM-25): the transfer restricted to direction changes — the
"u-rotation" sector — and whether the time kernel's action on it has a
floor uniform in m. That is the non-perturbative content, now isolated
to one sector.

Controls:
  C1  circle composition = native addition law on rational circle
      points; Stokes additive in the abelian subfabric, m = 2..8.
  C2  Lagrange identity exact on a rational grid; zero iff parallel.
  C3  degree-2 truncation exact (coefficients of eps^1, eps^2 match;
      eps^0 matches), m = 2..6.
  C4  parallel control: all cross products vanish, odd_2 = sum, even_2
      = 1 - sum of dot products.
  C5  E_na = 0 on parallel, > 0 on generic, rotation-invariant.
  C6  tamper: swapping the order of two non-commuting faces changes the
      second-order bivector by exactly -2 p_i x p_j (order is content,
      EMK-T2).
"""

from fractions import Fraction as F
import json
import os
import random
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import _dec, canonical_sha  # noqa: E402
from ymf1_chain_fabric import Quat, rational_unit, QONE  # noqa: E402


def cross(p, q):
    return (p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2],
            p[0] * q[1] - p[1] * q[0])


def dot(p, q):
    return sum(p[i] * q[i] for i in range(3))


def vadd(p, q):
    return tuple(p[i] + q[i] for i in range(3))


def vscale(p, c):
    return tuple(c * x for x in p)


def circle_point(t: F):
    """rational point on the unit circle: ((1-t^2)/(1+t^2), 2t/(1+t^2))."""
    return ((1 - t * t) / (1 + t * t), 2 * t / (1 + t * t))


def abelian_face(c: F, s: F, u):
    return Quat(c, s * u[0], s * u[1], s * u[2])


def rnd_vec():
    return tuple(F(random.randint(-6, 6), random.randint(1, 5)) for _ in range(3))


def product_faces(ps, eps: F):
    out = QONE
    for p in ps:
        out = out * Quat(1, eps * p[0], eps * p[1], eps * p[2])
    return out


def poly_coeffs(values, xs):
    """exact interpolation: coefficients of polynomial through (xs, values)."""
    n = len(xs)
    # solve Vandermonde exactly (Gaussian elimination in Q)
    A = [[xs[i] ** j for j in range(n)] + [values[i]] for i in range(n)]
    for col in range(n):
        piv = next(r for r in range(col, n) if A[r][col] != 0)
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [A[r][k] - f * A[col][k] for k in range(n + 1)]
    return [A[i][n] for i in range(n)]


def second_order(ps):
    odd = (F(0), F(0), F(0))
    even = F(1)
    for i, p in enumerate(ps):
        odd = vadd(odd, p)
        for j in range(i + 1, len(ps)):
            odd = vadd(odd, cross(p, ps[j]))
            even -= dot(p, ps[j])
    return odd, even


def e_na(ps):
    b = (F(0), F(0), F(0))
    for i in range(len(ps)):
        for j in range(i + 1, len(ps)):
            b = vadd(b, cross(ps[i], ps[j]))
    return dot(b, b)


def rotate(p, g: Quat):
    """conjugation by a unit quaternion: exact rational rotation."""
    q = g * Quat(0, p[0], p[1], p[2]) * g.inv()
    return (q.b, q.c, q.d)


def run():
    random.seed(20260822)
    # ---- T1 / C1 abelian subfabric = circular Euler orbit
    c1 = True
    u = (F(3, 7), F(2, 7), F(6, 7))             # 9+4+36 = 49 -> unit
    assert dot(u, u) == 1
    for m in range(2, 9):
        ts = [F(random.randint(-5, 5), random.randint(1, 4)) for _ in range(m)]
        faces = [abelian_face(*circle_point(t), u) for t in ts]
        # native addition law pairwise
        for i in range(m - 1):
            c1_, s1_ = circle_point(ts[i])
            c2_, s2_ = circle_point(ts[i + 1])
            prod = faces[i] * faces[i + 1]
            want = abelian_face(c1_ * c2_ - s1_ * s2_, s1_ * c2_ + c1_ * s2_, u)
            if prod != want:
                c1 = False
            if faces[i] * faces[i + 1] != faces[i + 1] * faces[i]:
                c1 = False
        # Stokes additive: total = angle sum, computed by iterated law
        c_tot, s_tot = F(1), F(0)
        for t in ts:
            c_, s_ = circle_point(t)
            c_tot, s_tot = c_tot * c_ - s_tot * s_, s_tot * c_ + c_tot * s_
        prod = QONE
        for f in faces:
            prod = prod * f
        if prod != abelian_face(c_tot, s_tot, u):
            c1 = False

    # ---- T2 / C2 Lagrange identity
    c2 = True
    for _ in range(40):
        p, q = rnd_vec(), rnd_vec()
        x = cross(p, q)
        if dot(x, x) != dot(p, p) * dot(q, q) - dot(p, q) ** 2:
            c2 = False
    p = rnd_vec()
    if dot(cross(p, vscale(p, F(-7, 3))), cross(p, vscale(p, F(-7, 3)))) != 0:
        c2 = False

    # ---- T3 / C3 exact degree-2 truncation via interpolation
    c3 = True
    for m in range(2, 7):
        ps = [rnd_vec() for _ in range(m)]
        xs = [F(k, 7) for k in range(1, m + 2)]          # degree <= m
        vals = [product_faces(ps, e) for e in xs]
        odd2, even2 = second_order(ps)
        for comp, want1, want2 in [("a", F(0), even2 - 1),
                                   ("b", odd2[0], None), ("c", odd2[1], None),
                                   ("d", odd2[2], None)]:
            coeffs = poly_coeffs([getattr(v, comp) for v in vals], xs)
            if comp == "a":
                if not (coeffs[0] == 1 and coeffs[1] == 0 and coeffs[2] == want2):
                    c3 = False
            else:
                # odd: eps^1 coefficient = sum p ; eps^2 = sum cross
                lin = sum((p[{"b": 0, "c": 1, "d": 2}[comp]] for p in ps), F(0))
                quad = want1 - lin
                if not (coeffs[0] == 0 and coeffs[1] == lin and coeffs[2] == quad):
                    c3 = False
    # ---- C4 parallel control
    base = rnd_vec()
    ps = [vscale(base, F(k, 3)) for k in range(1, 6)]
    odd2, even2 = second_order(ps)
    c4 = odd2 == tuple(sum((p[i] for p in ps), F(0)) for i in range(3)) and \
        even2 == 1 - sum(dot(ps[i], ps[j]) for i in range(5) for j in range(i + 1, 5))
    # ---- T4 / C5 non-abelian energy
    c5 = e_na(ps) == 0
    gen = [rnd_vec() for _ in range(5)]
    c5 = c5 and e_na(gen) > 0
    g = rational_unit(F(1, 2), F(-1, 3), F(2, 5))
    c5 = c5 and e_na([rotate(p, g) for p in gen]) == e_na(gen)
    # ---- C6 order tamper
    q1, q2 = gen[0], gen[1]
    o12, _ = second_order([q1, q2])
    o21, _ = second_order([q2, q1])
    diff = tuple(o12[i] - o21[i] for i in range(3))
    c6 = diff == vscale(cross(q1, q2), F(2)) and diff != (0, 0, 0)

    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM24_ABELIAN_SUBFABRIC_AND_NON_ABELIAN_RESIDUE",
        "claim_status": "exact_structure_of_the_commutator_sector__no_gap_claim",
        "theorems": {
            "T1_abelian_subfabric_is_circular_euler_orbit":
                "commuting faces compose by F00-E's native addition law; "
                "Stokes additive in the Euler coordinate (rotor chain)",
            "T2_non_abelian_residue_is_gram_defect":
                "|p x q|^2 = |p|^2|q|^2 - <p,q>^2; zero iff one subfabric",
            "T3_second_order_stokes":
                "odd_2 = sum p_i + sum_{i<j} p_i x p_j; even_2 = 1 - sum <p_i,p_j>; "
                "remainder degree >= 3 (exact interpolation)",
            "T4_soft_mode_is_abelian":
                "YM-23's soft mode lies in an exact abelian subfabric; "
                "E_na = |sum_{i<j} p_i x p_j|^2 vanishes there, positive "
                "generically, rotation-invariant",
        },
        "named_next": "YM-25: the direction-change (u-rotation) sector of the "
                      "transfer — the only place a weak-coupling gap can live",
        "witness_unit_direction_u": [str(x) for x in u],
        "controls": {
            "C1_circle_law_and_additive_stokes": bool(c1),
            "C2_lagrange_identity": bool(c2),
            "C3_degree2_truncation_exact": bool(c3),
            "C4_parallel_control": bool(c4),
            "C5_E_na_zero_on_abelian_positive_generic_invariant": bool(c5),
            "C6_order_tamper_is_2_cross": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM24_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM24.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    print("sha256:", sha)
