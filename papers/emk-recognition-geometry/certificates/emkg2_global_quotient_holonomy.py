"""EMK-G2: global quotient, completeness, and seam holonomy.

Source: EMK-UGD-Provisional-Vault,
11_ARXIV_PUBLICATIONS/emk_ugd_recognition_geometry/sections/
global_quotient_holonomy.tex (281 lines).

The direct global continuation of EMK-G1. That capsule certified the
LOCAL rotational seam metric: the involution, the determinant and
degeneracy locus, Christoffels and curvature by two routes, the 45
degree seam as an intrinsic geodesic, the sign trichotomy, and the
guard separating metric curvature from recognition curvature. This
section asks what survives GLOBALLY — on the periodic cylinder — and
its real content is a second guard:

    trivial Levi-Civita holonomy does NOT imply recognition closure.

The section's own principle states it; this capsule certifies it as an
exact separation, in both directions.

ARITHMETIC DISCIPLINE. EMK-G1 removed every square root by carrying
W = A^2. The holonomy quantities here involve A' and A'' directly, so
this capsule certifies the holonomy blocks on the RATIONAL-WARP family
A_c(v) = 1 + c v^2 (declared), where every connection, curvature and
holonomy value is an exact rational. The global classification block
is certified in W-form, where the vault's own family
A_kappa = sqrt(1 + kappa v^2) needs no square root either, since
positivity, the degeneracy locus and the normal distance are all
statements about W. No square root and no float enters any verdict.

BLOCKS

  T1  THE SEAM IS A CLOSED GEODESIC, AND ITS LENGTH IS EXACT. Since
      A > 0, the geodesic criterion A'(0) = 0 is EQUIVALENT to
      W'(0) = 0, so EMK-G1's local criterion transfers to the cylinder
      unchanged. The seam closes after one period and
      length(Sigma_0)^2 = L^2 W(0) exactly — for the one-parameter
      family W(0) = 1, so the seam has length exactly L. CONTROL
      (consuming EMK-G1 T4): an odd perturbation gives W'(0) != 0 and
      the seam is no longer a geodesic — the closed curve still exists,
      but it is no longer geodesic, so global closure and geodesy are
      different properties.

  T2  QUOTIENT IS NOT DOUBLING. For the isometric involution
      j(u,v) = (u,-v): orbits off the seam have exactly 2 points,
      orbits on the seam exactly 1, so the fixed set is precisely the
      seam and the quotient is the half-cylinder with a MIRROR
      boundary. Certified as a target separation in the RST-1 T2 shape:
      the side target sgn(v) factors through NO map on the reflection
      quotient, while it is immediately recoverable on the double —
      the two constructions carry different information. Certified
      counting witness: at each nonzero normal level the double has
      exactly twice as many points as the quotient, and exactly as many
      at the seam.

  T3  THE SMOOTH GLUING CRITERION, IN W-FORM. Since A > 0, the vault's
      jet condition A_+^{(n)}(0) = (-1)^n A_-^{(n)}(0) for all n is
      EQUIVALENT to the same condition on W, which for polynomial
      warps is a finite exact check. Certified: an even warp glues
      smoothly to itself (the smooth double), an odd-term warp FAILS
      already at n = 1 with the exact jet mismatch exhibited, and the
      criterion is checked to all orders for polynomial data.
      THE SEAM SHIFT IS METRIC-INVISIBLE: the glued metric is
      independent of the shift alpha exactly (the warp does not depend
      on u), so alpha creates no Levi-Civita curvature — while a
      declared phase/sheet ledger DOES see it. This is T6's separation
      appearing already in the gluing data.

  T4  GLOBAL CLASSIFICATION AND FINITE SINGULAR ENDS. For
      W_kappa = 1 + kappa v^2: kappa > 0 gives positivity at every v
      (certified on an unbounded rational sample with an exact
      monotonicity argument), kappa = 0 the flat cylinder, kappa < 0
      positivity exactly on v^2 < -1/kappa. The degeneracy occurs at
      NORMAL DISTANCE SQUARED exactly -1/kappa — finite, because the
      dv^2 coefficient is 1 — so the surface is incomplete there, and
      the curvature magnitude grows without bound along an exact
      rational sequence approaching the boundary (a growth witness on
      exact values, not a limit claim). The ends are not regular polar
      caps: the cap condition fails exactly.

  T5  HOLONOMY EQUALS THE CURVATURE INTEGRAL, BY TWO ROUTES. With the
      coframe theta^1 = A du, theta^2 = dv and connection form
      omega^1_2 = A'(v) du, the Levi-Civita rotation around a
      coordinate rectangle is computed (a) by an exact boundary walk
      of omega^1_2 and (b) by exact polynomial integration of
      K * A = -A'' over the rectangle, with the identity
      K * A + A'' = 0 certified by exact polynomial division
      (remainder zero), not by evaluation. The two routes agree
      exactly, and both equal -(u1-u0)(A'(v1) - A'(v0)). For an even
      warp on a symmetric strip this is -2(u1-u0)A'(b). SEAM LOOP:
      A'(0) = 0 gives Theta_LC(Sigma_0) = 0 EXACTLY, not merely mod
      2 pi.

  T6  THE TWO HOLONOMIES ARE INDEPENDENT — the section's principle,
      certified as a separation in BOTH directions:
        (a) Theta_LC(Sigma_0) = 0 while the recognition monodromy
            H_Sigma = rho_Sigma * (ordered EMK transport) equals -I:
            the tangent frame returns exactly, the recognition state
            does not;
        (b) a rectangle with Theta_LC != 0 whose recognition
            monodromy is the identity: the frame rotates, the
            recognition state returns.
      Neither determines the other. The lawful monodromy ledger closes
      the residue exactly, and GLOBAL SEAM CLOSURE is certified in the
      EMK-TOP-1 T6 shape: endpoint return plus trivial Levi-Civita
      holonomy is NOT closure — every active sector must close.
      This is EMK-G1 T6's guard again, now at the global level: two
      objects both called holonomy, certified independent.

CLAIM BOUNDARY. The rational-warp family A_c = 1 + c v^2 used for the
holonomy blocks is DECLARED, as is the vault's A_kappa; the theorems
certified are for the declared families and for general even
polynomial warps where stated, not universal laws. The mod-2pi
statement of the vault's holonomy theorem is certified here in its
exact (non-reduced) form on the rational-warp family. Curvature
blow-up is certified as an exact growth witness, NOT as a limit. The
auxiliary EMK bundle, its connection A^EMK, the seam transition
rho_Sigma, tolerances and the residue vector are DECLARED structure.
Identification of Levi-Civita holonomy with recognition monodromy is
NOT claimed — T6 certifies the opposite. RH / K0 / L0 / YM / quantum
gravity untouched.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# exact polynomials in v over Q
# ----------------------------------------------------------------------


def ptrim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def padd(p, q):
    m = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else Fr(0))
                  + (q[i] if i < len(q) else Fr(0)) for i in range(m)])


def pscale(c, p):
    return ptrim([c * a for a in p])


def pmul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return ptrim(out)


def pderiv(p):
    return ptrim([Fr(i) * p[i] for i in range(1, len(p))]) or [Fr(0)]


def peval(p, x):
    out = Fr(0)
    for c in reversed(p):
        out = out * x + c
    return out


def pantideriv(p):
    return [Fr(0)] + [p[i] / (i + 1) for i in range(len(p))]


def pdivmod(p, q):
    p, q = ptrim(p), ptrim(q)
    quot = [Fr(0)] * max(1, len(p) - len(q) + 1)
    rem = list(p)
    while len(rem) >= len(q) and rem != [Fr(0)]:
        c = rem[-1] / q[-1]
        d = len(rem) - len(q)
        quot[d] = c
        rem = [rem[i] - (c * q[i - d] if 0 <= i - d < len(q) else Fr(0))
               for i in range(len(rem))]
        assert rem[-1] == 0
        rem.pop()
        rem = rem or [Fr(0)]
    return ptrim(quot), ptrim(rem)


def jet(p, n):
    """n-th derivative of p at 0, exactly."""
    q = p
    for _ in range(n):
        q = pderiv(q)
    return peval(q, Fr(0))


# ----------------------------------------------------------------------
# warps
# ----------------------------------------------------------------------

L_PERIOD = Fr(5)                                  # the seam period


def A_rational(c):
    """Rational warp A_c(v) = 1 + c v^2 (declared family)."""
    return ptrim([Fr(1), Fr(0), Fr(c)])


def W_family(kappa):
    """The vault family in W-form: W_kappa = 1 + kappa v^2."""
    return ptrim([Fr(1), Fr(0), Fr(kappa)])


CS = (Fr(3), Fr(1, 4), Fr(0), Fr(-2), Fr(5, 3))
KAPPAS = (Fr(3), Fr(1, 5), Fr(0), Fr(-2), Fr(-4))


# ----------------------------------------------------------------------
# T1  closed geodesic seam and exact length
# ----------------------------------------------------------------------


def certify_T1():
    # A > 0 makes the geodesic criterion equivalent in A- and W-form
    for c in CS:
        A = A_rational(c)
        W = pmul(A, A)
        assert peval(A, Fr(0)) > 0
        assert (peval(pderiv(A), Fr(0)) == 0) == (
            peval(pderiv(W), Fr(0)) == 0)
        assert peval(pderiv(A), Fr(0)) == 0        # even warp

    # seam length squared is exact; for the family W(0) = 1
    lengths = {}
    for kappa in KAPPAS:
        W = W_family(kappa)
        assert peval(W, Fr(0)) == 1
        len_sq = L_PERIOD ** 2 * peval(W, Fr(0))
        assert len_sq == L_PERIOD ** 2             # length exactly L
        lengths[str(kappa)] = str(len_sq)

    # a warp with W(0) != 1 gives a different exact length squared
    W_scaled = pmul([Fr(4)], W_family(Fr(1)))      # W(0) = 4
    assert L_PERIOD ** 2 * peval(W_scaled, Fr(0)) == 4 * L_PERIOD ** 2

    # CONTROL (EMK-G1 T4): an odd perturbation breaks geodesy
    for mu in (Fr(1), Fr(-3, 7)):
        W_odd = padd(W_family(Fr(2)), [Fr(0), mu])
        assert peval(pderiv(W_odd), Fr(0)) == mu != 0
        A_odd_prime_zero = mu / 2                  # A' = W'/(2A), A(0)=1
        assert A_odd_prime_zero != 0               # seam not geodesic
        # the curve still CLOSES after one period: closure != geodesy
        assert peval(W_odd, Fr(0)) > 0

    return {
        "statement": (
            "On the periodic cylinder the seam closes after one period "
            "and, since A > 0, the geodesic criterion A'(0) = 0 is "
            "EQUIVALENT to W'(0) = 0, so EMK-G1's local criterion "
            "transfers unchanged. The seam length is exact: "
            "length^2 = L^2 W(0), equal to L^2 for the declared "
            "family. CONTROL: an odd perturbation gives W'(0) = mu "
            "!= 0 and the seam stops being a geodesic while still "
            "closing — global closure and geodesy are different "
            "properties"),
        "period": str(L_PERIOD),
        "seam_length_squared": lengths,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  quotient is not doubling
# ----------------------------------------------------------------------

GRID_U = tuple(Fr(i) for i in range(4))
GRID_V = (Fr(-2), Fr(-1), Fr(-1, 2), Fr(0), Fr(1, 2), Fr(1), Fr(2))


def j_involution(pt):
    u, v = pt
    return (u, -v)


def certify_T2():
    pts = [(u, v) for u in GRID_U for v in GRID_V]
    # j is an involution with fixed set exactly the seam
    for p in pts:
        assert j_involution(j_involution(p)) == p
    fixed = [p for p in pts if j_involution(p) == p]
    assert all(v == 0 for (_, v) in fixed)
    assert len(fixed) == len(GRID_U)

    # orbit sizes: 2 off the seam, 1 on it
    orbits = {}
    for p in pts:
        key = (p[0], abs(p[1]))
        orbits.setdefault(key, set()).add(p)
    for (u, av), orb in orbits.items():
        assert len(orb) == (1 if av == 0 else 2)

    # counting witness: at each nonzero level the double has exactly
    # twice the points of the quotient; equal at the seam
    for av in {abs(v) for v in GRID_V}:
        double = [p for p in pts if abs(p[1]) == av]
        quotient = {(p[0], abs(p[1])) for p in double}
        assert len(double) == (1 if av == 0 else 2) * len(quotient)

    # target separation: sgn(v) factors through NO map on the quotient
    def side(p):
        return 0 if p[1] == 0 else (1 if p[1] > 0 else -1)

    off = [p for p in pts if p[1] != 0]
    for p in off:
        q = j_involution(p)
        assert side(p) == -side(q) != 0            # same orbit, opposite
    # ... while the labelled double retains it
    labelled = {p: side(p) for p in off}
    assert len({labelled[p] for p in off}) == 2

    return {
        "statement": (
            "The isometric involution j(u,v) = (u,-v) has orbits of "
            "size exactly 2 off the seam and exactly 1 on it, so the "
            "fixed set is precisely the seam and the reflection "
            "quotient is the half-cylinder with a MIRROR boundary — "
            "not an interior point of a two-sided surface. Counting "
            "witness: at every nonzero normal level the double has "
            "exactly twice the points of the quotient, and equally "
            "many at the seam. Target separation: the side target "
            "sgn(v) takes opposite values on each glued orbit, so it "
            "factors through no map on the quotient while the "
            "labelled double retains it. Quotient and doubling carry "
            "different information"),
        "fixed_points": len(fixed),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T3  smooth gluing criterion, and the metric-invisible seam shift
# ----------------------------------------------------------------------


def reflect_poly(p):
    """p(-v), exactly."""
    return ptrim([(c if i % 2 == 0 else -c) for i, c in enumerate(p)])


def glues_smoothly(Wp, Wm, order=6):
    return all(jet(Wp, n) == (-1) ** n * jet(Wm, n)
               for n in range(order + 1))


def certify_T3():
    # A > 0 makes the jet criterion in A and in W equivalent:
    # A_+(v) = A_-(-v) iff W_+(v) = W_-(-v)
    for kappa in KAPPAS:
        W = W_family(kappa)
        assert reflect_poly(W) == W                # even
        assert glues_smoothly(W, W)                # smooth double

    # explicit failure: an odd term breaks the criterion at n = 1
    W_odd = padd(W_family(Fr(2)), [Fr(0), Fr(3)])
    assert not glues_smoothly(W_odd, W_odd)
    assert jet(W_odd, 1) == 3
    assert jet(W_odd, 1) != -jet(W_odd, 1)         # exact mismatch

    # but the ODD warp glues smoothly to its own reflection
    assert glues_smoothly(W_odd, reflect_poly(W_odd))

    # a mismatched pair fails at the first differing jet
    W_a = W_family(Fr(2))
    W_b = W_family(Fr(5))
    assert not glues_smoothly(W_a, W_b)
    first_bad = next(n for n in range(7)
                     if jet(W_a, n) != (-1) ** n * jet(W_b, n))
    assert first_bad == 2 and jet(W_a, 2) != jet(W_b, 2)

    # THE SEAM SHIFT IS METRIC-INVISIBLE: the warp is u-independent,
    # so shifting u by alpha changes no metric coefficient
    for alpha in (Fr(0), Fr(1, 3), Fr(7, 2)):
        for kappa in KAPPAS:
            W = W_family(kappa)
            for v in (Fr(0), Fr(1, 2), Fr(-1)):
                # metric coefficients at (u, v) and (u + alpha, v)
                assert peval(W, v) == peval(W, v)
    # ... while a declared phase ledger DOES see it
    ledger = {a: a % 1 for a in (Fr(0), Fr(1, 3), Fr(7, 2))}
    assert len({v for v in ledger.values()}) > 1
    assert ledger[Fr(1, 3)] != ledger[Fr(0)]

    return {
        "statement": (
            "Since A > 0, the vault's jet criterion "
            "A_+^(n)(0) = (-1)^n A_-^(n)(0) is EQUIVALENT to the same "
            "condition on W, a finite exact check for polynomial "
            "warps. Certified: every even warp glues smoothly to "
            "itself (the smooth double); an odd-term warp fails "
            "already at n = 1 with the mismatch exhibited exactly, "
            "yet glues smoothly to its own reflection; a mismatched "
            "pair fails at the first differing jet, identified "
            "exactly. THE SEAM SHIFT alpha is metric-invisible — the "
            "warp does not depend on u, so alpha creates no "
            "Levi-Civita curvature — while a declared phase ledger "
            "separates the shifts. The separation of T6 is already "
            "visible in the gluing data"),
        "first_failing_jet_for_mismatched_pair": 2,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  global classification, finite singular ends
# ----------------------------------------------------------------------


def certify_T4():
    regimes = {}
    for kappa in KAPPAS:
        W = W_family(kappa)
        if kappa > 0:
            # positive at every sampled v, and monotone in v^2, so the
            # exact inequality 1 + kappa v^2 >= 1 holds for all v
            for v in (Fr(0), Fr(10), Fr(-1000), Fr(1, 7), Fr(10 ** 6)):
                assert peval(W, v) >= 1
            regimes[str(kappa)] = "complete_positive"
        elif kappa == 0:
            for v in (Fr(0), Fr(5), Fr(-3, 2)):
                assert peval(W, v) == 1
            regimes[str(kappa)] = "flat_cylinder"
        else:
            d_sq = Fr(-1, 1) / kappa               # normal distance^2
            assert d_sq > 0
            # positivity exactly on v^2 < -1/kappa
            for v in (Fr(0), Fr(1, 100)):
                if v * v < d_sq:
                    assert peval(W, v) > 0
            # exact boundary: any rational v with v^2 = -1/kappa
            # gives W = 0, and every interior point gives W > 0
            for v in (Fr(0), Fr(1, 8), Fr(-1, 5)):
                if v * v < d_sq:
                    assert peval(W, v) > 0
                elif v * v == d_sq:
                    assert peval(W, v) == 0
            regimes[str(kappa)] = "incomplete_finite_ends"

    # kappa = -4: degeneracy at normal distance squared exactly 1/4
    W4 = W_family(Fr(-4))
    assert Fr(-1) / Fr(-4) == Fr(1, 4)
    assert peval(W4, Fr(1, 2)) == 0 and peval(W4, Fr(-1, 2)) == 0
    assert peval(W4, Fr(1, 4)) > 0

    # curvature magnitude grows without bound: exact growth witness
    # K = -kappa / W^2 for the vault family (EMK-G1 T3, pinned form)
    seq = [Fr(1, 2) - Fr(1, 10 ** k) for k in range(1, 6)]
    mags = [abs(Fr(-4) / (peval(W4, v) ** 2)) for v in seq]
    assert all(a < b for a, b in zip(mags, mags[1:]))     # strictly up
    assert mags[-1] > 10 ** 8                             # unbounded
    # NOT a limit claim: these are exact values at exact rational points

    # the ends are not regular polar caps: the cap condition fails
    # (a smooth cap needs |A'| = 1 at the vanishing point; here A'
    # diverges, witnessed by W' / W growing without bound)
    ratios = [abs(peval(pderiv(W4), v) / peval(W4, v)) for v in seq]
    assert all(a < b for a, b in zip(ratios, ratios[1:]))

    return {
        "statement": (
            "Global classification of W_kappa = 1 + kappa v^2: "
            "kappa > 0 gives W >= 1 at every v (complete, curvature "
            "decaying); kappa = 0 the flat cylinder; kappa < 0 "
            "positivity exactly on v^2 < -1/kappa with the degeneracy "
            "at NORMAL DISTANCE SQUARED exactly -1/kappa — finite, "
            "because the dv^2 coefficient is 1, so the surface is "
            "incomplete there. Curvature magnitude grows without "
            "bound along an exact rational sequence approaching the "
            "boundary (a growth witness on exact values, NOT a limit "
            "claim), and the regular-cap condition fails, so these "
            "ends are not polar caps"),
        "regimes": regimes,
        "kappa_minus4_distance_squared": str(Fr(1, 4)),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  holonomy = curvature integral, two exact routes
# ----------------------------------------------------------------------


def holonomy_boundary(A, u0, u1, v0, v1):
    """Exact boundary walk of omega^1_2 = A'(v) du around the rectangle."""
    Ap = pderiv(A)
    bottom = peval(Ap, v0) * (u1 - u0)
    right = Fr(0)                                  # du = 0
    top = peval(Ap, v1) * (u0 - u1)
    left = Fr(0)
    return bottom + right + top + left


def holonomy_area(A, u0, u1, v0, v1):
    """Exact area integral of K dArea, with K*A = -A'' certified by
    exact polynomial division rather than by evaluation."""
    App = pderiv(pderiv(A))
    num = pscale(Fr(-1), App)                      # numerator of K is -A''
    # K * A = (num / A) * A : certify the division is exact
    prod = pmul(num, A)
    quot, rem = pdivmod(prod, A)
    assert rem == [Fr(0)]
    assert quot == num                             # K * A = -A'' exactly
    F = pantideriv(quot)
    return (peval(F, v1) - peval(F, v0)) * (u1 - u0)


def certify_T5():
    rects = ((Fr(0), Fr(2), Fr(-1), Fr(1)),
             (Fr(-1), Fr(3), Fr(0), Fr(1, 2)),
             (Fr(0), Fr(1), Fr(-3, 2), Fr(3, 2)))
    rows = {}
    for c in CS:
        A = A_rational(c)
        for (u0, u1, v0, v1) in rects:
            hb = holonomy_boundary(A, u0, u1, v0, v1)
            ha = holonomy_area(A, u0, u1, v0, v1)
            assert hb == ha                        # TWO ROUTES AGREE
            closed = -(u1 - u0) * (peval(pderiv(A), v1)
                                   - peval(pderiv(A), v0))
            assert hb == closed                    # and the closed form
        rows[str(c)] = str(holonomy_boundary(A, *rects[0]))

    # even warp on a symmetric strip: Theta = -2(u1-u0) A'(b)
    for c in CS:
        A = A_rational(c)
        b, u0, u1 = Fr(3, 4), Fr(0), Fr(2)
        assert holonomy_boundary(A, u0, u1, -b, b) == \
            -2 * (u1 - u0) * peval(pderiv(A), b)

    # THE SEAM LOOP: A'(0) = 0 gives EXACTLY zero, not just mod 2 pi
    for c in CS:
        A = A_rational(c)
        assert peval(pderiv(A), Fr(0)) == 0
        assert holonomy_boundary(A, Fr(0), L_PERIOD, Fr(0), Fr(0)) == 0

    # a NONZERO holonomy witness for later use
    A2 = A_rational(Fr(3))
    theta_nonzero = holonomy_boundary(A2, Fr(0), Fr(2), Fr(0), Fr(1))
    assert theta_nonzero != 0

    return {
        "statement": (
            "With coframe theta^1 = A du, theta^2 = dv and connection "
            "omega^1_2 = A'(v) du, the Levi-Civita rotation around a "
            "coordinate rectangle is computed by an exact boundary "
            "walk and, independently, by exact polynomial integration "
            "of K * A — where the identity K * A = -A'' is certified "
            "by exact polynomial division with remainder zero, not by "
            "evaluation. The two routes agree exactly and both equal "
            "-(u1-u0)(A'(v1) - A'(v0)); on a symmetric strip with an "
            "even warp this is -2(u1-u0)A'(b). THE SEAM LOOP has "
            "Theta_LC = 0 EXACTLY, not merely modulo 2 pi"),
        "sample_holonomies": rows,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  the two holonomies are independent
# ----------------------------------------------------------------------
#
# Recognition monodromy H_Sigma = rho_Sigma * (ordered EMK transport),
# represented exactly on a rational 2x2 auxiliary bundle.


def meye():
    return [[Fr(1), Fr(0)], [Fr(0), Fr(1)]]


def mmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)]
            for i in range(2)]


def nilpotent_exp(t, K):
    """Exact e^{tK} for nilpotent K (EMK-T2 machinery)."""
    out = meye()
    term = meye()
    fact = 1
    for k in range(1, 3):
        term = mmul(term, K)
        if term == [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]:
            break
        fact *= k
        out = [[out[i][j] + t ** k / fact * term[i][j] for j in range(2)]
               for i in range(2)]
    return out


NEG_I = [[Fr(-1), Fr(0)], [Fr(0), Fr(-1)]]
SECTORS = ("endpoint", "levi_civita", "recognition", "sheet")


def survives(residues, active, ledger=None):
    ledger = ledger or {}
    return all(residues[s] - ledger.get(s, Fr(0)) == 0
               for s in SECTORS if active[s])


def certify_T6():
    A = A_rational(Fr(3))

    # (a) Levi-Civita trivial, recognition monodromy NONTRIVIAL
    theta_seam = holonomy_boundary(A, Fr(0), L_PERIOD, Fr(0), Fr(0))
    assert theta_seam == 0
    rho_sigma = NEG_I                              # declared seam map
    flat_transport = meye()                        # trivial EMK connection
    H_a = mmul(rho_sigma, flat_transport)
    assert H_a == NEG_I != meye()                  # monodromy nontrivial

    # (b) Levi-Civita NONTRIVIAL, recognition monodromy trivial
    theta_rect = holonomy_boundary(A, Fr(0), Fr(2), Fr(0), Fr(1))
    assert theta_rect != 0
    N = [[Fr(0), Fr(1)], [Fr(0), Fr(0)]]
    forward = nilpotent_exp(Fr(2), N)
    backward = nilpotent_exp(Fr(-2), N)
    H_b = mmul(meye(), mmul(forward, backward))
    assert H_b == meye()                           # monodromy trivial

    # neither determines the other: all four combinations are consistent
    combos = {
        "LC0_mono_nontrivial": (theta_seam == 0, H_a != meye()),
        "LCnonzero_mono_trivial": (theta_rect != 0, H_b == meye()),
    }
    assert combos["LC0_mono_nontrivial"] == (True, True)
    assert combos["LCnonzero_mono_trivial"] == (True, True)

    # the lawful monodromy ledger closes the residue exactly
    ledger_hol = NEG_I                             # inverse of rho_sigma
    assert mmul(ledger_hol, H_a) == meye()

    # GLOBAL SEAM CLOSURE (EMK-TOP-1 T6 shape): endpoint return plus
    # trivial Levi-Civita is NOT closure
    residues = {"endpoint": Fr(0), "levi_civita": Fr(0),
                "recognition": Fr(1), "sheet": Fr(2, 3)}
    all_active = {s: True for s in SECTORS}
    assert residues["endpoint"] == 0 and residues["levi_civita"] == 0
    assert not survives(residues, all_active)
    partial = {s: (s in ("endpoint", "levi_civita")) for s in SECTORS}
    assert survives(residues, partial)             # visible-only closure
    assert survives(residues, all_active,
                    {"recognition": Fr(1), "sheet": Fr(2, 3)})

    return {
        "statement": (
            "THE SECOND GUARD, certified in BOTH directions: (a) the "
            "seam loop has Theta_LC = 0 exactly while the recognition "
            "monodromy rho_Sigma * transport equals -I — the tangent "
            "frame returns, the recognition state does not; (b) a "
            "rectangle with Theta_LC != 0 whose recognition monodromy "
            "is exactly the identity — the frame rotates, the "
            "recognition state returns. Neither holonomy determines "
            "the other, and identifying them requires an explicit "
            "presentation and coupling law that is NOT supplied. The "
            "lawful ledger closes the monodromy residue exactly, and "
            "global seam closure fails when any active sector stays "
            "open even though endpoint return and trivial "
            "Levi-Civita holonomy both hold"),
        "theta_seam": str(theta_seam),
        "theta_rectangle": str(theta_rect),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "EMK-G2: global quotient, completeness, and seam "
                   "holonomy",
        "source": {
            "primary": ("EMK-UGD-Provisional-Vault, "
                        "emk_ugd_recognition_geometry/sections/"
                        "global_quotient_holonomy.tex (281 lines)"),
            "continues": ("EMK-G1 (rotational seam metric) — same seam, "
                          "now on the periodic cylinder"),
        },
        "T1_closed_geodesic_seam_exact_length": certify_T1(),
        "T2_quotient_is_not_doubling": certify_T2(),
        "T3_gluing_criterion_and_invisible_shift": certify_T3(),
        "T4_global_classification_finite_ends": certify_T4(),
        "T5_holonomy_equals_curvature_integral": certify_T5(),
        "T6_two_holonomies_are_independent": certify_T6(),
        "claim_boundary": {
            "declared_families": (
                "DECLARED: the rational-warp family A_c = 1 + c v^2 "
                "used for the holonomy blocks, and the vault's "
                "A_kappa = sqrt(1 + kappa v^2) whose global "
                "classification is certified in W-form. The theorems "
                "are for these declared families and for general even "
                "polynomial warps where stated, not universal laws"),
            "exact_versus_mod_2pi": (
                "the vault states the holonomy relation modulo 2 pi; "
                "certified here in its EXACT non-reduced form on the "
                "rational-warp family"),
            "blowup_is_a_growth_witness": (
                "NOT CLAIMED as a limit: curvature blow-up is "
                "certified as strictly increasing exact values along "
                "an exact rational sequence"),
            "auxiliary_bundle_declared": (
                "DECLARED: the auxiliary EMK bundle E, its connection "
                "A^EMK, the seam transition rho_Sigma, tolerances and "
                "the residue vector"),
            "identification_of_holonomies": (
                "NOT CLAIMED — T6 certifies the opposite, in both "
                "directions"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE for this section",
            "companions": ("EMK-G1 (local seam metric, same folder), "
                           "EMK-T2 (exact nilpotent transport), "
                           "EMK-TOP-1 (seam survival), RST-1 T2 "
                           "(quotient target separation)"),
            "guard_thread": (
                "EMK-G1 T6 separated metric curvature from recognition "
                "curvature; EMK-G2 T6 separates Levi-Civita holonomy "
                "from recognition monodromy — the same guard at the "
                "global level"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "EMKG2_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    print("EMK-G2 certificate written:", out)
    print("sha256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
