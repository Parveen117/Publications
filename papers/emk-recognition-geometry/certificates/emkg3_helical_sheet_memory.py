"""EMK-G3: the helical sheet-memory lift — flat is not memoryless.

Source: EMK-UGD-Provisional-Vault,
11_ARXIV_PUBLICATIONS/emk_ugd_recognition_geometry/sections/
helical_sheet_memory.tex (245 lines).

The door EMK-G2 left open. That capsule's closing remark says helical
lifts can now be added without confusing them with the base-surface
metric; this section adds them. Its central principle is sharper than
anything in the two preceding geometry capsules:

    A FLAT lifted connection may carry NONTRIVIAL global monodromy.
    Local curvature detects infinitesimal nonclosure; monodromy
    detects global return data. Conflating them erases exactly the
    memory sector the lift exists to carry.

The section is also careful about what the helix IS: not a curve drawn
in three-space, but the orbit structure of a deck transformation — a
mapping torus. This capsule certifies it that way, so no picture does
any work.

ARITHMETIC DISCIPLINE. The fibre is the minimal abelian model
F = R_phi x R_sigma x Z_k with the transition acting by TRANSLATION,
so every monodromy is exact rational/integer arithmetic and no
exponential is ever formed. Phase is carried as an EXPONENT modulo K
(UGD-1's discipline), so "visible phase return" is the exact condition
alpha = 0 mod K rather than a statement about 2 pi. Levi-Civita data
is consumed from the PINNED EMK-G2 module, not recomputed.

BLOCKS

  T1  THE MAPPING TORUS IS A BUNDLE, NOT A PICTURE. The deck
      transformation (u,v,f) -> (u+L, v, rho^{-1} f) generates a free
      action: certified exactly, no orbit point is fixed by a nonzero
      power, and every orbit meets the fundamental domain 0 <= u < L
      exactly once. The transition rho is a bijection with an exact
      inverse, and the monodromy after n circuits is rho^n EXACTLY,
      for positive, zero and negative n — so the bundle is well
      defined and its monodromy is the declared automorphism.

  T2  PITCH AND COMPENSATED COORDINATES. The compensated variables
      phi_hat = phi - p_phi u and sigma_hat = sigma - p_sigma u, with
      pitches p_phi = alpha/L and p_sigma = beta/L, are EXACTLY
      invariant under the deck transformation on a rational grid,
      while the uncompensated phi and sigma are exactly NOT — the
      failure equals alpha and beta respectively. The pitch is fibre
      advance per unit visible seam length, certified as an exact
      rational.

  T3  FLAT IS NOT MEMORYLESS — the principle, as an exact separation.
      For the constant-pitch connection A = p du the abelian curvature
      d A = (du b - dv a) du ^ dv is EXACTLY ZERO as a bivariate
      polynomial identity, and every contractible rectangle loop has
      EXACTLY ZERO holonomy — yet one circuit of the noncontractible
      base circle gives holonomy EXACTLY alpha = p L, nonzero. Local
      flatness and global monodromy are certified independent: also
      exhibited is a non-flat connection whose curvature is exactly
      -1 with nonzero rectangle holonomy, so the two quantities move
      independently in both directions.

  T4  VISIBLE RETURN WITH LIFTED NONRETURN. With phase an exponent mod
      K, visible return is exactly alpha = 0 mod K, while full lifted
      return additionally requires beta = 0 and q = 0. Certified
      exhaustively on a parameter grid: every combination with
      alpha = 0 mod K and (beta, q) != (0,0) has visible return and
      lifted nonreturn, and lifted return holds EXACTLY when all three
      close. After n circuits the fibre displacement is exactly
      (n alpha, n beta, n q), so a nonzero sheet increment can never
      be undone by more circuits unless n = 0.

  T5  THE FOUR RETURN CLASSES, REALIZED. The section's corollary
      names four classes; a classifier decides them by exact
      arithmetic with a witness of each — exact return (all residues
      zero), lawfully transported return (nonzero advances equal
      their declared ledgers), memory-bearing return (visible sectors
      closed, an admitted sheet entry nonzero and FULLY ledgered), and
      open obstruction (an active unledgered entry). The four are
      certified mutually exclusive and exhaustive on the witness set,
      and the memory-bearing class is certified DISTINCT from exact
      return: the ledger is what makes it lawful, and the sheet entry
      is still there.

  T6  NO PROPER SUBSET OF SECTORS FORCES CLOSURE. Consuming the
      PINNED EMK-G2 seam holonomy: base-point return holds, the
      Levi-Civita seam holonomy is exactly zero, and the lifted
      curvature is exactly flat — and the transition is STILL open
      because the helical residue is nonzero. Certified in the
      strongest available form: for EVERY proper subset of the five
      declared sectors there is a residue vector closing exactly that
      subset and no more, so no proper subset can force closure.
      Also certified: a structure-preserving projection maps closure
      to closure, while a NON-FAITHFUL projection that drops the
      sheet sector maps an OPEN lift to a CLOSED image — the
      converse of the projection theorem fails, exactly.

CLAIM BOUNDARY. The fibre model, the transition rho, the pitches, the
EMK connection A^EMK, the ledgers, the tolerances and the thermodynamic
presentation map are DECLARED structure. The non-abelian lift, the
distributional/cut-localized sheet term, and any continuous relaxation
of the sheet sector are NOT certified. The identification of alpha,
beta, q with phase lag, scale drift or hysteretic branch count is NOT
claimed — the section itself says the mapping torus provides
bookkeeping geometry and does not manufacture a physical law. Phase
"visible return" is certified in the exponent-mod-K form, not as a
statement about 2 pi. RH / K0 / L0 / YM / quantum gravity untouched.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g2 = _load("emkg2_for_g3", "emkg2_global_quotient_holonomy.py")

# ----------------------------------------------------------------------
# declared model constants
# ----------------------------------------------------------------------

L = g2.L_PERIOD                      # the base period, from EMK-G2
KMOD = 6                             # phase exponents live in Z/KMOD

ALPHA = 3                            # phase advance per period (exponent)
BETA = Fr(2, 5)                      # support/scale transport
Q = 2                                # sheet-memory increment


# ----------------------------------------------------------------------
# T1  the mapping torus
# ----------------------------------------------------------------------


def rho(f, alpha=ALPHA, beta=BETA, q=Q):
    """The elementary sheet transition, acting by translation."""
    phi, sigma, k = f
    return ((phi + alpha) % KMOD, sigma + beta, k + q)


def rho_inv(f, alpha=ALPHA, beta=BETA, q=Q):
    phi, sigma, k = f
    return ((phi - alpha) % KMOD, sigma - beta, k - q)


def rho_power(f, n, alpha=ALPHA, beta=BETA, q=Q):
    phi, sigma, k = f
    return ((phi + n * alpha) % KMOD, sigma + n * beta, k + n * q)


def deck(pt, n):
    """(u, v, f) -> (u + nL, v, rho^{-n} f)."""
    u, v, f = pt
    return (u + n * L, v, rho_power(f, -n))


FIBRES = tuple((p, Fr(s), k)
               for p in (0, 2, 5) for s in (0, 1, -3) for k in (0, 1, -2))


def certify_T1():
    # rho is a bijection with an exact inverse
    for f in FIBRES:
        assert rho_inv(rho(f)) == f
        assert rho(rho_inv(f)) == f
    images = {rho(f) for f in FIBRES}
    assert len(images) == len(FIBRES)                # injective

    # monodromy after n circuits is exactly rho^n
    for f in FIBRES:
        acc = f
        for n in range(1, 5):
            acc = rho(acc)
            assert acc == rho_power(f, n)
        acc = f
        for n in range(1, 4):
            acc = rho_inv(acc)
            assert acc == rho_power(f, -n)
        assert rho_power(f, 0) == f

    # the deck action is FREE: no point is fixed by a nonzero power
    pts = [(Fr(u), Fr(v), f) for u in (0, 1) for v in (0, 2)
           for f in FIBRES[:6]]
    for p in pts:
        for n in range(-3, 4):
            if n != 0:
                assert deck(p, n) != p               # free
                assert deck(p, n)[0] != p[0]         # u moves by nL

    # every orbit meets the fundamental domain 0 <= u < L exactly once
    for p in pts:
        hits = [n for n in range(-4, 5)
                if Fr(0) <= deck(p, n)[0] < L]
        assert len(hits) == 1

    # composition: deck(deck(p, m), n) = deck(p, m + n)
    for p in pts[:4]:
        for m in (-2, 0, 3):
            for n in (-1, 1, 2):
                assert deck(deck(p, m), n) == deck(p, m + n)

    return {
        "statement": (
            "The helix is certified as the orbit structure of a deck "
            "transformation, not as a drawn curve. The transition rho "
            "is a bijection with an exact inverse; the monodromy after "
            "n circuits is EXACTLY rho^n for positive, zero and "
            "negative n; the deck action is FREE (no point fixed by "
            "any nonzero power) and composes additively; and every "
            "orbit meets the fundamental domain 0 <= u < L exactly "
            "once. The mapping torus is therefore a well-defined "
            "bundle whose monodromy is the declared automorphism"),
        "period": str(L),
        "transition": {"alpha": ALPHA, "beta": str(BETA), "q": Q,
                       "K": KMOD},
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  pitch and compensated coordinates
# ----------------------------------------------------------------------

P_PHI = Fr(ALPHA) / L
P_SIGMA = BETA / L


def compensated(u, phi, sigma):
    return (phi - P_PHI * u, sigma - P_SIGMA * u)


def certify_T2():
    # exact invariance under the deck transformation, on a rational grid
    for u in (Fr(0), Fr(1, 3), Fr(-2), Fr(7, 2)):
        for phi in (Fr(0), Fr(5, 4), Fr(-1)):
            for sigma in (Fr(0), Fr(3), Fr(-2, 7)):
                before = compensated(u, phi, sigma)
                after = compensated(u + L, phi + ALPHA, sigma + BETA)
                assert after == before               # invariant, exactly

    # the UNCOMPENSATED coordinates are exactly not invariant, and the
    # failure is exactly alpha and beta
    for u in (Fr(0), Fr(1, 3)):
        phi, sigma = Fr(2), Fr(1, 2)
        assert (phi + ALPHA) - phi == ALPHA != 0
        assert (sigma + BETA) - sigma == BETA != 0

    # pitch is fibre advance per unit visible seam length
    assert P_PHI * L == ALPHA
    assert P_SIGMA * L == BETA
    assert P_PHI == Fr(3, 5) and P_SIGMA == Fr(2, 25)

    # a wrong pitch fails to compensate, by an exact amount
    wrong = P_PHI + Fr(1, 4)
    u, phi = Fr(0), Fr(0)
    lhs = (phi + ALPHA) - wrong * (u + L)
    rhs = phi - wrong * u
    assert lhs - rhs == ALPHA - wrong * L == -Fr(5, 4)

    return {
        "statement": (
            "The compensated variables phi - p_phi u and "
            "sigma - p_sigma u are EXACTLY invariant under the deck "
            "transformation on a rational grid, while the "
            "uncompensated coordinates fail by exactly alpha and "
            "beta. The pitches p_phi = alpha/L and p_sigma = beta/L "
            "are exact rationals measuring fibre advance per unit "
            "visible seam length, and a wrong pitch fails to "
            "compensate by an exactly computed amount"),
        "p_phi": str(P_PHI),
        "p_sigma": str(P_SIGMA),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# bivariate polynomials in (u, v) over Q, for the flatness block
# ----------------------------------------------------------------------


def bp(terms):
    """{(i, j): c} with c the coefficient of u^i v^j; zeros dropped."""
    return {k: v for k, v in terms.items() if v != 0}


def bp_du(p):
    return bp({(i - 1, j): c * i for (i, j), c in p.items() if i > 0})


def bp_dv(p):
    return bp({(i, j - 1): c * j for (i, j), c in p.items() if j > 0})


def bp_sub(p, q):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, Fr(0)) - c
    return bp(out)


def bp_eval(p, u, v):
    return sum(c * u ** i * v ** j for (i, j), c in p.items())


def bp_int_du(p, u0, u1, v):
    """Exact integral of p du at fixed v."""
    tot = Fr(0)
    for (i, j), c in p.items():
        tot += c * v ** j * (u1 ** (i + 1) - u0 ** (i + 1)) / (i + 1)
    return tot


def bp_int_dv(p, v0, v1, u):
    tot = Fr(0)
    for (i, j), c in p.items():
        tot += c * u ** i * (v1 ** (j + 1) - v0 ** (j + 1)) / (j + 1)
    return tot


def abelian_curvature(a, b):
    """F = d(a du + b dv) = (du b - dv a) du ^ dv."""
    return bp_sub(bp_du(b), bp_dv(a))


def loop_holonomy(a, b, u0, u1, v0, v1):
    """Exact circulation of a du + b dv around a coordinate rectangle."""
    return (bp_int_du(a, u0, u1, v0)
            + bp_int_dv(b, v0, v1, u1)
            + bp_int_du(a, u1, u0, v1)
            + bp_int_dv(b, v1, v0, u0))


# ----------------------------------------------------------------------
# T3  flat is not memoryless
# ----------------------------------------------------------------------


def certify_T3():
    # the constant-pitch connection A = p du
    a_flat = bp({(0, 0): P_PHI})
    b_flat = bp({})
    F = abelian_curvature(a_flat, b_flat)
    assert F == {}                                   # EXACTLY flat

    # every contractible rectangle loop has exactly zero holonomy
    rects = ((Fr(0), Fr(2), Fr(-1), Fr(1)),
             (Fr(-1), Fr(3), Fr(0), Fr(1, 2)),
             (Fr(1, 2), Fr(5, 2), Fr(-2), Fr(3)))
    for r in rects:
        assert loop_holonomy(a_flat, b_flat, *r) == 0

    # ONE CIRCUIT OF THE NONCONTRACTIBLE CIRCLE: exactly alpha
    for v in (Fr(0), Fr(1), Fr(-3, 4)):
        circuit = bp_int_du(a_flat, Fr(0), L, v)
        assert circuit == P_PHI * L == ALPHA != 0
    # n circuits give exactly n alpha
    for n in (2, 5, -3):
        assert bp_int_du(a_flat, Fr(0), n * L, Fr(0)) == n * ALPHA

    # the OTHER direction: a non-flat connection with nonzero curvature
    a_curved = bp({(0, 0): P_PHI, (0, 1): Fr(1)})    # p + v
    F2 = abelian_curvature(a_curved, b_flat)
    assert F2 == {(0, 0): Fr(-1)}                    # curvature exactly -1
    hol = loop_holonomy(a_curved, b_flat, *rects[0])
    assert hol != 0
    # and it equals the exact curvature flux over the rectangle
    u0, u1, v0, v1 = rects[0]
    flux = Fr(-1) * (u1 - u0) * (v1 - v0)
    assert hol == flux

    # so the two quantities move independently in both directions
    assert (F == {}) and (bp_int_du(a_flat, Fr(0), L, Fr(0)) != 0)
    assert (F2 != {}) and (loop_holonomy(a_curved, b_flat,
                                         *rects[0]) != 0)

    return {
        "statement": (
            "THE PRINCIPLE, as an exact separation. For the "
            "constant-pitch connection A = p du the abelian curvature "
            "is EXACTLY ZERO as a bivariate polynomial identity and "
            "every contractible rectangle loop has EXACTLY ZERO "
            "holonomy — yet one circuit of the noncontractible base "
            "circle gives holonomy exactly alpha = p L != 0, and n "
            "circuits give exactly n alpha. A non-flat connection with "
            "curvature exactly -1 and nonzero rectangle holonomy "
            "matching its exact flux completes the separation: local "
            "curvature and global monodromy are independent. Flat does "
            "NOT mean memoryless"),
        "flat_curvature": "0",
        "flat_circuit_holonomy": str(ALPHA),
        "curved_curvature": "-1",
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  visible return with lifted nonreturn
# ----------------------------------------------------------------------


def visible_returns(alpha, n=1):
    """Visible phase return: the exponent closes modulo K."""
    return (n * alpha) % KMOD == 0


def lifted_returns(alpha, beta, q, n=1):
    return (visible_returns(alpha, n) and n * beta == 0 and n * q == 0)


def certify_T4():
    alphas = (0, 2, 3, 6, KMOD)
    betas = (Fr(0), Fr(2, 5), Fr(-1))
    qs = (0, 1, -2)

    separating = 0
    for al in alphas:
        for be in betas:
            for q in qs:
                vis = visible_returns(al)
                lift = lifted_returns(al, be, q)
                # lifted return implies visible return, always
                assert (not lift) or vis
                # and lifted return holds exactly when all three close
                assert lift == (al % KMOD == 0 and be == 0 and q == 0)
                if vis and not lift:
                    separating += 1
    assert separating > 0

    # the headline witness: phase closes, sheet and support do not
    assert visible_returns(KMOD) and not lifted_returns(KMOD, BETA, Q)
    assert visible_returns(0) and not lifted_returns(0, Fr(0), Q)
    assert lifted_returns(0, Fr(0), 0)

    # after n circuits the displacement is exactly (n a, n b, n q)
    f0 = (0, Fr(0), 0)
    for n in (1, 3, 7, -2):
        f = rho_power(f0, n)
        assert f == ((n * ALPHA) % KMOD, n * BETA, n * Q)
    # a nonzero sheet increment is never undone by more circuits
    for n in range(-6, 7):
        if n != 0:
            assert rho_power(f0, n)[2] == n * Q != 0

    return {
        "statement": (
            "With phase carried as an exponent modulo K, visible "
            "return is exactly alpha = 0 mod K while full lifted "
            "return additionally requires beta = 0 and q = 0 — "
            "certified exhaustively on a parameter grid, with lifted "
            "return always implying visible return and never the "
            "converse. The headline witness: the phase closes while "
            "support and sheet do not. After n circuits the fibre "
            "displacement is exactly (n alpha, n beta, n q), so a "
            "nonzero sheet increment can never be undone by further "
            "circuits unless n = 0: sheet memory is not a phase that "
            "eventually wraps"),
        "grid_separating_cases": separating,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  the four return classes
# ----------------------------------------------------------------------


def classify_return(res, ledger, admitted):
    """res = (phase, support, sheet, emk) residues before ledgering."""
    phase, support, sheet, emk = res
    lp, ls, lk, le = ledger
    net = (phase - lp, support - ls, sheet - lk, emk - le)
    if all(x == 0 for x in res):
        return "exact_return"
    if all(x == 0 for x in net):
        # everything nonzero was matched by a declared ledger
        if sheet != 0 and admitted["sheet"]:
            return "memory_bearing_return"
        return "lawfully_transported_return"
    return "open_obstruction"


def certify_T5():
    admitted = {"sheet": True}
    witnesses = {
        "exact_return": ((0, Fr(0), 0, Fr(0)), (0, Fr(0), 0, Fr(0))),
        "lawfully_transported_return": (
            (2, Fr(3, 4), 0, Fr(0)), (2, Fr(3, 4), 0, Fr(0))),
        "memory_bearing_return": (
            (0, Fr(0), 3, Fr(0)), (0, Fr(0), 3, Fr(0))),
        "open_obstruction": (
            (0, Fr(0), 3, Fr(1, 2)), (0, Fr(0), 3, Fr(0))),
    }
    for expected, (res, led) in witnesses.items():
        assert classify_return(res, led, admitted) == expected

    # the four classes are mutually exclusive on the witness set
    labels = [classify_return(r, l, admitted)
              for (r, l) in witnesses.values()]
    assert len(set(labels)) == 4

    # memory-bearing is DISTINCT from exact return: the sheet entry is
    # still there, the ledger is what makes it lawful
    res_mem, led_mem = witnesses["memory_bearing_return"]
    assert res_mem[2] != 0                       # memory present
    assert classify_return(res_mem, led_mem, admitted) != "exact_return"
    # remove the ledger and the same state becomes an open obstruction
    assert classify_return(res_mem, (0, Fr(0), 0, Fr(0)),
                           admitted) == "open_obstruction"
    # de-admit the sheet sector and it reads as lawful transport
    assert classify_return(res_mem, led_mem,
                           {"sheet": False}) == \
        "lawfully_transported_return"

    return {
        "statement": (
            "The section's four return classes are decided by exact "
            "arithmetic with a witness of each, mutually exclusive on "
            "the witness set. MEMORY-BEARING RETURN is certified "
            "distinct from exact return: the sheet entry is still "
            "nonzero and the declared ledger is precisely what makes "
            "it lawful — removing the ledger turns the same state into "
            "an open obstruction, and de-admitting the sheet sector "
            "reclassifies it as ordinary lawful transport. The class "
            "depends on the declared admissions, not on the visible "
            "state alone"),
        "classes": sorted(witnesses),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  no proper subset of sectors forces closure
# ----------------------------------------------------------------------

SECTORS = ("base_point", "levi_civita", "helical", "rtc_curvature",
           "thermodynamic")


def closed_on(residues, active):
    return all(residues[s] == 0 for s in SECTORS if active[s])


def certify_T6():
    # consume the PINNED EMK-G2 Levi-Civita seam holonomy
    A = g2.A_rational(Fr(3))
    theta_seam = g2.holonomy_boundary(A, Fr(0), L, Fr(0), Fr(0))
    assert theta_seam == 0

    # the lifted connection is exactly flat (T3) ...
    F = abelian_curvature(bp({(0, 0): P_PHI}), bp({}))
    assert F == {}

    # ... base point returns, and yet the helical residue is open
    residues = {"base_point": Fr(0), "levi_civita": Fr(theta_seam),
                "helical": Fr(Q), "rtc_curvature": Fr(0),
                "thermodynamic": Fr(0)}
    all_active = {s: True for s in SECTORS}
    assert residues["base_point"] == 0
    assert residues["levi_civita"] == 0
    assert residues["rtc_curvature"] == 0
    assert not closed_on(residues, all_active)

    # THE STRONG FORM: for EVERY proper subset of sectors there is a
    # residue vector closing exactly that subset and no more
    proper = 0
    for r in range(len(SECTORS)):
        for S in combinations(SECTORS, r):
            res = {s: (Fr(0) if s in S else Fr(1)) for s in SECTORS}
            active_S = {s: (s in S) for s in SECTORS}
            assert closed_on(res, active_S)          # that subset closes
            assert not closed_on(res, all_active)    # the whole does not
            proper += 1
    assert proper == 2 ** len(SECTORS) - 1           # all proper subsets

    # PROJECTION: structure-preserving sends closure to closure
    def project(res, faithful):
        out = dict(res)
        if not faithful:
            out["helical"] = Fr(0)                   # drops sheet memory
        return out

    closed_res = {s: Fr(0) for s in SECTORS}
    assert closed_on(project(closed_res, True), all_active)
    assert closed_on(project(closed_res, False), all_active)

    # ... but a NON-FAITHFUL projection maps an OPEN lift to a CLOSED
    # image: the converse of the projection theorem fails, exactly
    assert not closed_on(residues, all_active)
    assert closed_on(project(residues, False), all_active)
    assert not closed_on(project(residues, True), all_active)

    return {
        "statement": (
            "Consuming the pinned EMK-G2 seam holonomy: base-point "
            "return holds, Levi-Civita holonomy is exactly zero, RTC "
            "curvature is zero and the lifted connection is exactly "
            "flat — and the transition is STILL open, because the "
            "helical residue is nonzero. Certified in the strongest "
            "available form: for EVERY one of the 31 proper subsets "
            "of the five declared sectors there is a residue vector "
            "closing exactly that subset while the whole stays open, "
            "so no proper subset can force closure. And a "
            "structure-preserving projection maps closure to closure "
            "while a NON-FAITHFUL projection dropping the sheet "
            "sector maps an OPEN lift to a CLOSED image — the "
            "converse of the projection theorem fails exactly"),
        "proper_subsets_checked": 2 ** len(SECTORS) - 1,
        "sectors": list(SECTORS),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "EMK-G3: helical sheet-memory lift — flat is not "
                   "memoryless",
        "source": {
            "primary": ("EMK-UGD-Provisional-Vault, "
                        "emk_ugd_recognition_geometry/sections/"
                        "helical_sheet_memory.tex (245 lines)"),
            "continues": ("EMK-G2, whose closing remark states that "
                          "helical lifts can now be added without "
                          "confusing them with the base-surface "
                          "metric; the Levi-Civita seam holonomy is "
                          "consumed from that pinned module"),
        },
        "T1_mapping_torus_is_a_bundle": certify_T1(),
        "T2_pitch_and_compensated_coordinates": certify_T2(),
        "T3_flat_is_not_memoryless": certify_T3(),
        "T4_visible_return_lifted_nonreturn": certify_T4(),
        "T5_four_return_classes": certify_T5(),
        "T6_no_proper_subset_forces_closure": certify_T6(),
        "claim_boundary": {
            "declared_structure": (
                "DECLARED: the fibre model F = R_phi x R_sigma x Z_k, "
                "the transition rho, the pitches, the EMK connection "
                "A^EMK, the ledgers, the tolerances and the "
                "thermodynamic presentation map"),
            "not_certified": (
                "NOT CERTIFIED: the non-abelian lift, the "
                "cut-localized or distributional sheet term, and any "
                "continuous relaxation of the sheet sector"),
            "physical_identification": (
                "NOT CLAIMED: that alpha, beta, q are phase lag, "
                "scale drift or hysteretic branch count — the section "
                "itself says the mapping torus supplies bookkeeping "
                "geometry and does not manufacture a physical law"),
            "phase_convention": (
                "visible phase return is certified in the "
                "exponent-mod-K form (UGD-1 discipline), NOT as a "
                "statement about 2 pi; no root of unity is ever "
                "formed"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE for this section",
            "companions": ("EMK-G1 (local seam metric), EMK-G2 "
                           "(global quotient and holonomy, PINNED), "
                           "EMK-TOP-1 (seam survival), UGD-1 (phase "
                           "exponents), EMK-T2 (ordered transport)"),
            "guard_thread": (
                "EMK-G1 T6 metric vs recognition curvature; EMK-G2 T6 "
                "Levi-Civita holonomy vs recognition monodromy; "
                "EMK-G3 T3 local flatness vs global monodromy — the "
                "same guard three times, each at a higher level"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "EMKG3_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    print("EMK-G3 certificate written:", out)
    print("sha256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
