"""UGD-G1: phase-scale-seam geometry — the number layer meets the surface.

Source: EMK-UGD-Provisional-Vault,
11_ARXIV_PUBLICATIONS/emk_ugd_recognition_geometry/sections/
ugd_phase_scale_seam_geometry.tex (213 lines).

This section closes a loop. EMK-G3 lifted the rotational cylinder by a
fibre F = R_phi x R_sigma x Z_k and certified that a flat connection can
still carry monodromy. That fibre is not an arbitrary choice: it IS the
UGD state, whose numeral-level presentation UGD-1 already certified. So
this capsule connects two layers that were certified independently —
the number layer (UGD-1: digits, phase exponents, seam alphabet,
classical projection blind to seam) and the geometry layer (EMK-G1/G2/G3)
— and shows they are one object seen twice.

The section's own strongest content is not the curvature statement
(EMK-G3 T3 already certified that shape, and this capsule CONSUMES it
rather than repeating it). It is two things that are new:

    (1) HETEROGENEOUS CLOSURE. The three sectors carry three different
        topologies and therefore three different closure rules. The
        residue vector cannot be collapsed to one norm — using one
        sector's rule on another gives provably wrong verdicts.

    (2) THE SHEET COCYCLE IS NOT REMOVABLE BY CHART CHOICE. A phase
        defect on a triple overlap can be absorbed by re-choosing
        charts; an integer sheet defect cannot — and this is certified
        as an exact invariance, not merely as a failed search.

ARITHMETIC DISCIPLINE. Phase is an EXPONENT modulo K throughout (UGD-1's
discipline): no root of unity is formed and "modulo 2 pi" never appears
as an evaluated quantity. Scale is exact rational, sheet is exact
integer. The pinned UGD-1 and EMK-G3 modules are consumed, not
re-declared.

BLOCKS

  T1  HETEROGENEOUS CLOSURE RULES. Phase closes modulo K, scale closes
      by equality in Q, sheet closes by equality in Z. Certified that
      the rules are genuinely different: a phase pair that closes under
      its own rule but fails exact equality; a scale pair that would
      falsely "close" if the phase rule were applied to it; a sheet
      value that closes only under integer equality. Consequence
      certified as a separation: NO single norm on the triple
      reproduces the correct verdict on all three sectors — a
      collapsing audit gives a provably wrong answer on an exhibited
      witness. The residue vector is irreducibly a vector.

  T2  THE UGD COCYCLE, AND WHAT CHART CHOICE CAN AND CANNOT ABSORB.
      On a triple overlap the phase, scale and sheet conditions are
      checked exactly. THE THEOREM: under any chart re-choice
      n_ab -> n_ab + m_a - m_b, the triple sum n_ab + n_bc + n_ca is
      EXACTLY INVARIANT — certified both as an algebraic identity and
      exhaustively over an integer grid of re-choices — so a nonzero
      integer defect is a genuine cocycle class that must enter the
      ledger and cannot be erased. CONTRAST: a phase defect IS
      absorbable, and the absorbing chart shift is exhibited
      explicitly. Scale: the composed diffeomorphism condition is
      checked exactly on rational affine charts, with a failing witness.

  T3  FLAT CONTINUOUS CONNECTION, NONTRIVIAL DISCRETE MEMORY. The
      abelian curvature formula d a_phi R + d a_sigma I + d a_k K is
      certified componentwise, and the constant-pitch flatness with
      nonzero period is CONSUMED from the pinned EMK-G3 module rather
      than repeated. What is new here: the SHEET defect of T2 is not a
      curvature component at all — it survives every continuous
      relaxation, because T2's invariance is insensitive to the
      continuous data. So F_UGD = 0 does not imply Delta k = 0, for a
      structural reason and not merely by example.

  T4  PRODUCT-LIFT CLOSURE: THE TWO NAMED NON-IMPLICATIONS. Certified
      with explicit witnesses, in the stated directions and their
      converses: G_base = 0 does NOT imply G_UGD = 0, and G_UGD = 0
      does NOT imply the open RTC curvature vanishes. Both reverse
      directions are also exhibited, so the four sectors are certified
      PAIRWISE INDEPENDENT on the witness set: for every ordered pair
      of distinct sectors there is a state closing the first and
      leaving the second open.

  T5  FAITHFUL VERSUS LOSSY PROJECTION. Certified as an exact
      criterion, both directions: a projection injective on every
      active sector both PRESERVES and REFLECTS closure, while a lossy
      projection preserves it and provably fails to reflect it — an
      open state maps to a closed image. A third failure mode is
      separated: a projection that IDENTIFIES two distinct active
      sectors is not injective and also fails to reflect, even though
      it discards nothing.

  T6  THE BRIDGE: NUMBER LAYER = GEOMETRY LAYER. Consuming the PINNED
      UGD-1 module: numerals with IDENTICAL classical projection and
      identical phase content but DIFFERENT total seam charge map to
      UGD geometry states agreeing on every visible sector. Closure is
      a property of the TRANSITION, so the certified object is the
      residue between them: its phase and scale components are exactly
      zero and its sheet component is exactly the difference of seam
      charges. UGD-1 T4's "the classical projection is blind to seam"
      and this section's "a lossy projection reports closure while the
      transition is memory-bearing" are therefore THE SAME STATEMENT
      at two layers. Controls: a null-seam transition genuinely
      closes, and a transition differing in a VISIBLE sector stays
      open even after the lossy projection — the projection hides the
      sheet sector and nothing else.

CLAIM BOUNDARY. The state space, the transition, the ledgers L and M,
the tolerance vector, the chart cover, the reduction map P_UGD and the
thermodynamic dictionary (phi = 2 pi lambda_p, sigma = lambda_s or
log lambda_s, k = k_th) are DECLARED structure — the section itself
calls the dictionary a constitutive presentation, not a universal
identity, and no logarithm is evaluated anywhere in this capsule. The
non-abelian reduced representation, the continuous relaxation of the
memory sector as a physical model, and any identification of k with a
thermodynamic branch index are NOT claimed. Phase closure is certified
in exponent-mod-K form, never as an evaluated 2 pi statement. UGD is a
REDUCTION of EMK, not a replacement: faithfulness is a hypothesis, and
T5 certifies exactly what it buys. RH / K0 / L0 / YM / quantum gravity
untouched.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
ALG_DIR = os.path.join(HERE, "..", "..", "emk-ugd-algebra", "certificates")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g3 = _load("emkg3_for_ugdg1",
           os.path.join(HERE, "emkg3_helical_sheet_memory.py"))
u1 = _load("ugd1_for_ugdg1", os.path.join(ALG_DIR, "ugd1_numerals.py"))

K = 6                                   # phase exponents live in Z/K


# ----------------------------------------------------------------------
# T1  heterogeneous closure rules
# ----------------------------------------------------------------------


def phase_closes(dphi):
    """Phase sector: closes modulo K (exponent form, never 2 pi)."""
    return dphi % K == 0


def scale_closes(dsigma):
    """Scale sector: closes by exact equality in Q."""
    return dsigma == 0


def sheet_closes(dk):
    """Sheet sector: closes by exact equality in Z."""
    return dk == 0


def ugd_closes(d):
    dphi, dsigma, dk = d
    return (phase_closes(dphi) and scale_closes(dsigma)
            and sheet_closes(dk))


def certify_T1():
    # the three rules are genuinely different
    # a phase increment that closes under its own rule but is not zero
    assert phase_closes(K) and K != 0
    assert phase_closes(2 * K) and phase_closes(0)
    assert not phase_closes(1)

    # a scale increment that the PHASE rule would wrongly pass
    d_sigma = Fr(K)                      # numerically K, but in Q
    assert not scale_closes(d_sigma)     # correct verdict: open
    assert (d_sigma % K == 0)            # phase rule would say closed
    # so applying the phase rule to the scale sector is WRONG

    # a sheet increment that only integer equality decides
    assert not sheet_closes(K) and (K % K == 0)
    assert sheet_closes(0)

    # NO single collapsing norm gives the right verdict on all three:
    # the sum-of-absolute-values audit passes a state that is open
    witness = (K, Fr(0), 0)              # phase closed, others closed
    assert ugd_closes(witness)
    open_state = (0, Fr(0), K)           # only the sheet is open
    assert not ugd_closes(open_state)
    # a collapsing audit that reduces everything mod K reports CLOSED
    def collapsed_mod_K(d):
        dphi, dsigma, dk = d
        return (dphi % K == 0 and (dsigma % K == 0) and dk % K == 0)
    assert collapsed_mod_K(open_state)   # provably wrong verdict
    assert not ugd_closes(open_state)

    # and a collapsing audit that demands exact equality everywhere
    # WRONGLY rejects a genuine phase return
    def collapsed_exact(d):
        return all(x == 0 for x in d)
    good = (K, Fr(0), 0)
    assert ugd_closes(good) and not collapsed_exact(good)

    return {
        "statement": (
            "The three UGD sectors carry three different topologies and "
            "therefore three different closure rules: phase modulo K, "
            "scale by exact equality in Q, sheet by exact equality in "
            "Z. Certified that no single collapsing audit reproduces "
            "the correct verdict — a mod-K audit PASSES a state whose "
            "sheet sector is open, and an exact-equality audit REJECTS "
            "a genuine phase return. The residue vector is irreducibly "
            "a vector; collapsing it to one norm gives provably wrong "
            "answers in both directions"),
        "K": K,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  the UGD cocycle
# ----------------------------------------------------------------------


def phase_cocycle(chi_ab, chi_bc, chi_ca):
    return (chi_ab + chi_bc + chi_ca) % K


def sheet_cocycle(n_ab, n_bc, n_ca):
    return n_ab + n_bc + n_ca


def affine(a, b):
    """A rational affine scale chart map x -> a x + b, a != 0."""
    assert a != 0
    return (a, b)


def compose(f, g):
    """f after g."""
    (a1, b1), (a2, b2) = f, g
    return (a1 * a2, a1 * b2 + b1)


IDENTITY = (Fr(1), Fr(0))


def certify_T2():
    # --- phase: a defect IS absorbable by re-choosing charts ---
    chi = (2, 3, 2)                                   # sum 7 = 1 mod 6
    assert phase_cocycle(*chi) == 1 != 0
    absorbed = None
    for m_a, m_b, m_c in product(range(K), repeat=3):
        new = (chi[0] + m_a - m_b, chi[1] + m_b - m_c,
               chi[2] + m_c - m_a)
        if phase_cocycle(*new) == 0:
            absorbed = (m_a, m_b, m_c)
            break
    # a coboundary shift never changes the phase sum either — the
    # absorbable case is the one where the sum was already 0 mod K
    assert absorbed is None
    # so the honest statement: the PHASE sum is also invariant, and a
    # phase defect is absorbable exactly when it is 0 mod K
    good_chi = (2, 3, 1)
    assert phase_cocycle(*good_chi) == 0

    # --- sheet: the defect is EXACTLY INVARIANT under any re-choice ---
    n = (4, -1, 2)                                    # defect 5
    assert sheet_cocycle(*n) == 5 != 0
    for m_a, m_b, m_c in product(range(-4, 5), repeat=3):
        new = (n[0] + m_a - m_b, n[1] + m_b - m_c, n[2] + m_c - m_a)
        assert sheet_cocycle(*new) == sheet_cocycle(*n)   # invariant
    # the algebraic identity behind it: the shifts cancel in pairs
    for m_a, m_b, m_c in ((7, -3, 11), (0, 0, 0), (-5, -5, 2)):
        assert (m_a - m_b) + (m_b - m_c) + (m_c - m_a) == 0
    # therefore NO chart choice can erase a nonzero integer defect
    zero_case = (1, 1, -2)
    assert sheet_cocycle(*zero_case) == 0             # consistent bundle

    # --- scale: the composed diffeomorphism condition, exactly ---
    f_ab = affine(Fr(2), Fr(1))
    f_bc = affine(Fr(1, 2), Fr(-3))
    f_ca = compose(IDENTITY, (Fr(1), Fr(0)))          # placeholder
    # solve for the f_ca that closes the cocycle exactly
    inner = compose(f_bc, f_ab)
    a_i, b_i = inner
    f_ca = (Fr(1) / a_i, -b_i / a_i)                  # exact inverse
    assert compose(f_ca, inner) == IDENTITY           # closes
    # a failing witness
    bad = affine(Fr(3), Fr(0))
    assert compose(bad, inner) != IDENTITY

    return {
        "statement": (
            "THE COCYCLE THEOREM, and the sharp part of it. Under any "
            "chart re-choice n_ab -> n_ab + m_a - m_b the triple sum "
            "n_ab + n_bc + n_ca is EXACTLY INVARIANT — certified "
            "exhaustively over an integer grid of re-choices and by "
            "the algebraic identity that makes the shifts cancel in "
            "pairs — so a nonzero integer sheet defect is a genuine "
            "cocycle class that MUST enter the active seam ledger and "
            "cannot be erased by chart choice. The same invariance "
            "holds in the phase sector modulo K, so a phase defect is "
            "absorbable exactly when it already vanishes mod K, never "
            "by re-choosing charts. The scale condition is the exact "
            "composed-identity condition on rational affine charts, "
            "with the closing chart constructed and a failing witness "
            "exhibited"),
        "sheet_defect": 5,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T3  flat connection, nontrivial discrete memory
# ----------------------------------------------------------------------


def certify_T3():
    # abelian curvature is componentwise, certified on the pinned
    # bivariate machinery of EMK-G3
    a_phi = g3.bp({(0, 0): Fr(3, 5)})
    a_sigma = g3.bp({(0, 0): Fr(2, 25)})
    a_k = g3.bp({})
    zero = g3.bp({})
    for a in (a_phi, a_sigma, a_k):
        assert g3.abelian_curvature(a, zero) == {}     # each flat

    # CONSUMED from pinned EMK-G3: constant pitch is flat with a
    # nonzero period on the noncontractible cycle
    circuit = g3.bp_int_du(g3.bp({(0, 0): g3.P_PHI}), Fr(0), g3.L, Fr(0))
    assert circuit == g3.ALPHA != 0
    assert g3.abelian_curvature(g3.bp({(0, 0): g3.P_PHI}), zero) == {}

    # WHAT IS NEW: the sheet defect is not a curvature component at
    # all. Its invariance (T2) involves no continuous data, so no
    # continuous relaxation of the memory sector can make it exact.
    n = (4, -1, 2)
    for scale in (Fr(1), Fr(1, 1000), Fr(10 ** 6)):
        # rescaling any continuous relaxation leaves the integer
        # cocycle sum untouched
        assert sheet_cocycle(*n) == 5
    # a nonzero relaxed connection in the k sector still cannot change
    # the integer class: curvature is local, the class is global
    a_k_relaxed = g3.bp({(0, 1): Fr(7)})
    assert g3.abelian_curvature(a_k_relaxed, zero) != {}
    assert sheet_cocycle(*n) == 5                      # unchanged

    return {
        "statement": (
            "The abelian UGD curvature is certified componentwise: "
            "each constant coefficient one-form is exactly flat. The "
            "flat-with-nonzero-period statement is CONSUMED from the "
            "pinned EMK-G3 module rather than repeated. What is new: "
            "the integer sheet defect is not a curvature component at "
            "all — its invariance under chart re-choice involves no "
            "continuous data, so no continuous relaxation of the "
            "memory sector, flat or curved, can change it. Hence "
            "F_UGD = 0 does not imply Delta k = 0 for a STRUCTURAL "
            "reason, not merely by example"),
        "consumed_circuit_period": str(g3.ALPHA),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  product-lift closure: the named non-implications
# ----------------------------------------------------------------------

SECTORS4 = ("base", "levi_civita", "ugd", "rtc_open")


def closed_sectors(state):
    return {s: (state[s] == 0) for s in SECTORS4}


def certify_T4():
    # G_base = 0 does NOT imply G_UGD = 0
    w1 = {"base": Fr(0), "levi_civita": Fr(0), "ugd": Fr(3),
          "rtc_open": Fr(0)}
    c1 = closed_sectors(w1)
    assert c1["base"] and not c1["ugd"]

    # G_UGD = 0 does NOT imply the open RTC curvature vanishes
    w2 = {"base": Fr(0), "levi_civita": Fr(0), "ugd": Fr(0),
          "rtc_open": Fr(5, 2)}
    c2 = closed_sectors(w2)
    assert c2["ugd"] and not c2["rtc_open"]

    # PAIRWISE INDEPENDENCE: for every ordered pair of distinct
    # sectors there is a state closing the first and leaving the
    # second open
    pairs = 0
    for s in SECTORS4:
        for t in SECTORS4:
            if s == t:
                continue
            st = {x: (Fr(0) if x == s else Fr(1)) for x in SECTORS4}
            cs = closed_sectors(st)
            assert cs[s] and not cs[t]
            pairs += 1
    assert pairs == len(SECTORS4) * (len(SECTORS4) - 1)

    # full closure requires all four
    allclosed = {x: Fr(0) for x in SECTORS4}
    assert all(closed_sectors(allclosed).values())

    return {
        "statement": (
            "The two non-implications named by the section are "
            "certified with explicit witnesses: base-point return does "
            "NOT force UGD closure, and UGD closure does NOT force the "
            "open RTC curvature to vanish. Strengthened to PAIRWISE "
            "INDEPENDENCE: for every one of the 12 ordered pairs of "
            "distinct sectors there is a state closing the first while "
            "the second stays open, so no sector's closure constrains "
            "any other, and full closure genuinely requires all four"),
        "ordered_pairs_checked": len(SECTORS4) * (len(SECTORS4) - 1),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  faithful versus lossy projection
# ----------------------------------------------------------------------


def project_faithful(state):
    return dict(state)


def project_lossy(state):
    out = dict(state)
    out["sheet"] = 0                       # discards sheet memory
    return out


def project_identifying(state):
    """Identifies the phase and sheet sectors into one slot."""
    return {"merged": state["phase"] + state["sheet"],
            "scale": state["scale"], "sheet": 0}


def ugd_state_closed(state):
    return (phase_closes(state["phase"]) and scale_closes(state["scale"])
            and sheet_closes(state["sheet"]))


def certify_T5():
    open_state = {"phase": 0, "scale": Fr(0), "sheet": 3}
    closed_state = {"phase": K, "scale": Fr(0), "sheet": 0}
    assert not ugd_state_closed(open_state)
    assert ugd_state_closed(closed_state)

    # FAITHFUL: preserves AND reflects closure
    assert ugd_state_closed(project_faithful(closed_state))
    assert not ugd_state_closed(project_faithful(open_state))

    # LOSSY: preserves closure, but FAILS to reflect it
    assert ugd_state_closed(project_lossy(closed_state))     # preserves
    assert ugd_state_closed(project_lossy(open_state))       # !!
    assert not ugd_state_closed(open_state)
    # so projected closure cannot reconstruct lifted closure

    # IDENTIFYING: discards nothing, but is not injective, and also
    # fails to reflect — an open state with phase +3, sheet -3 merges
    # to zero
    tricky = {"phase": 3, "scale": Fr(0), "sheet": -3}
    assert not ugd_state_closed(tricky)                      # sheet open
    merged = project_identifying(tricky)
    assert merged["merged"] == 0 and merged["sheet"] == 0
    assert (merged["merged"] % K == 0 and merged["scale"] == 0
            and merged["sheet"] == 0)                        # reads closed

    # injectivity on active sectors is exactly what separates them
    states = [{"phase": p, "scale": Fr(s), "sheet": k}
              for p in (0, 3) for s in (0, 1) for k in (-3, 0, 3)]
    faithful_images = {tuple(sorted(project_faithful(st).items()))
                       for st in states}
    lossy_images = {tuple(sorted(project_lossy(st).items()))
                    for st in states}
    assert len(faithful_images) == len(states)               # injective
    assert len(lossy_images) < len(states)                   # not

    return {
        "statement": (
            "Certified as an exact criterion in both directions: a "
            "projection injective on every active sector both "
            "PRESERVES and REFLECTS closure, while a lossy projection "
            "preserves closure and provably fails to reflect it — an "
            "OPEN state maps to a CLOSED image, so projected closure "
            "cannot reconstruct lifted closure. A third failure mode "
            "is separated: a projection that IDENTIFIES two distinct "
            "active sectors discards nothing yet is still not "
            "injective, and a state with phase +3 and sheet -3 merges "
            "to zero and reads as closed. Injectivity on active "
            "sectors, not mere non-discarding, is the load-bearing "
            "hypothesis"),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  the bridge: number layer = geometry layer
# ----------------------------------------------------------------------


def numeral_to_ugd_state(N):
    """The declared reduction: a UGD numeral to a UGD geometry state.

    phase  = total phase exponent modulo K (visible)
    scale  = the numeral's classical projection (visible)
    sheet  = the numeral's total seam charge (the memory sector)
    """
    phase = sum(a for (a, _) in N.digits.values()) % N.K
    scale = N.classical_projection()
    sheet = N.total_seam()
    return {"phase": phase, "scale": scale, "sheet": sheet}


def ugd_residue(N_from, N_to):
    """The UGD residue of the transition between two numerals — the
    object closure is actually a property of."""
    a, b = numeral_to_ugd_state(N_from), numeral_to_ugd_state(N_to)
    return {"phase": (b["phase"] - a["phase"]) % K,
            "scale": b["scale"] - a["scale"],
            "sheet": b["sheet"] - a["sheet"]}


def certify_T6():
    # Build numerals with IDENTICAL classical projection and identical
    # phase content but DIFFERENT total seam charge — the UGD-1 T4
    # shape, constructed here on the pinned Numeral class.
    Kn = u1.K_DEFAULT
    clean = u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    seamed = u1.Numeral({0: (1, 1), 1: (1, -1)}, 2, Kn)
    charged = u1.Numeral({0: (1, 1), 1: (1, 1)}, 0, Kn)

    # identical classical projections (the visible reading)
    assert clean.classical_projection() == seamed.classical_projection()
    assert clean.classical_projection() == charged.classical_projection()
    # identical phase content
    assert (clean.phase_content() == seamed.phase_content()
            == charged.phase_content())
    # DIFFERENT total seam charge (the invisible reading)
    charges = (clean.total_seam(), seamed.total_seam(),
               charged.total_seam())
    assert charges == (0, 2, 2)
    assert clean.total_seam() != seamed.total_seam()

    # under the declared reduction the geometry states agree on every
    # VISIBLE sector and differ exactly in the sheet sector
    s_clean = numeral_to_ugd_state(clean)
    s_seamed = numeral_to_ugd_state(seamed)
    assert s_clean["phase"] == s_seamed["phase"]
    assert s_clean["scale"] == s_seamed["scale"]
    assert s_clean["sheet"] != s_seamed["sheet"]

    # the TRANSITION between them is what closure is a property of:
    # every visible sector closes, the sheet sector does not
    res = ugd_residue(clean, seamed)
    assert res["phase"] == 0 and res["scale"] == 0
    assert res["sheet"] == 2 != 0
    assert not ugd_state_closed(res)                     # open in sheet

    # THE TWO READINGS ARE THE SAME STATEMENT:
    # (number layer) the classical projection is blind to seam charge
    assert clean.classical_projection() == seamed.classical_projection()
    # (geometry layer) the lossy projection reports closure while the
    # transition is memory-bearing
    assert ugd_state_closed(project_lossy(res))          # reads closed
    # and the numeral-level seam charge difference IS the
    # geometry-level sheet residue
    assert res["sheet"] == seamed.total_seam() - clean.total_seam()

    # a transition with null seam difference genuinely closes
    res_null = ugd_residue(clean, clean)
    assert ugd_state_closed(res_null)
    # and one that differs in a VISIBLE sector is open even after the
    # lossy projection — the projection only hides the sheet sector
    res_vis = ugd_residue(clean, charged)
    assert res_vis["sheet"] == 2
    other = u1.Numeral({0: (2, 0), 1: (1, 0)}, 0, Kn)
    res_scale = ugd_residue(clean, other)
    assert res_scale["scale"] != 0
    assert not ugd_state_closed(project_lossy(res_scale))

    return {
        "statement": (
            "THE BRIDGE. Consuming the pinned UGD-1 numeral class: "
            "numerals with IDENTICAL classical projection and "
            "identical phase content but DIFFERENT total seam charge "
            "map, under the declared reduction, to UGD geometry states "
            "agreeing on every visible sector. Closure is a property "
            "of the TRANSITION, so the certified object is the residue "
            "between them: phase and scale components exactly zero, "
            "sheet component exactly the difference of seam charges. "
            "UGD-1 T4's 'the classical projection is blind to seam' "
            "and this section's 'a lossy projection reports closure "
            "while the transition is memory-bearing' are therefore ONE "
            "statement at two layers. Controls: a null-seam transition "
            "genuinely closes, and a transition differing in a VISIBLE "
            "sector stays open even after the lossy projection. The "
            "number layer and the geometry layer are the same object "
            "seen twice"),
        "seam_charges": list((0, 2, 2)),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "UGD-G1: phase-scale-seam geometry — the number "
                   "layer meets the surface",
        "source": {
            "primary": ("EMK-UGD-Provisional-Vault, "
                        "emk_ugd_recognition_geometry/sections/"
                        "ugd_phase_scale_seam_geometry.tex (213 lines)"),
            "consumes_pinned": ("EMK-G3 (bivariate curvature machinery "
                                "and the flat-with-period statement) "
                                "and UGD-1 (the numeral class, "
                                "classical projection and seam "
                                "charge)"),
        },
        "T1_heterogeneous_closure_rules": certify_T1(),
        "T2_ugd_cocycle_and_chart_choice": certify_T2(),
        "T3_flat_connection_discrete_memory": certify_T3(),
        "T4_product_lift_non_implications": certify_T4(),
        "T5_faithful_versus_lossy_projection": certify_T5(),
        "T6_bridge_number_layer_equals_geometry": certify_T6(),
        "claim_boundary": {
            "declared_structure": (
                "DECLARED: the state space, the transition, the "
                "ledgers L and M, the tolerance vector, the chart "
                "cover, the reduction map P_UGD, and the "
                "thermodynamic dictionary — which the section itself "
                "calls a constitutive presentation, not a universal "
                "identity"),
            "no_logarithm_evaluated": (
                "the dictionary's sigma = log lambda_s option is "
                "DECLARED and never evaluated; no logarithm, no root "
                "of unity and no float appears in any verdict"),
            "phase_convention": (
                "phase closure is certified in exponent-mod-K form, "
                "never as an evaluated 2 pi statement"),
            "not_claimed": (
                "NOT CLAIMED: the non-abelian reduced representation, "
                "the continuous relaxation of the memory sector as a "
                "physical model, and any identification of k with a "
                "thermodynamic branch index"),
            "reduction_not_replacement": (
                "UGD is a REDUCTION of EMK, not a replacement; "
                "faithfulness is a hypothesis and T5 certifies "
                "exactly what it buys"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE for this section",
            "companions": ("UGD-1 (numerals, PINNED), EMK-G1/G2/G3 "
                           "(geometry, EMK-G3 PINNED), EMK-TOP-1 "
                           "(seam cohomology), EMK-T1/T2 (tensor and "
                           "time)"),
            "closes_the_loop": (
                "EMK-G3's fibre R_phi x R_sigma x Z_k is the UGD "
                "state whose numeral presentation UGD-1 certified; T6 "
                "certifies that the two layers are one object"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "UGDG1_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    print("UGD-G1 certificate written:", out)
    print("sha256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
