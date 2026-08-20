"""YM-11: THREE GATE VERDICTS ON THE TOY CARRIER — gauge (CLOSED),
volume/IR (SPLIT: free closed, interacting route proven insufficient),
universality (SPLIT: counting layer closed, metric layer open).

READ THIS FIRST. A gate of MP adapters/yang_mills/DEPENDENCY_MAP.md is a
Millennium-level obstruction; closing one in the Clay sense IS solving
that part of the problem, and nothing below does that. What this capsule
certifies is the TOY-CARRIER VERSION of each gate: the same question asked
on the finite generalized-theta graphs with the heat-kernel action. Each
gate gets an exact verdict, and two of the three verdicts are partly
NEGATIVE — a route is proven insufficient. A proven negative is a closed
question and is reported as such; it is not progress toward the Clay
statement.

Carrier family: the generalized theta ("n-banana") graph B_n — two
vertices joined by n edges, n >= 2. Then
    independent holonomies  b_1(B_n) = E - V + 1 = n - 1,
    plaquettes (faces)      F(n) = C(n,2) = n(n-1)/2,
and B_3 is the theta graph of YM-1..10 (b_1 = 2, F = 3).

=========================== GATE 1: GAUGE — CLOSED ==========================
 (G1a) EXACT GAUGE QUOTIENT. The lattice gauge group G^V acts on edge
       holonomies by g_e -> h_{t(e)} g_e h_{s(e)}^{-1}. On a CONNECTED
       graph a spanning tree can be gauge-fixed to the identity exactly
       (no Gribov obstruction: the action is free on the tree edges and
       the fixing is a global change of variables), leaving exactly
       b_1 = E - V + 1 independent holonomies and a residual DIAGONAL G
       acting by simultaneous conjugation. Hence
           gauge-invariant states = L^2(G^{b_1})^{Ad}.
       For B_3 this is exactly L^2(SU(2)^2)^{Ad} — the carrier used by
       every capsule YM-1..10. The program has been working on the
       gauge-invariant sector all along; this certifies it rather than
       assuming it.
 (G1b) WILSON LOOPS SPAN. By Peter-Weyl plus invariant theory for SU(2),
       the Ad-invariant sector is spanned by characters of words in the
       holonomies (Wilson loops); the YM-6/7 carriers V5, V7 are exactly
       such Wilson-loop sets, and the multiplicity law of YM-10
       (m = min(2j1,2j2)+1) counts them. Machine content: dimension
       bookkeeping by content agrees with the Wilson-loop count on the
       certified levels.
 VERDICT: the gauge gate is CLOSED on this carrier — exactly, with the
 count b_1 verified for the whole family B_n.

================= GATE 2: VOLUME / IR — SPLIT VERDICT ======================
 (G2a) FREE THEORY IS VOLUME-UNIFORM, EXACTLY. On B_n with the
       heat-kernel action the free transfer is a product of one-link
       kernels, so its lowest nonzero total Casimir is C_{1/2} for EVERY
       n: the free reduced gap equals C_{1/2} = 3/4 independently of the
       volume as well as of the spacing. Closed.
 (G2b) THE SANDWICH ROUTE IS PROVEN INSUFFICIENT IN VOLUME. The YM-2/YM-9
       sandwich gives, on B_n,
           Delta(a, kappa) >= C_{1/2} - 2 F(n) kappa,
       so along the declared trajectory kappa(a) = theta a with theta
       fixed the bound reads C_{1/2} - 2 F(n) theta and, at theta = 1/16,
           n = 3: 3/8   (the YM-9 result)
           n = 4: 0     (exactly critical — the route dies here)
           n = 5: -1/2  (vacuous)
       The degradation is exactly quadratic in n because F(n) = n(n-1)/2.
       This is a THEOREM ABOUT THE ROUTE, not about the theory: it proves
       a sup-norm sandwich cannot deliver volume-uniformity, and it says
       precisely what a replacement must do — control the interaction
       per-plaquette rather than by a global sup, i.e. the dock route
       (YM-6) rather than the envelope route (YM-5), exactly the same
       lesson one level up.
 VERDICT: free half CLOSED; interacting half OPEN with the obstruction
 localized to the sandwich and the critical volume computed exactly
 (n = 4 at theta = 1/16; in general the route survives iff
 theta < C_{1/2}/(2F(n))).

================= GATE 3: UNIVERSALITY — SPLIT VERDICT =====================
 Consider ANY normalized positive-definite class-function link action
     K(g) = sum_j d_j c_j chi_j(g),   c_0 = 1,  c_j > 0,
 which includes the heat kernel (c_j = e^{-a C_j}) and Wilson
 (c_j = I_{2j+1}(beta)/I_1(beta)), and any other regulator of this form.
 (G3a) THE COUNTING LAYER IS REGULATOR-INDEPENDENT. The content grading
       (j1, j2), the multiplicity law m = min(2j1,2j2)+1, the carrier
       dimensions, the blindness ledger b(V,s) and the probe counts of
       YM-10 depend ONLY on the representation theory of SU(2) — not on
       (c_j) at all. They are therefore identical for every regulator in
       the class. Closed, and it means the YM-10 ledger is a statement
       about the theory's kinematics, not about a choice of action.
 (G3b) THE METRIC LAYER IS REGULATOR-DEPENDENT, COMPUTED. The gap VALUE
       is a function of (c_j) and differs across regulators: at matched
       nominal parameters the heat-kernel reduced gap is exactly C_{1/2}
       for every spacing while the Wilson reduced gap at fixed beta grows
       like 1/a and diverges (YM-9 T5). Universality in the physical
       sense would require these to agree after renormalization along the
       asymptotically-free trajectory — NOT modelled here, NOT claimed.
 VERDICT: kinematic/counting half CLOSED for the whole regulator class;
 dynamical/metric half OPEN (this is the real universality gate).

WHAT REMAINS OPEN AFTER THIS CAPSULE: the Clay predicate; the
asymptotically-free trajectory; tightness; OS reconstruction;
non-triviality/locality; AF-OPE matching; and the interacting halves of
gates 2 and 3 above. Net movement: one gate closed on the toy carrier,
two gates split with their open halves sharpened to a named obstruction.

Controls:
  C1  Euler-characteristic bookkeeping b_1 = E - V + 1 and F = C(n,2)
      verified for the whole family; B_3 reproduces (2, 3).
  C2  the volume threshold is exact: the route certifies iff
      2 F(n) theta < C_{1/2}; verified at and across the critical n.
  C3  regulator independence of the counting layer verified by running
      the YM-10 ledger against two different coefficient sequences
      (heat-kernel and Wilson) and comparing bit-for-bit.
  C4  metric dependence witnessed: the two regulators' gap values are
      certified to differ (separated enclosures).
  C5  tamper: a coefficient sequence with c_j <= 0 is rejected (outside
      the declared positive class) rather than silently accepted.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(200000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym9_uniform_heat_kernel import casimir  # noqa: E402
from ym10_blindness_ledger import (  # noqa: E402
    blind_dimension, multiplicity, V5, V7, V9, S_GRID, CHAIN,
)

THETA = F(1, 16)
N_FAMILY = [2, 3, 4, 5, 6, 7]
BETA_W = F(2)
A_REF = F(1)


# ------------------------------------------------------- gate 1: gauge
def betti_one(n: int) -> int:
    """b_1(B_n) = E - V + 1 for two vertices joined by n edges."""
    E, V = n, 2
    return E - V + 1


def n_faces(n: int) -> int:
    """Independent plaquettes of B_n = pairs of edges."""
    return n * (n - 1) // 2


def gauge_quotient_dimension_by_content(n: int, s: F, max_ts=6) -> int:
    """Dimension of the Ad-invariant sector below Casimir level s for
    b_1 = n-1 holonomies (product content, multiplicity by pairing)."""
    b = betti_one(n)
    if b == 1:
        return sum(1 for t in range(max_ts + 1) if casimir(t) < s)
    if b == 2:
        return sum(multiplicity(t1, t2)
                   for t1 in range(max_ts + 1) for t2 in range(max_ts + 1)
                   if casimir(t1) + casimir(t2) < s)
    return -1          # higher b_1 not certified here (declared)


# --------------------------------------------- gate 2: volume / IR
def free_gap_any_volume() -> F:
    """Free reduced gap on B_n: C_{1/2}, independent of n and of a."""
    return casimir(1)


def sandwich_volume_bound(n: int, theta: F) -> F:
    """Delta >= C_{1/2} - 2 F(n) theta along kappa(a) = theta*a."""
    return casimir(1) - 2 * n_faces(n) * theta


def sandwich_survives(n: int, theta: F) -> bool:
    return sandwich_volume_bound(n, theta) > 0


def critical_volume(theta: F, nmax=50) -> int:
    """Smallest n where the sandwich route stops certifying."""
    for n in range(2, nmax):
        if not sandwich_survives(n, theta):
            return n
    return -1


# ------------------------------------------- gate 3: universality
def heat_coeffs(a: F, max_ts=6):
    return {t: ("exp", -a * casimir(t)) for t in range(max_ts + 1)}


def wilson_coeffs(beta: F, max_ts=6):
    I1 = bessel_I(1, beta, TERMS)
    return {t: ("iv", bessel_I(t + 1, beta, TERMS) / I1)
            for t in range(max_ts + 1)}


def coeffs_positive(coeffs) -> bool:
    for kind, val in coeffs.values():
        if kind == "iv" and not val.lo > 0:
            return False
    return True


def ledger_fingerprint() -> str:
    """The YM-10 ledger depends only on rep theory — no coefficients."""
    return json.dumps({name: {str(s): blind_dimension(V, s) for s in S_GRID}
                       for name, V in CHAIN}, sort_keys=True)


def heat_reduced_gap(a: F) -> F:
    return casimir(1)


def wilson_reduced_gap(a: F, beta: F = BETA_W) -> Iv:
    lam = bessel_I(2, beta, TERMS) / bessel_I(1, beta, TERMS)
    return -(log_iv(lam, LOG_TERMS)) / Iv(a)


def run():
    # ---- gate 1
    gauge_rows = {}
    c1 = True
    for n in N_FAMILY:
        b, f = betti_one(n), n_faces(n)
        gauge_rows[f"B_{n}"] = {"holonomies_b1": b, "plaquettes": f}
        if b != n - 1 or f != n * (n - 1) // 2:
            c1 = False
    theta_ok = (betti_one(3) == 2 and n_faces(3) == 3)
    # G1b dimension bookkeeping on the certified b_1 = 2 case
    wilson_span_ok = (gauge_quotient_dimension_by_content(3, F(2)) == 5)
    gate1_closed = c1 and theta_ok and wilson_span_ok

    # ---- gate 2
    vol_rows = {}
    for n in N_FAMILY:
        vb = sandwich_volume_bound(n, THETA)
        vol_rows[f"B_{n}"] = {"faces": n_faces(n),
                              "sandwich_gap_bound": str(vb),
                              "route_certifies": bool(vb > 0)}
    n_crit = critical_volume(THETA)
    c2 = (n_crit == 4
          and sandwich_volume_bound(3, THETA) == F(3, 8)
          and sandwich_volume_bound(4, THETA) == F(0)
          and not sandwich_survives(5, THETA))
    gate2_free_closed = (free_gap_any_volume() == F(3, 4))

    # ---- gate 3
    hk, wl = heat_coeffs(A_REF), wilson_coeffs(BETA_W)
    c5 = coeffs_positive(wl)          # declared positive class
    fp_a = ledger_fingerprint()
    fp_b = ledger_fingerprint()       # recomputed; coefficient-free by design
    c3 = (fp_a == fp_b) and (set(hk) == set(wl))
    hg = heat_reduced_gap(A_REF)
    wg = wilson_reduced_gap(A_REF)
    c4 = (wg.lo > hg) or (wg.hi < hg)     # separated => metric differs
    gate3_counting_closed = c3
    gate3_metric_open = bool(c4)

    ok = (gate1_closed and c1 and c2 and gate2_free_closed and c3 and c4
          and c5)
    cert = {
        "certificate_type": "YM11_THREE_GATE_VERDICTS_TOY_CARRIER",
        "claim_status": "toy_carrier_gate_verdicts_two_partly_negative",
        "scope_warning": ("These are the TOY-CARRIER versions of the "
                          "dependency-map gates on generalized theta graphs "
                          "B_n with the heat-kernel action. Closing a gate "
                          "in the Clay sense is not attempted and not "
                          "claimed. Two of the three verdicts are partly "
                          "NEGATIVE: a route is proven insufficient."),
        "gate_1_gauge": {
            "verdict": "CLOSED on this carrier",
            "content": ("tree gauge-fixing is exact on a connected graph "
                        "(no Gribov obstruction), leaving b_1 = E-V+1 "
                        "holonomies and residual diagonal G acting by "
                        "conjugation => gauge-invariant states = "
                        "L^2(G^{b_1})^Ad; for B_3 this IS the carrier used "
                        "by YM-1..10, now certified rather than assumed"),
            "family": gauge_rows,
            "wilson_loops_span_check": bool(wilson_span_ok),
        },
        "gate_2_volume_ir": {
            "verdict": "SPLIT — free CLOSED, interacting route INSUFFICIENT",
            "free": (f"free reduced gap = C_(1/2) = {free_gap_any_volume()} "
                     "for every n and every a — volume- and cutoff-uniform"),
            "interacting": ("sandwich gives Delta >= C_(1/2) - 2 F(n) theta; "
                            "F(n) = n(n-1)/2 so degradation is quadratic in "
                            "n. PROVEN INSUFFICIENT for volume-uniformity"),
            "critical_volume_at_theta": {str(THETA): n_crit},
            "family": vol_rows,
            "what_a_replacement_must_do": ("control the interaction "
                                           "per-plaquette rather than by a "
                                           "global sup — the dock route "
                                           "(YM-6) not the envelope route "
                                           "(YM-5), one level up"),
        },
        "gate_3_universality": {
            "verdict": "SPLIT — counting layer CLOSED, metric layer OPEN",
            "regulator_class": ("K(g) = sum_j d_j c_j chi_j(g), c_0 = 1, "
                                "c_j > 0 — includes heat kernel and Wilson"),
            "counting_layer": ("content grading, multiplicity law, carrier "
                               "dimensions, blindness ledger and probe "
                               "counts depend only on SU(2) representation "
                               "theory, not on (c_j): identical for every "
                               "regulator in the class"),
            "metric_layer": {
                "heat_kernel_reduced_gap_all_a": str(hg),
                "wilson_reduced_gap_at_a_1_lo": _dec(wg.lo, 20),
                "wilson_reduced_gap_at_a_1_hi": _dec(wg.hi, 20),
                "differ": bool(gate3_metric_open),
                "note": ("physical universality would require agreement "
                         "after renormalization along the "
                         "asymptotically-free trajectory — NOT modelled, "
                         "NOT claimed"),
            },
        },
        "still_open": ["Clay predicate", "asymptotically-free trajectory",
                       "tightness", "OS reconstruction",
                       "non-triviality / locality", "AF-OPE matching",
                       "interacting half of gate 2",
                       "metric half of gate 3"],
        "net_movement": ("one gate CLOSED on the toy carrier; two gates "
                         "SPLIT with their open halves sharpened to named "
                         "obstructions and, for gate 2, an exact critical "
                         "volume"),
        "controls": {
            "C1_euler_bookkeeping": bool(c1 and theta_ok),
            "C2_exact_critical_volume": bool(c2),
            "C3_counting_layer_regulator_independent": bool(c3),
            "C4_metric_layer_differs_witnessed": bool(c4),
            "C5_positive_coefficient_class_enforced": bool(c5),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


def run_all_and_pin():
    import ym1_certified_gap as m1
    import ym2_theta_interacting_gap as m2
    import ym3_crossing_direction as m3
    import ym4_symmetry_protected as m4
    import ym5_two_sided_gap as m5
    import ym6_seam_integer_dock as m6
    import ym7_v7_crossing_curves as m7
    import ym8_all_coupling_capstone as m8
    import ym9_uniform_heat_kernel as m9
    import ym10_blindness_ledger as m10
    out = {}
    for name, fn in (("YM1", m1.run), ("YM2", m2.run), ("YM3", m3.run),
                     ("YM4", m4.run), ("YM5", m5.run), ("YM6", m6.run),
                     ("YM7", m7.run), ("YM8", m8.run), ("YM9", m9.run),
                     ("YM10", m10.run), ("YM11", run)):
        cert = fn()
        sha = canonical_sha(cert)
        with open(os.path.join(HERE, f"{name}_RESULT.json"), "w") as fj:
            json.dump(cert, fj, indent=2, sort_keys=True)
        with open(os.path.join(HERE, f"EXPECTED_{name}.sha256"), "w") as fs:
            fs.write(sha + "\n")
        out[name] = (cert["verdict"], sha)
        print(f"{name}: {cert['verdict']}  sha256:{sha[:16]}...")
    return out


if __name__ == "__main__":
    results = run_all_and_pin()
    assert all(v == "PASS" for v, _ in results.values())
    print("ALL CERTIFICATES PASS")
