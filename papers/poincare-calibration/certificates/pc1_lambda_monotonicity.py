"""PC-1: PERELMAN'S LAMBDA-MONOTONICITY, EXACTLY REALIZED ON THE SU(2)
HOMOGENEOUS CARRIER — a calibration capsule for the Poincare adapter.

PURPOSE — CALIBRATION, NOT DISCOVERY. The Poincare conjecture is a proved
theorem (Hamilton-Perelman, 2003). Nothing here claims any part of it, and
nothing here is new mathematics. This capsule exists to test whether the
adapter framework can carry a genuine piece of the Ricci-flow architecture
EXACTLY, without hiding any step inside a decorative closure predicate —
the stated purpose of MP adapters/poincare/adapter.tex v0.2 ("Poincare
Calibration Adapter"), whose dictionary this capsule instantiates.

Carrier (deliberately the same group as the YM line): left-invariant
metrics on SU(2), diagonal in a Milnor frame, g = diag(A, B, C), A,B,C > 0
("Berger spheres" when two entries agree; round S^3 when all three do).
Ricci flow on this carrier is an exact rational ODE system:
    A' = 2[(B-C)^2 - A^2]/(BC),  and cyclically,
i.e. A' = -2 r_1 with the Ricci eigenvalues
    r_1 = [A^2 - (B-C)^2]/(BC),  r_2 = [B^2 - (C-A)^2]/(CA),
    r_3 = [C^2 - (A-B)^2]/(AB).

CERTIFIED (all exact — symbolic rational-function identities plus exact
rational arithmetic; no floats in any verdict):

 (T1) LAMBDA REDUCES TO SCALAR CURVATURE. Perelman's functional
        lambda(g) = inf { Int (R + |grad f|^2) e^{-f} dV : Int e^{-f} dV = 1 }
                  = min spec( -4 Delta_g + R )
      satisfies lambda(g) = R(g) on ANY homogeneous metric: R is constant,
      -4 Delta is a nonnegative operator whose kernel is the constants, so
      the bottom of the spectrum is attained at the constant function and
      equals R. Machine content: R is verified constant on the carrier and
      equal to the trace of the Ricci eigenvalues,
        R = [2(AB+BC+CA) - (A^2+B^2+C^2)]/(ABC).

 (T2) EXACT MONOTONICITY IDENTITY (this is Perelman's lambda-monotonicity
      on this carrier). As a rational-function identity in Q(A,B,C):
        dR/dt  -  2 |Ric|^2  =  0        (verified symbolically, remainder
                                          identically zero)
      with the explicit SUM-OF-SQUARES witness
        |Ric|^2 = (r_1/A)^2 + (r_2/B)^2 + (r_3/C)^2 >= 0,
      each summand a square of a rational function. Hence
        d lambda / dt = 2 |Ric|^2 >= 0
      along Ricci flow, with equality iff Ric = 0 (impossible on SU(2) with
      positive-definite g). Perelman's monotonicity is therefore EXACT and
      strict here, and its proof on this carrier is a polynomial identity —
      no analysis, no approximation, no declared tolerance.

 (T3) EXACT ROUND SOLUTION AND EXTINCTION. For A=B=C=a: a(t) = a_0 - 2t,
      R(t) = 3/(a_0 - 2t), extinction time t* = a_0/2 EXACTLY (rational),
      lambda strictly increasing and blowing up at t*. Verified as exact
      rational arithmetic at rational sample times.

 (T4) EXACT ROUNDING OF BERGER SPHERES. For A=B=x, C=z, the anisotropy
      u = x/z obeys the exact identity
        du/dt = -4 (x - z) / (x z),
      so sign(du/dt) = -sign(u - 1) EXACTLY: u moves monotonically toward
      1 from either side. The carrier rounds out — the shape-convergence
      half of the Ricci-flow picture, proved here by one rational identity.

 (T5) SEAM READOUT (thermodynamics/02 instance, same grammar as YM-6).
      lambda(t) = min spec(-4 Delta + R) is a bottom-of-spectrum threshold
      object and T2 makes it strictly increasing, so the threshold count
      k(t, mu) = #{ spec > mu } can only flow in ONE direction along the
      flow: seam crossings are irreversible. The YM line met the same
      grammar with a constant count (no flow); here the count flows
      monotonically. Same primitive, two different regimes.

THE HONEST REMAINDER — what this capsule does NOT touch, restated from the
adapter's obstruction tower:
  - surgery (neck/cap data at singular seams) — NOT addressed
  - non-collapsing / kappa-solutions / canonical neighbourhoods — NOT
  - singularity classification beyond this homogeneous carrier — NOT
  - extinction in finite time for general simply-connected 3-manifolds,
    topological reconstruction from the flow + surgery ledger — NOT
  - the Poincare conjecture — NOT CLAIMED, in any part.
Hamilton-Perelman remains a PINNED NAMED DEPENDENCY (CIRC-1 discipline:
classical results are cited as witnesses, never rederived natively). The
homogeneous carrier is exactly the case where the hard content (surgery,
collapse) is absent — that is why it is a calibration and not a proof.

Controls:
  C1  identity tamper: perturbing the R formula makes the T2 remainder
      nonzero (the identity is verified, not asserted).
  C2  flow tamper: perturbing the ODE right-hand side breaks T2.
  C3  round/general consistency: the general formulas specialize to the
      exact round solution at A=B=C.
  C4  strictness: |Ric|^2 > 0 at sampled rational metrics (equality would
      require Ric = 0).
  C5  collapse boundary recorded honestly: as C -> 0 with A=B fixed the
      carrier leaves the adapter's admissible region (collapse is an
      obstruction-tower entry, not a certified regime); sampled values are
      reported, no claim is made there.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(200000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import sympy as sp  # noqa: E402

_YM = os.path.join(os.path.dirname(HERE), "..",
                   "yang-mills-certified-benchmark", "certificates")
sys.path.insert(0, os.path.abspath(_YM))
from ym1_certified_gap import _dec, canonical_sha  # noqa: E402

ANCHOR = ("Hamilton-Perelman Ricci flow with surgery (2003) and Perelman's "
          "monotonicity of the lambda-functional: PINNED NAMED DEPENDENCY, "
          "cited as classical witness, never rederived natively (CIRC-1). "
          "This capsule reproduces the lambda-monotonicity statement "
          "exactly on the SU(2) homogeneous carrier only.")

OBSTRUCTION_TOWER_UNTOUCHED = [
    "surgery / neck-cap data at singular seams",
    "non-collapsing, kappa-solutions, canonical neighbourhoods",
    "singularity classification off the homogeneous carrier",
    "finite-time extinction for general simply-connected 3-manifolds",
    "topological reconstruction from flow + surgery ledger",
    "the Poincare conjecture itself",
]


# ------------------------------------------------- exact symbolic layer
def symbols_and_flow(tamper_R=False, tamper_flow=False):
    A, B, C = sp.symbols("A B C", positive=True)
    r1 = (A**2 - (B - C)**2) / (B * C)
    r2 = (B**2 - (C - A)**2) / (C * A)
    r3 = (C**2 - (A - B)**2) / (A * B)
    R = r1 / A + r2 / B + r3 / C
    if tamper_R:
        R = R + sp.Rational(1, 7) / (A * B * C)
    dA, dB, dC = -2 * r1, -2 * r2, -2 * r3
    if tamper_flow:
        dA = dA * sp.Rational(9, 10)
    return (A, B, C), (r1, r2, r3), R, (dA, dB, dC)


def monotonicity_remainder(tamper_R=False, tamper_flow=False):
    """dR/dt - 2|Ric|^2 as a simplified rational function (0 iff identity)."""
    (A, B, C), (r1, r2, r3), R, (dA, dB, dC) = symbols_and_flow(
        tamper_R, tamper_flow)
    dR = sp.diff(R, A) * dA + sp.diff(R, B) * dB + sp.diff(R, C) * dC
    ric2 = (r1 / A)**2 + (r2 / B)**2 + (r3 / C)**2
    return sp.simplify(sp.together(dR - 2 * ric2)), ric2


def R_closed_form():
    (A, B, C), _, R, _ = symbols_and_flow()
    return sp.factor(sp.simplify(R))


def berger_anisotropy_identity():
    """du/dt for u = x/z on A=B=x, C=z; returns simplified expression."""
    x, z = sp.symbols("x z", positive=True)
    r1 = (x**2 - (x - z)**2) / (x * z)
    r3 = z**2 / x**2
    dx, dz = -2 * r1, -2 * r3
    u = x / z
    return sp.factor(sp.simplify(sp.diff(u, x) * dx + sp.diff(u, z) * dz)), \
        (x, z)


# ------------------------------------------------- exact rational layer
def ric_eigs(A: F, B: F, C: F):
    return ((A * A - (B - C) ** 2) / (B * C),
            (B * B - (C - A) ** 2) / (C * A),
            (C * C - (A - B) ** 2) / (A * B))


def scalar_R(A: F, B: F, C: F) -> F:
    r1, r2, r3 = ric_eigs(A, B, C)
    return r1 / A + r2 / B + r3 / C


def ric_norm_sq(A: F, B: F, C: F) -> F:
    r1, r2, r3 = ric_eigs(A, B, C)
    return (r1 / A) ** 2 + (r2 / B) ** 2 + (r3 / C) ** 2


def round_solution(a0: F, t: F):
    """Exact round-sphere solution: a(t) = a0 - 2t, lambda = R = 3/a."""
    a = a0 - 2 * t
    return a, (F(3) / a if a != 0 else None)


def run():
    # T1 machine content: closed form of R and its trace consistency
    R_expr = R_closed_form()
    A, B, C = sp.symbols("A B C", positive=True)
    R_ref = (2 * (A * B + B * C + C * A) - (A**2 + B**2 + C**2)) / (A * B * C)
    t1 = sp.simplify(R_expr - R_ref) == 0

    # T2 exact identity + SOS witness
    rem, ric2 = monotonicity_remainder()
    t2 = sp.simplify(rem) == 0
    sos_terms = 3   # (r1/A)^2 + (r2/B)^2 + (r3/C)^2

    # C1 / C2 tampers must break it
    rem_tr, _ = monotonicity_remainder(tamper_R=True)
    rem_tf, _ = monotonicity_remainder(tamper_flow=True)
    c1 = sp.simplify(rem_tr) != 0
    c2 = sp.simplify(rem_tf) != 0

    # T3 exact round solution, sampled at rational times
    a0 = F(1)
    round_rows = {}
    prev = None
    round_monotone = True
    for t in [F(0), F(1, 10), F(1, 4), F(2, 5), F(49, 100)]:
        a, lam = round_solution(a0, t)
        round_rows[str(t)] = {"a": str(a), "lambda": _dec(lam, 20)}
        if prev is not None and not (lam > prev):
            round_monotone = False
        prev = lam
    t_star = a0 / 2
    # C3 consistency: general formulas at A=B=C=a reproduce R = 3/a
    c3 = all(scalar_R(a, a, a) == F(3) / a
             for a in [F(1), F(3, 7), F(11, 5)])

    # T4 Berger rounding identity and exact sign law
    du_expr, (xs, zs) = berger_anisotropy_identity()
    du_ref = -4 * (xs - zs) / (xs * zs)
    t4 = sp.simplify(du_expr - du_ref) == 0
    berger_rows = {}
    sign_law_ok = True
    for (x, z) in [(F(2), F(1)), (F(1), F(2)), (F(5, 4), F(1)),
                   (F(1), F(9, 10))]:
        u = x / z
        du = -4 * (x - z) / (x * z)
        berger_rows[f"u={u}"] = {"du_dt": _dec(du, 20),
                                 "sign_matches_-(u-1)":
                                     (du > 0) == (u < 1)}
        if (du > 0) != (u < 1):
            sign_law_ok = False

    # C4 strictness: |Ric|^2 > 0 on sampled metrics
    samples = [(F(1), F(1), F(1)), (F(2), F(1), F(1)), (F(3), F(2), F(1)),
               (F(1), F(1), F(1, 4))]
    c4 = all(ric_norm_sq(*s) > 0 for s in samples)
    lam_and_rate = {f"({a},{b},{c})": {"lambda": _dec(scalar_R(a, b, c), 15),
                                       "d_lambda_dt":
                                           _dec(2 * ric_norm_sq(a, b, c), 15)}
                    for (a, b, c) in samples}

    # C5 collapse boundary recorded, not claimed
    collapse_rows = {}
    for z in [F(1, 10), F(1, 100), F(1, 1000)]:
        collapse_rows[f"C={z}"] = {"lambda": _dec(scalar_R(F(1), F(1), z), 10),
                                   "status": "outside certified regime "
                                             "(collapse = obstruction-tower "
                                             "entry)"}

    ok = (t1 and t2 and c1 and c2 and c3 and t4 and sign_law_ok
          and round_monotone and c4)
    cert = {
        "certificate_type": "PC1_PERELMAN_LAMBDA_MONOTONICITY_SU2_CALIBRATION",
        "claim_status": "calibration_only_exact_on_homogeneous_carrier",
        "purpose": ("Calibration of the Poincare adapter against the accepted "
                    "Hamilton-Perelman architecture. NOT a proof, NOT new "
                    "mathematics, NO part of the Poincare conjecture "
                    "claimed."),
        "anchor_pinned": ANCHOR,
        "theorems": {
            "T1_lambda_equals_R_on_homogeneous":
                "lambda(g) = min spec(-4 Delta + R) = R on any homogeneous "
                "metric; R = [2(AB+BC+CA)-(A^2+B^2+C^2)]/(ABC) verified",
            "T2_exact_monotonicity_identity":
                "dR/dt - 2|Ric|^2 = 0 identically in Q(A,B,C); SOS witness "
                f"with {sos_terms} squares => d lambda/dt >= 0, strict here",
            "T3_exact_round_solution":
                f"a(t) = a0 - 2t, lambda = 3/a, extinction t* = a0/2 "
                f"= {t_star} exactly for a0 = {a0}",
            "T4_exact_berger_rounding":
                "du/dt = -4(x-z)/(xz) for u = x/z; sign(du/dt) = -sign(u-1) "
                "exactly => monotone rounding toward the round sphere",
            "T5_seam_readout":
                "lambda is a bottom-of-spectrum threshold object, strictly "
                "increasing => threshold crossings are irreversible "
                "(one-directional seam flow; thermodynamics/02 instance, "
                "same primitive as the YM seam-integer dock)",
        },
        "round_flow_samples": round_rows,
        "berger_anisotropy_samples": berger_rows,
        "lambda_and_rate_samples": lam_and_rate,
        "collapse_boundary_recorded": collapse_rows,
        "honest_remainder": {
            "obstruction_tower_untouched": OBSTRUCTION_TOWER_UNTOUCHED,
            "poincare_conjecture": "NOT CLAIMED",
            "note": ("The homogeneous carrier is precisely the case where "
                     "surgery and collapse are absent; that is why this is a "
                     "calibration, not a proof."),
        },
        "controls": {
            "C1_R_tamper_breaks_identity": bool(c1),
            "C2_flow_tamper_breaks_identity": bool(c2),
            "C3_round_specialization_consistent": bool(c3),
            "C4_strict_positivity_of_Ric_norm": bool(c4),
            "C5_collapse_boundary_recorded_not_claimed": True,
            "T1_R_closed_form_verified": bool(t1),
            "T2_identity_remainder_zero": bool(t2),
            "T4_berger_identity_verified": bool(t4),
            "round_lambda_strictly_increasing": bool(round_monotone),
            "berger_sign_law_exact": bool(sign_law_ok),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "PC1_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_PC1.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"])
    print("T2 remainder zero:", cert["controls"]["T2_identity_remainder_zero"])
    print("sha256:", sha)
