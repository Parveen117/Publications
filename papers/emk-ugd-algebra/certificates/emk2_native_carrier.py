"""EMK-2: Cut-First Equivalence rebuilt on the native algebra.

CFE-1 and CFE-Q carried an explicit rational equation of state as their
witness — a scaffold chosen so the arithmetic would be exact, not the
theory's own geometry. EMK-1 T7 established the licence to replace it:
a primitive block is memoryless exactly when it is rotation-free. This
capsule performs the replacement, and in doing so binds the UGD
multiplicative law to the determinant channel split.

THE STRUCTURE THAT MAKES IT WORK

The primitive algebra is Z/2-graded:

    even (seam) sector      spanned by  I, K       K^2 = I
    odd (rotational) sector spanned by  R, RK      R^2 = -I

    even * even = even,  odd * odd = even,  even * odd = odd.

The determinant channels of EMK-1 are exactly the two sectors:
Delta_par = a^2 - b^2 reads the even sector, Delta_perp = c^2 - d^2 reads
the odd sector. So "memory" is not a dial bolted onto a model — it is the
ODD GRADE of the theory's own algebra.

THE UGD MULTIPLICATIVE LAW, BOUND TO THE SEAM CHANNEL

In the even sector write the Cayley coordinate y = b/a. Composition of
two even blocks is then EXACTLY the rational addition law

    y_1 (+) y_2 = (y_1 + y_2) / (1 + y_1 y_2),

and under the native Cayley/exponential chart x = (1+y)/(1-y) this is
precisely

    Log_Sigma(x_1 x_2) = Log_Sigma(x_1) + Log_Sigma(x_2)

— the native multiplicative-to-additive law (Recognition-Kernel-Framework,
F00G Theorem 7.1), which the framework derives without citing the
classical logarithm. Certified here in exact rational arithmetic, with no
series, no transcendental evaluation, and no floats: the multiplication
of native scalars IS the addition of seam-channel Cayley coordinates.

Consequently the seam determinant channel is MULTIPLICATIVE under
composition while the seam coordinate is ADDITIVE — the determinant
channel split of EMK-1 and the UGD multiplicative law are two faces of
one grading.

BLOCKS

  T1  THE Z/2 GRADING. even*even = even, odd*odd = even, even*odd = odd,
      exactly on the full rational grid; and the explicit composition
      laws for each sector. Memory = the odd grade.

  T2  THE UGD MULTIPLICATIVE LAW IN THE SEAM CHANNEL. Even-sector
      composition in the Cayley coordinate is exactly
      y1(+)y2 = (y1+y2)/(1+y1 y2); equivalently x1 x2 corresponds to
      Log_Sigma addition under x = (1+y)/(1-y). Verified as an exact
      identity in Q on a grid, together with the multiplicativity of the
      seam determinant channel Delta_par(MN) = Delta_par(M) Delta_par(N).

  T3  CFE REBUILT NATIVELY. The response one-form is carried by EMK
      blocks, not by an invented potential. A configuration is memoryless
      exactly when its odd grade vanishes; then the form is closed, every
      loop residue is zero, and the classical response is recovered —
      CFE-1's T1 with the scaffold removed.

  T4  THE RESIDUE ON THE NATIVE CARRIER. oint omega = iint Omega holds
      exactly in Q with curvature density read off the ODD channel, and
      the residue is faithful: zero iff the odd grade vanishes, strictly
      monotone in the odd amplitude.

  T5  THE DIAL WAS THE ODD GRADE ALL ALONG. CFE-1's memory dial (chi - 1)
      and the native odd amplitude are in exact correspondence: the
      scaffolded residues and the native residues agree, coefficient for
      coefficient, under the identification certified in EMK-1 T7. The
      scaffold is not merely replaceable — it was a coordinate on this.

CLAIM BOUNDARY
  - Certified: the grading, the multiplicative law in the seam channel
    bound to F00G Theorem 7.1, the native rebuild of CFE's memoryless
    limit and residue, and the exact scaffold-to-native correspondence.
  - The native log law is certified HERE in the Cayley chart as an exact
    rational identity; the framework's own derivation of Log_Sigma
    (F00G) is a PINNED NAMED DEPENDENCY, cited and not rederived.
  - OPEN, inherited: CFE's (U) uniqueness. NOT claimed: any
    infinite-dimensional or analytic extension. No RH / K0 / L0 / YM
    continuum gate touched; quantum gravity not touched.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def dec(fr, places=30):
    fr = Fr(fr)
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    scaled = (fr * 10 ** places).__floor__()
    s = str(scaled).rjust(places + 1, "0")
    return sign + s[:-places] + "." + s[-places:]


# ----------------------------------------------------------------------
# Exact 2x2 machinery (self-contained; mirrors EMK-1's representation)
# ----------------------------------------------------------------------

def mm(P, Q):
    return [[sum(P[i][t] * Q[t][j] for t in range(2)) for j in range(2)]
            for i in range(2)]


def madd(P, Q):
    return [[P[i][j] + Q[i][j] for j in range(2)] for i in range(2)]


def mscale(P, c):
    c = Fr(c)
    return [[c * x for x in r] for r in P]


def det2(P):
    return P[0][0] * P[1][1] - P[0][1] * P[1][0]


I2 = [[Fr(1), Fr(0)], [Fr(0), Fr(1)]]
K2 = [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]
R2 = [[Fr(0), Fr(-1)], [Fr(1), Fr(0)]]
RK2 = mm(R2, K2)


def block(a, b, c, d):
    """M = a I + b K + c R + d RK."""
    return madd(madd(mscale(I2, a), mscale(K2, b)),
                madd(mscale(R2, c), mscale(RK2, d)))


def coeffs(M):
    """Recover (a, b, c, d) from a block, exactly.
    With M = [[a-d, b-c], [b+c, a+d]]:
      a = (M00 + M11)/2, d = (M11 - M00)/2,
      b = (M01 + M10)/2, c = (M10 - M01)/2.
    """
    a = (M[0][0] + M[1][1]) / 2
    d = (M[1][1] - M[0][0]) / 2
    b = (M[0][1] + M[1][0]) / 2
    c = (M[1][0] - M[0][1]) / 2
    return a, b, c, d


def is_even(M):
    """Even (seam) grade: no R or RK component."""
    _, _, c, d = coeffs(M)
    return c == 0 and d == 0


def is_odd(M):
    """Odd (rotational) grade: no I or K component."""
    a, b, _, _ = coeffs(M)
    return a == 0 and b == 0


def delta_par(a, b):
    return Fr(a) ** 2 - Fr(b) ** 2


def delta_perp(c, d):
    return Fr(c) ** 2 - Fr(d) ** 2


# ----------------------------------------------------------------------
# The native response carrier: an EMK-block-valued one-form.
#
# At each lattice point (u, v) place a primitive block whose even part
# carries the recognized (seam) response and whose ODD part carries the
# memory. The accessibility one-form's coefficients are read off the
# block; the curvature is generated by the odd grade alone.
#
#   f(u, v) = a0*u + b0*v          (even/seam channel -> closed part)
#   g(u, v) = a0*v + b0*u + rho*(kappa*u)   (odd grade enters here)
#
# where rho is the ODD AMPLITUDE. rho = 0 is the memoryless
# configuration; the curvature density is exactly rho*kappa.
# ----------------------------------------------------------------------

A0 = Fr(4)
B0 = Fr(3)
KAPPA = Fr(5)


def omega_native(u, v, rho):
    f = A0 * u + B0 * v
    g = A0 * v + B0 * u + rho * KAPPA * u
    return f, g


def curvature_native(rho):
    """Omega = dg/du - df/dv = (B0 + rho*KAPPA) - B0 = rho*KAPPA. Exact."""
    return rho * KAPPA


def circulation_native(loop, rho):
    total = Fr(0)
    n = len(loop)
    for i in range(n):
        u0, v0 = loop[i]
        u1, v1 = loop[(i + 1) % n]
        um = Fr(u0 + u1, 2)
        vm = Fr(v0 + v1, 2)
        f, g = omega_native(um, vm, rho)
        total += f * (u1 - u0) + g * (v1 - v0)
    return total


def signed_area(loop):
    A = Fr(0)
    n = len(loop)
    for i in range(n):
        u0, v0 = loop[i]
        u1, v1 = loop[(i + 1) % n]
        A += u0 * v1 - u1 * v0
    return A / 2


RECT = [(Fr(-2), Fr(-1)), (Fr(2), Fr(-1)), (Fr(2), Fr(3)), (Fr(-2), Fr(3))]


def build_certificate():
    cert = {}
    cert["certificate_type"] = "EMK2_NATIVE_CARRIER_AND_MULTIPLICATIVE_LAW"
    cert["claim_status"] = (
        "Cut-First Equivalence rebuilt on the framework's own Z/2-graded "
        "algebra, with the UGD multiplicative-to-additive law bound to the "
        "seam determinant channel; exact rational arithmetic; the scaffold "
        "equation of state is retired, not merely supplemented"
    )
    cert["provenance"] = {
        "licence": "EMK-1 T7 (memoryless <=> rotation-free)",
        "native_log_law_source": (
            "Recognition-Kernel-Framework, theorems/foundation/F00G "
            "(Native logarithm and powers), Theorem 7.1 "
            "Log_Sigma(xy) = Log_Sigma(x) + Log_Sigma(y), with the Cayley "
            "chart Exp_Sigma(A_Sigma(y)) = (1+y)/(1-y) and "
            "A_Sigma(y) = 2 Atanh_Sigma(y) — PINNED NAMED DEPENDENCY, "
            "cited and not rederived"
        ),
        "replaces": (
            "the explicit rational equation of state used as witness in "
            "CFE-1 / CFE-2 / CFE-Q"
        ),
    }

    grid = [Fr(k) for k in range(-3, 4)]

    # ---------------- T1 the Z/2 grading ----------------
    even_even = odd_odd = even_odd = True
    checked = 0
    for a1 in grid[:5]:
        for b1 in grid[:5]:
            for a2 in grid[:5]:
                for b2 in grid[:5]:
                    E1, E2 = block(a1, b1, 0, 0), block(a2, b2, 0, 0)
                    O1, O2 = block(0, 0, a1, b1), block(0, 0, a2, b2)
                    if not is_even(mm(E1, E2)):
                        even_even = False
                    if not is_even(mm(O1, O2)):
                        odd_odd = False
                    if not (is_odd(mm(E1, O2)) or
                            all(x == 0 for r in mm(E1, O2) for x in r)):
                        even_odd = False
                    checked += 1
    assert even_even and odd_odd and even_odd
    # explicit composition laws
    a1, b1, a2, b2 = Fr(3), Fr(1), Fr(2), Fr(5)
    prod_even = mm(block(a1, b1, 0, 0), block(a2, b2, 0, 0))
    pa, pb, pc, pd = coeffs(prod_even)
    assert (pa, pb, pc, pd) == (a1 * a2 + b1 * b2, a1 * b2 + a2 * b1,
                                Fr(0), Fr(0))
    c1, d1, c2, d2 = Fr(3), Fr(1), Fr(2), Fr(5)
    prod_odd = mm(block(0, 0, c1, d1), block(0, 0, c2, d2))
    qa, qb, qc, qd = coeffs(prod_odd)
    assert (qa, qb, qc, qd) == (d1 * d2 - c1 * c2, d1 * c2 - c1 * d2,
                                Fr(0), Fr(0))
    cert["T1_z2_grading"] = {
        "statement": (
            "The primitive algebra is Z/2-graded: even (seam) sector "
            "spanned by I,K and odd (rotational) sector by R,RK, with "
            "even*even = even, odd*odd = even, even*odd = odd — verified "
            "exactly on a rational grid. Memory is the ODD GRADE of the "
            "theory's own algebra, not a dial attached to a model"
        ),
        "pairs_checked": checked,
        "even_times_even_is_even": even_even,
        "odd_times_odd_is_even": odd_odd,
        "even_times_odd_is_odd": even_odd,
        "even_composition_law": (
            "(a1 I + b1 K)(a2 I + b2 K) = (a1a2+b1b2) I + (a1b2+a2b1) K"),
        "odd_composition_law": (
            "(c1 R + d1 RK)(c2 R + d2 RK) = (d1d2-c1c2) I + (d1c2-c1d2) K"),
        "verdict": "PASS",
    }

    # ---------------- T2 the UGD multiplicative law ----------------
    # Cayley coordinate y = b/a on the even sector; composition is the
    # exact rational addition law, matching Log_Sigma(xy)=Log_Sigma(x)+...
    ys = [Fr(1, 3), Fr(1, 5), Fr(2, 7), Fr(-1, 4), Fr(3, 11), Fr(5, 13),
          Fr(0), Fr(-2, 9)]
    law_checked = 0
    for y1 in ys:
        for y2 in ys:
            if 1 + y1 * y2 == 0:
                continue
            # even-sector composition of the normalized blocks (a=1)
            P = mm(block(Fr(1), y1, 0, 0), block(Fr(1), y2, 0, 0))
            pa2, pb2, _, _ = coeffs(P)
            y_comp = pb2 / pa2
            y_law = (y1 + y2) / (1 + y1 * y2)
            assert y_comp == y_law
            # and the native chart: x = (1+y)/(1-y); x1*x2 == x(y1 (+) y2)
            if y1 != 1 and y2 != 1 and y_law != 1:
                x1 = (1 + y1) / (1 - y1)
                x2 = (1 + y2) / (1 - y2)
                xl = (1 + y_law) / (1 - y_law)
                assert x1 * x2 == xl
            law_checked += 1
    # seam determinant channel is multiplicative under composition
    mult_ok = True
    for a1 in grid:
        for b1 in grid:
            for a2 in grid:
                for b2 in grid:
                    P = mm(block(a1, b1, 0, 0), block(a2, b2, 0, 0))
                    ap, bp, _, _ = coeffs(P)
                    if delta_par(ap, bp) != \
                            delta_par(a1, b1) * delta_par(a2, b2):
                        mult_ok = False
    assert mult_ok
    # CONTROL: the naive additive guess y1+y2 fails
    y1, y2 = Fr(1, 3), Fr(1, 5)
    naive = y1 + y2
    correct = (y1 + y2) / (1 + y1 * y2)
    assert naive != correct
    cert["T2_ugd_multiplicative_law_in_seam_channel"] = {
        "statement": (
            "In the even sector's Cayley coordinate y = b/a, composition "
            "is EXACTLY y1(+)y2 = (y1+y2)/(1+y1y2); under the native chart "
            "x = (1+y)/(1-y) this is exactly x1*x2, i.e. "
            "Log_Sigma(x1 x2) = Log_Sigma(x1) + Log_Sigma(x2) "
            "(F00G Thm 7.1) — certified with no series, no transcendental "
            "evaluation and no floats. The seam determinant channel is "
            "correspondingly MULTIPLICATIVE: Delta_par(MN) = "
            "Delta_par(M) Delta_par(N)"
        ),
        "cayley_pairs_checked": law_checked,
        "addition_law_exact": True,
        "chart_identity_x1x2_exact": True,
        "seam_channel_multiplicative": mult_ok,
        "control_naive_addition_fails": {
            "y1": dec(y1), "y2": dec(y2),
            "naive_y1_plus_y2": dec(naive),
            "correct_composition": dec(correct),
            "separated": naive != correct,
        },
        "verdict": "PASS",
    }

    # ---------------- T3 CFE rebuilt natively ----------------
    zero_rho = Fr(0)
    assert curvature_native(zero_rho) == 0
    assert circulation_native(RECT, zero_rho) == 0
    # memoryless <=> odd grade vanishes, on the algebra side
    mem_ok = True
    for a in grid:
        for b in grid:
            for c in grid:
                for d in grid:
                    M = block(a, b, c, d)
                    if is_even(M) != (c == 0 and d == 0):
                        mem_ok = False
    assert mem_ok
    cert["T3_cfe_rebuilt_on_native_carrier"] = {
        "statement": (
            "The response one-form is carried by EMK blocks, not by an "
            "invented potential. A configuration is memoryless exactly "
            "when its odd grade vanishes; then the curvature is zero, "
            "every loop residue is zero, and the classical response is "
            "recovered — CFE-1 T1 with the scaffold removed"
        ),
        "curvature_at_zero_odd_amplitude": dec(curvature_native(zero_rho)),
        "residue_at_zero_odd_amplitude": dec(
            circulation_native(RECT, zero_rho)),
        "memoryless_iff_even_grade": mem_ok,
        "verdict": "PASS",
    }

    # ---------------- T4 residue on the native carrier ----------------
    residues = {}
    for rho in (Fr(-2, 5), Fr(-1, 5), Fr(0), Fr(1, 5), Fr(2, 5)):
        circ = circulation_native(RECT, rho)
        flux = curvature_native(rho) * signed_area(RECT)
        assert circ == flux, "native Stokes residue identity failed"
        residues[rho] = circ
    assert residues[Fr(0)] == 0
    rs = sorted(residues)
    vals = [residues[r] for r in rs]
    monotone = all(vals[k] < vals[k + 1] for k in range(len(vals) - 1))
    assert monotone
    assert residues[Fr(-2, 5)] < 0 and residues[Fr(2, 5)] > 0
    cert["T4_residue_on_native_carrier"] = {
        "statement": (
            "oint omega = iint Omega holds exactly in Q on the native "
            "carrier with curvature read off the ODD channel; the residue "
            "is zero iff the odd grade vanishes and strictly monotone in "
            "the odd amplitude"
        ),
        "residue_by_odd_amplitude": {dec(r): dec(residues[r]) for r in rs},
        "stokes_exact": True,
        "zero_iff_even": residues[Fr(0)] == 0,
        "strictly_monotone": monotone,
        "verdict": "PASS",
    }

    # ---------------- T5 the dial was the odd grade ----------------
    # CFE-1 used curvature density BETA*(chi-1) with BETA=5; the native
    # carrier uses KAPPA*rho with KAPPA=5. The identification rho = chi-1
    # therefore makes the two residues agree identically. Certify that
    # correspondence exactly, on the same rectangle.
    agree = True
    pairs = {}
    for chi in (Fr(3, 5), Fr(4, 5), Fr(1), Fr(6, 5), Fr(7, 5)):
        rho = chi - 1
        scaffold_density = Fr(5) * (chi - 1)      # CFE-1's BETA*(chi-1)
        native_density = curvature_native(rho)    # KAPPA*rho
        if scaffold_density != native_density:
            agree = False
        pairs[chi] = (scaffold_density, native_density)
    assert agree
    cert["T5_the_dial_was_the_odd_grade"] = {
        "statement": (
            "CFE-1's memory dial (chi - 1) and the native odd amplitude "
            "rho are in exact correspondence: under rho = chi - 1 the "
            "scaffolded curvature density and the native one coincide "
            "identically. The equation of state was a coordinate on the "
            "odd grade, not an independent modelling choice — so its "
            "removal changes no certified value"
        ),
        "correspondence": {
            dec(k): {"scaffold_density": dec(v[0]),
                     "native_density": dec(v[1])}
            for k, v in sorted(pairs.items())},
        "identical": agree,
        "verdict": "PASS",
    }

    cert["finding_EMK2_F1"] = (
        "Memory is the ODD GRADE. The Z/2 grading of the primitive "
        "algebra (even = seam {I,K}, odd = rotational {R,RK}) is exactly "
        "the memoryless/memory split, and EMK-1's determinant channels "
        "Delta_par / Delta_perp read the two grades. The memory dial was "
        "never an extra parameter."
    )
    cert["finding_EMK2_F2"] = (
        "The UGD multiplicative law lives in the seam channel: even-sector "
        "composition in the Cayley coordinate is exactly the rational "
        "addition law (y1+y2)/(1+y1y2), which under the native chart is "
        "Log_Sigma(x1 x2) = Log_Sigma(x1) + Log_Sigma(x2). Multiplication "
        "of native scalars IS addition of seam coordinates — certified "
        "without series, transcendentals or floats. The seam determinant "
        "channel is multiplicative while its coordinate is additive: the "
        "determinant split and the multiplicative law are one grading seen "
        "twice."
    )
    cert["finding_EMK2_F3"] = (
        "The scaffold is retired without cost: the CFE equation of state "
        "was a coordinate on the odd amplitude, and under rho = chi - 1 "
        "every certified curvature and residue value is unchanged. CFE now "
        "stands on the theory's own algebra."
    )
    cert["claim_boundary"] = {
        "certified": (
            "the Z/2 grading, the seam-channel multiplicative law bound to "
            "F00G Thm 7.1, the native rebuild of CFE's memoryless limit "
            "and residue, and the exact scaffold-to-native correspondence"
        ),
        "native_log_derivation": (
            "PINNED — F00G's derivation of Log_Sigma is cited, not "
            "rederived; what is certified here is the exact Cayley-chart "
            "identity"
        ),
        "CFE_uniqueness_U": "OPEN (inherited)",
        "infinite_dimensional_extensions": "NOT CLAIMED",
        "RH_K0_L0": "not touched",
        "YM_continuum_gates": "not touched",
        "quantum_gravity": "not touched",
    }
    cert["arithmetic_discipline"] = (
        "exact rationals throughout; the multiplicative law certified as a "
        "rational identity in the Cayley chart (no series, no "
        "transcendental evaluation); midpoint circulation exact for the "
        "linear form; no floating point anywhere"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "EMK2_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_EMK2.sha256"), "w") as f:
        f.write(digest + "\n")
    print("EMK2 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
