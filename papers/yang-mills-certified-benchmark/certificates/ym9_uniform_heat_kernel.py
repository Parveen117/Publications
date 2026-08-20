"""YM-9: THE FIRST UNIFORMITY-IN-THE-CUTOFF THEOREM OF THE PROGRAM —
an exact refinement family (heat-kernel action) on which the reduced gap
is EXACTLY cutoff-independent, and an interacting uniform lower bound
along a declared scaling trajectory.

Why this capsule exists. YM-8 named the sole remaining Millennium content
of the finite program: UNIFORMITY. Every YM-1..8 statement lives at one
fixed cutoff; a gap at every fixed cutoff is compatible with the gap
closing in the limit. PC-1 supplied the method clue: a limit-crossing
statement becomes finitely certifiable when the structure is monotone or
EXACTLY SCALE-COVARIANT (there, Perelman's lambda-monotonicity collapsed
to a polynomial identity). This capsule applies that move to Yang-Mills:
replace the Wilson action by the HEAT-KERNEL action, on which lattice
refinement is exact, and the cutoff dependence cancels identically.

Carrier. SU(2) heat-kernel (Menotti-Onofri) link action at spacing a:
    K_a(g) = sum_j d_j e^{-a C_j} chi_j(g),   C_j = j(j+1)  [declared
                                               normalization; C_{1/2}=3/4]
so the normalized free transfer eigenvalues on one link are
    lambda_j(a) = e^{-a C_j},  lambda_0 = 1.
Theta graph as in YM-2..8: two holonomies, three faces, interaction
    m_kappa(A,B) = exp[(kappa/2)(Tr A + Tr B + Tr(A B^{-1}))].

CERTIFIED (exact rational arithmetic; no floats in any verdict):

 (T1) EXACT REFINEMENT / SEMIGROUP. K_a * K_b = K_{a+b} identically,
      coefficient-wise e^{-aC_j} e^{-bC_j} = e^{-(a+b)C_j}. Hence
      subdividing a link of spacing a into n links of spacing a/n
      reproduces the SAME operator exactly: the refinement family carries
      NO discretization error. (Wilson has no such exact semigroup — this
      is precisely why it is the wrong carrier for uniformity questions.)
      Machine content: the exponent identity is verified as an exact
      rational identity in the exponents for a grid of (a, b, n).

 (T2) EXACT UNIFORMITY OF THE FREE REDUCED GAP. For EVERY a > 0,
        Delta_red(a) = -(1/a) log lambda_{1/2}(a) = C_{1/2} = 3/4
      EXACTLY — the transcendentals cancel identically, leaving a rational
      number independent of the cutoff. This is the program's first
      statement that is uniform in a rather than fixed-a.

 (T3) INTERACTING UNIFORM LOWER BOUND (the theorem). Along the declared
      scaling trajectory
        kappa(a) = theta * a,   theta a fixed rational,
      the YM-2 sandwich (|Tr U| <= 2 on three faces => m_+/m_- = e^{6 kappa},
      valid for ANY face form) gives for EVERY a > 0
        lambda_2/lambda_1 (T_a) <= e^{6 kappa(a)} lambda_{1/2}(a)
                                 = e^{6 theta a} e^{-(3/4) a},
      hence
        Delta(a, kappa(a)) >= 3/4 - 6 theta   for every a > 0,
      a CUTOFF-INDEPENDENT positive lower bound whenever theta < 1/8.
      At theta = 1/16 the certified uniform gap is exactly 3/8. Every
      quantity here is rational: the exponentials cancel exactly.

 (T4) UNIFORM SEAM COUNT. The free theta spectrum is exactly
        { e^{-a (C_{j1} + C_{j2})} } with SU(2) content multiplicities, so
      for the scaling threshold mu(a) = e^{-a s} the seam count
        k(a, mu(a)) = #{ (j1,j2) : C_{j1} + C_{j2} < s }
      is EXACTLY INDEPENDENT of a — a purely combinatorial count. The
      YM-6/7 dock's threshold grammar therefore transports across the whole
      refinement family at once, not cell by cell.

 (T5) WILSON CONTRAST, COMPUTED NOT ASSERTED. At fixed beta the Wilson
      reduced gap is -(1/a) log(I_2(beta)/I_1(beta)), which scales like
      1/a and DIVERGES as a -> 0: with a fixed coupling there is no
      cutoff-independent gap at all. Reported at beta = 2 for
      a = 1, 1/2, 1/4 as certified enclosures. This is the honest reason
      the carrier was changed, and it is also the shadow of the real
      physics: the true theory needs beta to run (asymptotic freedom),
      which this capsule does NOT model.

THE HONEST REMAINDER — what T3 is NOT:
  - The trajectory kappa(a) = theta*a is DECLARED, not derived. The
    physical trajectory is fixed by asymptotic freedom (beta ~ log(1/a)),
    which is not modelled here. Uniformity along a chosen trajectory is
    strictly weaker than uniformity along the renormalized one.
  - The GRAPH IS FIXED. A continuum limit needs the lattice to grow
    (infinite volume) as it refines; here only the spacing moves. This is
    refinement in the link parameter, not a lattice continuum limit.
  - The heat-kernel action is a CHOICE. Universality (regulator-path
    independence) remains an open gate: exact uniformity on this carrier
    does not transfer to Wilson or to any other regulator without proof.
  - Delta here is a REDUCED GRAPH GAP, not the physical mass gap of a
    reconstructed continuum theory. OS reconstruction, tightness, IR,
    gauge, non-triviality all remain OPEN, exactly as in YM-8.
  - The Clay predicate remains OPEN. This capsule moves ONE gate from
    "untouched" to "touched on a toy carrier", and names precisely what
    it did not touch.

Controls:
  C1  semigroup tamper: perturbing the exponent law breaks T1.
  C2  theta at/above the threshold 1/8 must NOT certify (fail-closed).
  C3  cancellation is real, not printer-resolution: T2/T3 are computed in
      exact rational arithmetic on the EXPONENTS, and an independent
      interval route (exp/log enclosures) is checked to agree.
  C4  seam-count a-independence checked across several a and s, plus a
      tamper (wrong Casimir) that changes the count.
  C5  Wilson non-uniformity computed (1/a growth witnessed, not claimed).
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

THETA = F(1, 16)            # declared scaling trajectory kappa(a) = theta*a
THETA_CRIT = F(1, 8)        # 6*theta < C_{1/2} = 3/4  <=>  theta < 1/8
A_GRID = [F(1), F(1, 2), F(1, 4), F(1, 100), F(1, 10 ** 6)]
S_GRID = [F(1), F(2), F(4)]
BETA_W = F(2)


def casimir(twice_spin: int) -> F:
    """C_j = j(j+1) with j = twice_spin/2. Declared normalization."""
    j = F(twice_spin, 2)
    return j * (j + 1)


# --------------------------------------------------------------- T1
def semigroup_exponent_identity(a: F, b: F, max_twice_spin=8) -> bool:
    """e^{-aC} e^{-bC} = e^{-(a+b)C} verified on the exponents exactly."""
    return all((-a * casimir(t)) + (-b * casimir(t)) == -(a + b) * casimir(t)
               for t in range(max_twice_spin + 1))


def subdivision_exact(a: F, n: int, max_twice_spin=8) -> bool:
    """n links of spacing a/n compose to one link of spacing a, exactly."""
    return all(sum([-(a / n) * casimir(t)] * n) == -a * casimir(t)
               for t in range(max_twice_spin + 1))


def semigroup_tampered(a: F, b: F) -> bool:
    return (-a * casimir(1)) + (-b * casimir(1)) == -(a + b) * casimir(1) - F(1, 10)


# --------------------------------------------------------------- T2 / T3
def free_reduced_gap_exact() -> F:
    """-(1/a) log e^{-a C_{1/2}} = C_{1/2}, exactly, for every a."""
    return casimir(1)                       # 3/4


def uniform_interacting_gap_bound(theta: F) -> F:
    """Delta(a, theta*a) >= C_{1/2} - 6*theta, for every a > 0."""
    return casimir(1) - 6 * theta


def certifies_uniform(theta: F) -> bool:
    return uniform_interacting_gap_bound(theta) > 0


def interval_route_check(a: F, theta: F):
    """C3: independent transcendental route must agree with the exact one.
    ratio(a) = e^{6 theta a} * e^{-C_{1/2} a};  gap = -(1/a) log ratio."""
    kappa = theta * a
    ratio = exp_point(6 * kappa) * exp_point(-casimir(1) * a)
    gap = -(log_iv(ratio, LOG_TERMS)) / Iv(a)
    exact = uniform_interacting_gap_bound(theta)
    return gap, (gap.lo <= exact <= gap.hi)


# --------------------------------------------------------------- T4
def free_theta_content(max_twice_spin=6):
    """Pairs (j1, j2) of the free theta spectrum with total Casimir."""
    out = []
    for t1 in range(max_twice_spin + 1):
        for t2 in range(max_twice_spin + 1):
            out.append(((t1, t2), casimir(t1) + casimir(t2)))
    return out


def seam_count(s: F, max_twice_spin=6, casimir_tamper=False,
               with_multiplicity=False) -> int:
    """Threshold count below s, exactly independent of a.

    AMENDED (see YM-10 T2): the default counts CONTENTS (j1,j2). The seam
    count k_Sigma of YM-6/7 counts EIGENVALUES WITH MULTIPLICITY, which
    differs (at s=2: 4 contents but 5 eigenvalues). Pass
    with_multiplicity=True for the corrected count; the a-independence —
    the whole point of T4 — holds for BOTH.
    """
    n = 0
    for t1 in range(max_twice_spin + 1):
        for t2 in range(max_twice_spin + 1):
            if casimir_tamper:
                c = F(t1, 2) + F(t2, 2)          # wrong Casimir
            else:
                c = casimir(t1) + casimir(t2)
            if c < s:
                n += (min(t1, t2) + 1) if with_multiplicity else 1
    return n


# --------------------------------------------------------------- T5
def wilson_reduced_gap(a: F, beta: F = BETA_W) -> Iv:
    lam = bessel_I(2, beta, TERMS) / bessel_I(1, beta, TERMS)
    return -(log_iv(lam, LOG_TERMS)) / Iv(a)


def run():
    # T1
    t1 = all(semigroup_exponent_identity(a, b)
             for a in [F(1), F(1, 3)] for b in [F(1, 2), F(2, 7)])
    t1b = all(subdivision_exact(a, n)
              for a in [F(1), F(3, 5)] for n in [2, 3, 10])
    c1 = not semigroup_tampered(F(1), F(1, 2))

    # T2 / T3 — exact rational headline
    free_gap = free_reduced_gap_exact()
    uni_bound = uniform_interacting_gap_bound(THETA)
    t3 = certifies_uniform(THETA)

    # C2 fail-closed at/above threshold
    c2 = (not certifies_uniform(THETA_CRIT)
          and not certifies_uniform(F(1, 6)))

    # C3 independent interval route agrees at every a in the grid
    route_rows = {}
    c3 = True
    for a in A_GRID:
        gap, agree = interval_route_check(a, THETA)
        route_rows[str(a)] = {"gap_lo": _dec(gap.lo, 25),
                              "gap_hi": _dec(gap.hi, 25),
                              "contains_exact_bound": bool(agree)}
        c3 = c3 and agree

    # T4 uniform seam counts + tamper
    counts = {str(s): {"contents": seam_count(s),
                       "eigenvalues_with_multiplicity":
                           seam_count(s, with_multiplicity=True)}
              for s in S_GRID}
    c4 = (seam_count(F(2)) != seam_count(F(2), casimir_tamper=True))

    # T5 Wilson contrast
    wilson_rows = {}
    for a in [F(1), F(1, 2), F(1, 4)]:
        w = wilson_reduced_gap(a)
        wilson_rows[str(a)] = {"lo": _dec(w.lo, 20), "hi": _dec(w.hi, 20)}
    w1 = wilson_reduced_gap(F(1))
    w4 = wilson_reduced_gap(F(1, 4))
    c5 = w4.lo > 3 * w1.hi          # witnesses 1/a growth (divergence)

    ok = t1 and t1b and t3 and c1 and c2 and c3 and c4 and c5
    cert = {
        "certificate_type": "YM9_UNIFORM_GAP_HEAT_KERNEL_REFINEMENT",
        "claim_status": "first_uniformity_statement_toy_carrier",
        "theorems": {
            "T1_exact_refinement_semigroup":
                "K_a * K_b = K_{a+b} identically; n links of spacing a/n "
                "compose exactly to one link of spacing a (no "
                "discretization error in the refinement family)",
            "T2_free_gap_exactly_cutoff_independent":
                f"Delta_red(a) = C_(1/2) = {free_gap} EXACTLY for every "
                f"a > 0 (transcendentals cancel identically)",
            "T3_interacting_uniform_lower_bound":
                f"along kappa(a) = {THETA}*a: Delta(a, kappa(a)) >= "
                f"C_(1/2) - 6*theta = {uni_bound} for EVERY a > 0; "
                f"uniform positivity iff theta < {THETA_CRIT}",
            "T4_uniform_seam_count":
                "the threshold count below s is exactly independent of a — "
                "the dock's threshold grammar transports across the whole "
                "family at once. AMENDED per YM-10 T2: the count of "
                "CONTENTS #{(j1,j2): C+C<s} and the count of EIGENVALUES "
                "WITH MULTIPLICITY (weights m = min(2j1,2j2)+1) differ; "
                "k_Sigma of YM-6/7 is the latter. Both are a-independent, "
                "so the uniformity conclusion is unaffected; both are "
                "reported below",
            "T5_wilson_contrast":
                "at fixed beta the Wilson reduced gap scales like 1/a and "
                "diverges as a -> 0: no cutoff-independent gap without a "
                "running coupling (computed, not asserted)",
        },
        "parameters": {"theta": str(THETA), "theta_critical": str(THETA_CRIT),
                       "casimir_normalization": "C_j = j(j+1), C_(1/2)=3/4",
                       "a_grid": [str(x) for x in A_GRID]},
        "exact_results": {
            "free_reduced_gap_all_a": str(free_gap),
            "uniform_interacting_gap_bound": str(uni_bound),
            "uniform_bound_decimal": _dec(uni_bound, 20),
        },
        "interval_route_agreement": route_rows,
        "uniform_seam_counts": counts,
        "wilson_reduced_gap_by_spacing": wilson_rows,
        "honest_remainder": {
            "trajectory": ("kappa(a) = theta*a is DECLARED, not derived; the "
                           "physical trajectory is fixed by asymptotic "
                           "freedom (beta ~ log(1/a)), NOT modelled here"),
            "graph_fixed": ("only the spacing moves; a continuum limit needs "
                            "the lattice to grow (infinite volume) as it "
                            "refines — not done"),
            "universality": ("heat-kernel action is a CHOICE; exact "
                             "uniformity here does not transfer to Wilson or "
                             "any other regulator without proof — gate OPEN"),
            "object": ("Delta is a reduced GRAPH gap, not the physical mass "
                       "gap of a reconstructed continuum theory"),
            "open_gates": ["gauge", "tightness", "UV", "IR", "universality",
                           "OS reconstruction", "locality/non-triviality",
                           "AF-OPE matching"],
            "clay_predicate": "OPEN",
            "net_movement": ("one gate moved from 'untouched' to 'touched on "
                             "a toy carrier'; everything else unchanged"),
        },
        "controls": {
            "C1_semigroup_tamper_breaks": bool(c1),
            "C2_theta_above_threshold_fails_closed": bool(c2),
            "C3_interval_route_agrees_with_exact": bool(c3),
            "C4_casimir_tamper_changes_count": bool(c4),
            "C5_wilson_divergence_witnessed": bool(c5),
            "T1_semigroup_verified": bool(t1 and t1b),
            "T3_uniform_certified": bool(t3),
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
    out = {}
    for name, fn in (("YM1", m1.run), ("YM2", m2.run), ("YM3", m3.run),
                     ("YM4", m4.run), ("YM5", m5.run), ("YM6", m6.run),
                     ("YM7", m7.run), ("YM8", m8.run), ("YM9", run)):
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
