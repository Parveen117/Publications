"""LAM-2: Derived seam phase — the counterpoint to LAM1-F1.

LAM-1 certified that the legacy Lambda-FE phase e^{i phi} was DECLARED.
LAM-2 certifies the derived alternative and the impossibility half:

T1  GAUSS MAGNITUDE, EXACT. For every primitive Dirichlet character chi
    mod q on the grid, |tau(chi)|^2 = q as an exact identity in the
    cyclotomic ring Z[zeta_L] (cyclic-convolution arithmetic + reduction
    mod Phi_L; no floats, no complex numbers). Control: the imprimitive
    character mod 9 (induced from mod 3) gives |tau|^2 = 0, not 9.

T1b QUARTER-TURN INSTANCE. tau(chi_4)^2 = -4 exactly, i.e. the derived
    seam phase tau/sqrt(q) at q = 4 satisfies the cut relation
    iota^2 = -1. The derived phase lands on the quarter turn.

T2  THE DERIVATION STEP, EXACT. The separated-Gauss identity
        sum_{a mod q} chibar(a) e^{2 pi i a n / q} = chi(n) tau(chibar)
    verified for EVERY residue n and EVERY primitive chi on the grid —
    including n with gcd(n, q) > 1, where both sides vanish (that
    vanishing IS primitivity). This identity is precisely the step that
    converts a multiplicative twist into additive characters that the
    Poisson/theta seam can transport: it is where a genuine twist's FE
    phase is DERIVED. Control: for the imprimitive character mod 9 the
    identity FAILS at an explicit residue (exact nonzero witness).

T3  THE DERIVED PHASE, CERTIFIED ANALYTICALLY. The odd twisted theta
    seam at q = 4:  theta_chi(1/x) = (tau(chi)/(i sqrt(q))) x^{3/2}
    theta_chi(x) with derived multiplier exactly 1, certified at x = 2
    by two-sided directed rational enclosures: theta_chi(1/2) overlaps
    sqrt(8) * theta_chi(2), widths < 1e-20. Control: multiplier tamper
    (sqrt(2) in place of sqrt(8)) separates.

T4  NO CONDUCTOR FOR THE LEGACY TWIST. For every modulus q <= 60 there
    is an explicit congruent pair n == m (mod q) with
    Omega(n) != Omega(m): chi_phi(n) = e^{i phi Omega(n)} is periodic
    for NO modulus, hence is not a Dirichlet character of any
    conductor, hence the Gauss-seam derivation of T1/T2 is UNAVAILABLE
    to it. Finding LAM2-F3: the legacy declaration (LAM1-F1) was
    forced, not an oversight.

Claim boundary: all classical statements (Gauss sums, twisted theta
functional equation) are pinned classical anchors; this capsule
verifies exact finite instances and one certified analytic instance.
N1/N2/N3, K0/L0/RH, YM continuum gates remain OPEN. No claims about
zeros of any L-function. The native-N2 content established here is the
TEMPLATE (seam phase = finite arithmetic on the cut's residue classes,
derivable) plus the impossibility half for the legacy twist.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))


def _load_lam1():
    spec = importlib.util.spec_from_file_location(
        "lam1", os.path.join(HERE, "lam1_seam_interface.py"))
    mod = importlib.util.module_from_spec(spec)
    sys.modules.setdefault("lam1", mod)
    spec.loader.exec_module(mod)
    return mod


lam1 = _load_lam1()

# ----------------------------------------------------------------------
# Exact cyclotomic-ring arithmetic: elements of Z[zeta_L] as integer
# vectors in Z[x]/(x^L - 1); equality decided by reduction mod Phi_L.
# ----------------------------------------------------------------------

def vec_new(L):
    return [0] * L


def vec_add_monomial(v, exponent, coeff=1):
    v[exponent % len(v)] += coeff


def vec_mul(a, b):
    L = len(a)
    out = [0] * L
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    out[(i + j) % L] += x * y
    return out


def vec_conj(a):
    L = len(a)
    out = [0] * L
    for i, x in enumerate(a):
        out[(-i) % L] += x
    return out


def vec_sub(a, b):
    return [x - y for x, y in zip(a, b)]


def vec_scale(a, c):
    return [c * x for x in a]


def vec_shift(a, k):
    L = len(a)
    out = [0] * L
    for i, x in enumerate(a):
        out[(i + k) % L] += x
    return out


def reduce_elem(v):
    """Canonical form of v in Z[zeta_L]: remainder mod Phi_L."""
    L = len(v)
    _, r = lam1.poly_divmod_exact(v, lam1.cyclotomic(L))
    return lam1.poly_trim(r)


def elem_is_zero(v):
    return reduce_elem(v) == [0]


def elem_as_int(v):
    r = reduce_elem(v)
    assert len(r) == 1, "element is not a rational integer"
    return r[0]


# ----------------------------------------------------------------------
# Dirichlet characters on cyclic unit groups (prime q; q = 4; q = 9)
# ----------------------------------------------------------------------

def unit_group(q):
    return [a for a in range(1, q) if gcd(a, q) == 1]


def primitive_root(q):
    units = unit_group(q)
    order = len(units)
    for g in range(2, q):
        if gcd(g, q) != 1:
            continue
        seen, x = set(), 1
        for _ in range(order):
            x = (x * g) % q
            seen.add(x)
        if len(seen) == order:
            return g
    raise AssertionError(f"no primitive root mod {q}")


def dlog_table(q, g):
    table, x = {}, 1
    for k in range(len(unit_group(q))):
        table[x] = k
        x = (x * g) % q
    return table


class CyclicCharacter:
    """chi_j on the cyclic group (Z/q)^* of order m0, embedded in
    Z[zeta_L] with L = lcm(q, m0)."""

    def __init__(self, q, j):
        self.q = q
        self.j = j
        self.g = primitive_root(q)
        self.dlog = dlog_table(q, self.g)
        self.m0 = len(self.dlog)          # group order
        L = q * self.m0 // gcd(q, self.m0)
        self.L = L
        self.step_mult = L // self.m0     # zeta_{m0} = zeta_L^{step}
        self.step_add = L // q            # zeta_q = zeta_L^{step}

    def value_exponent(self, a):
        """chi(a) = zeta_L^{e}; None means chi(a) = 0."""
        a %= self.q
        if gcd(a, self.q) != 1:
            return None
        return (self.j * self.dlog[a] % self.m0) * self.step_mult

    def conj_value_exponent(self, a):
        e = self.value_exponent(a)
        return None if e is None else (-e) % self.L

    def order(self):
        return self.m0 // gcd(self.j, self.m0)

    def is_primitive(self):
        q = self.q
        if q in (3, 4, 5, 7, 11, 13):     # prime or q=4: nontrivial => primitive
            return self.j % self.m0 != 0
        if q == 9:                        # trivial on 1+3Z <=> j in {0, 3}
            return self.j % self.m0 not in (0, 3)
        raise AssertionError("grid modulus only")

    def gauss_sum_vec(self, conjugate=False):
        v = vec_new(self.L)
        for a in unit_group(self.q):
            e = (self.conj_value_exponent(a) if conjugate
                 else self.value_exponent(a))
            vec_add_monomial(v, e + a * self.step_add)
        return v


def grid_characters():
    out = []
    for q in (3, 4, 5, 7, 9, 11, 13):
        m0 = len(unit_group(q))
        for j in range(1, m0):
            out.append(CyclicCharacter(q, j))
    return out


# ----------------------------------------------------------------------
# T3: certified twisted theta seam at q = 4
# ----------------------------------------------------------------------

def chi4_sign(n):
    if n % 2 == 0:
        return 0
    return 1 if n % 4 == 1 else -1


def twisted_theta_iv(x, N=14):
    """Enclosure of theta_chi4(x) = 2 sum_{n>=1} chi4(n) n e^{-pi n^2 x/4},
    exact rational x > 0, directed rational intervals."""
    pi_lo, pi_hi = lam1.pi_brackets()
    s_lo = Fr(0)
    s_hi = Fr(0)
    for n in range(1, N + 1):
        sgn = chi4_sign(n)
        if sgn == 0:
            continue
        e_lo, e_hi = lam1.exp_neg_iv(pi_lo * n * n * x / 4,
                                     pi_hi * n * n * x / 4)
        t_lo, t_hi = n * e_lo, n * e_hi
        if sgn > 0:
            s_lo += t_lo
            s_hi += t_hi
        else:
            s_lo -= t_hi
            s_hi -= t_lo
    # tail: |n e^{-pi x n^2/4}| for n >= N+1; ratio of consecutive
    # magnitudes <= 2 e^{-pi x (2N+3)/4} =: r < 1/2 (asserted)
    head_hi = (N + 1) * lam1.exp_neg_brackets(pi_lo * x * (N + 1) ** 2 / 4)[1]
    r_hi = 2 * lam1.exp_neg_brackets(pi_lo * x * (2 * N + 3) / 4)[1]
    assert r_hi < Fr(1, 2)
    tail = head_hi / (1 - r_hi)
    lo = 2 * s_lo - 2 * tail
    hi = 2 * s_hi + 2 * tail
    return lam1.round_out(lo, hi)


# ----------------------------------------------------------------------
# Certificate build
# ----------------------------------------------------------------------

def build_certificate():
    cert = {}
    cert["certificate_type"] = "LAM2_DERIVED_SEAM_PHASE"
    cert["claim_status"] = (
        "derived-phase template certified on finite grid + one certified "
        "analytic seam; impossibility half certified for the legacy "
        "winding twist; no zero claims"
    )
    cert["anchor_pinned"] = (
        "Classical statements — Gauss sums of Dirichlet characters, "
        "separated-Gauss identity for primitive characters, odd twisted "
        "theta functional equation theta_chi(1/x) = (tau(chi)/(i sqrt q)) "
        "x^{3/2} theta_chibar(x): PINNED NAMED DEPENDENCIES, cited as "
        "classical witnesses (CIRC-1). This capsule verifies exact finite "
        "instances in Z[zeta_L] and one certified analytic instance."
    )

    chars = grid_characters()
    primitive = [c for c in chars if c.is_primitive()]
    imprimitive = [c for c in chars if not c.is_primitive()]
    assert len(imprimitive) == 1 and imprimitive[0].q == 9 \
        and imprimitive[0].j == 3
    chi9_impr = imprimitive[0]

    # ---------------- T1 ----------------
    t1_count = 0
    for ch in primitive:
        tau = ch.gauss_sum_vec()
        mag = elem_as_int(vec_mul(tau, vec_conj(tau)))
        assert mag == ch.q, f"|tau|^2 = {mag} != q = {ch.q} at (q,j)=" \
            f"({ch.q},{ch.j})"
        t1_count += 1
    tau_impr = chi9_impr.gauss_sum_vec()
    mag_impr = elem_as_int(vec_mul(tau_impr, vec_conj(tau_impr)))
    assert mag_impr == 0
    cert["T1_gauss_magnitude_exact"] = {
        "statement": (
            "|tau(chi)|^2 = q as an exact identity in Z[zeta_L] for every "
            "primitive character on the grid (cyclic convolution + "
            "reduction mod Phi_L; no floats, no complex numbers)"
        ),
        "grid_moduli": [3, 4, 5, 7, 9, 11, 13],
        "primitive_characters_checked": t1_count,
        "verdict": "PASS",
    }
    cert["controls_C1_imprimitive_magnitude"] = {
        "q": 9, "j": 3,
        "induced_from": "quadratic character mod 3",
        "tau_times_conj_tau": mag_impr,
        "expected_if_primitive": 9,
        "separated": True,
    }

    # ---------------- T1b ----------------
    chi_4 = next(c for c in primitive if c.q == 4)
    tau4 = chi_4.gauss_sum_vec()
    tau4_sq = elem_as_int(vec_mul(tau4, tau4))
    assert tau4_sq == -4
    cert["T1b_quarter_turn_instance"] = {
        "statement": (
            "tau(chi_4)^2 = -4 exactly: the derived seam phase "
            "tau/sqrt(q) at q=4 satisfies the cut relation iota^2 = -1 — "
            "the derived phase IS the quarter turn (remark-level native "
            "hook, exact)"
        ),
        "tau_squared": tau4_sq,
        "verdict": "PASS",
    }

    # ---------------- T2 ----------------
    t2_checks = 0
    for ch in primitive:
        tau_bar = ch.gauss_sum_vec(conjugate=True)
        for n in range(ch.q):
            # RHS of the derivation step: sum_a chibar(a) zeta_q^{an}
            rhs = vec_new(ch.L)
            for a in unit_group(ch.q):
                vec_add_monomial(
                    rhs, ch.conj_value_exponent(a) + a * n * ch.step_add)
            # LHS: chi(n) tau(chibar)
            e_n = ch.value_exponent(n)
            lhs = vec_new(ch.L) if e_n is None else vec_shift(tau_bar, e_n)
            assert elem_is_zero(vec_sub(rhs, lhs)), \
                f"derivation step fails at (q,j,n)=({ch.q},{ch.j},{n})"
            t2_checks += 1
    # control: imprimitive chi mod 9 fails at an explicit residue
    fail_witness = None
    tau_bar9 = chi9_impr.gauss_sum_vec(conjugate=True)
    for n in range(chi9_impr.q):
        rhs = vec_new(chi9_impr.L)
        for a in unit_group(9):
            vec_add_monomial(
                rhs, chi9_impr.conj_value_exponent(a)
                + a * n * chi9_impr.step_add)
        e_n = chi9_impr.value_exponent(n)
        lhs = (vec_new(chi9_impr.L) if e_n is None
               else vec_shift(tau_bar9, e_n))
        if not elem_is_zero(vec_sub(rhs, lhs)):
            fail_witness = n
            break
    assert fail_witness is not None, "imprimitive control failed to separate"
    cert["T2_derivation_step_exact"] = {
        "statement": (
            "Separated-Gauss identity sum_a chibar(a) zeta_q^{an} = "
            "chi(n) tau(chibar) verified exactly for EVERY residue n and "
            "every primitive chi on the grid, including gcd(n,q)>1 where "
            "both sides vanish (the vanishing IS primitivity). This is "
            "the exact step that converts a multiplicative twist into "
            "additive characters transportable by the Poisson/theta "
            "seam — where a genuine twist's FE phase is DERIVED"
        ),
        "identities_checked": t2_checks,
        "verdict": "PASS",
    }
    cert["controls_C2_imprimitive_derivation_fails"] = {
        "q": 9, "j": 3, "first_failing_residue_n": fail_witness,
        "separated": True,
    }

    # ---------------- T3 ----------------
    th_half = twisted_theta_iv(Fr(1, 2))
    th_two = twisted_theta_iv(Fr(2))
    assert th_half[0] > 0 and th_two[0] > 0
    s8 = lam1.sqrt_brackets(8)
    rhs = lam1.iv_mul_pos(s8, th_two)
    assert lam1.iv_overlap(th_half, rhs), "twisted theta seam separates"
    budget = Fr(1, 10 ** 20)
    assert lam1.iv_width(th_half) < budget
    assert lam1.iv_width(rhs) < budget
    s2 = lam1.sqrt_brackets(2)
    rhs_bad = lam1.iv_mul_pos(s2, th_two)
    assert not lam1.iv_overlap(th_half, rhs_bad), \
        "multiplier tamper failed to separate"
    gap = th_half[0] - rhs_bad[1]
    assert gap > 0
    cert["T3_derived_phase_certified_analytically"] = {
        "statement": (
            "Odd twisted theta seam at q=4 with DERIVED multiplier "
            "tau(chi_4)/(i sqrt 4) = 1: theta_chi(1/2) overlaps "
            "sqrt(8) * theta_chi(2), two-sided directed rational "
            "enclosures, no floats in any verdict path"
        ),
        "theta_half_enclosure": [lam1.dec(th_half[0]), lam1.dec(th_half[1])],
        "sqrt8_times_theta2_enclosure": [lam1.dec(rhs[0]), lam1.dec(rhs[1])],
        "width_theta_half": lam1.dec(lam1.iv_width(th_half), 40),
        "width_rhs": lam1.dec(lam1.iv_width(rhs), 40),
        "width_budget": "1e-20",
        "overlap": True,
        "verdict": "PASS",
    }
    cert["controls_C3_multiplier_tamper_separates"] = {
        "tamper": "sqrt(2) in place of sqrt(8)",
        "separation_gap_lower_bound": lam1.dec(gap, 20),
        "separated": True,
    }

    # ---------------- T4 ----------------
    witnesses = {}
    max_search = 0
    for q in range(2, 61):
        found = None
        for n in range(2, 4000):
            if lam1.omega_total(n) != lam1.omega_total(n + q):
                found = (n, n + q)
                break
        assert found is not None, f"no conductor witness found for q={q}"
        witnesses[q] = found
        max_search = max(max_search, found[0])
    sample = {str(q): list(witnesses[q]) for q in (2, 4, 12, 30, 60)}
    cert["T4_no_conductor_for_legacy_twist"] = {
        "statement": (
            "For every modulus q <= 60 there is an explicit congruent "
            "pair n == m (mod q) with Omega(n) != Omega(m): the legacy "
            "winding twist chi_phi(n) = e^{i phi Omega(n)} is periodic "
            "for NO modulus, hence is not a Dirichlet character of any "
            "conductor, hence the Gauss-seam derivation of T1/T2 is "
            "unavailable to it"
        ),
        "moduli_covered": "2..60",
        "max_witness_n": max_search,
        "sample_witnesses_n_m": sample,
        "verdict": "PASS",
    }

    cert["finding_LAM2_F1"] = (
        "Derived-phase template certified: for genuine (primitive-"
        "character) twists the functional-equation seam phase is finite "
        "arithmetic on the cut's residue classes — |tau|^2 = q exact "
        "(T1), the derivation step exact at every residue (T2), and the "
        "derived multiplier confirmed by a certified analytic seam (T3). "
        "A derived phase is load-bearing: tampering the multiplier "
        "separates certified enclosures."
    )
    cert["finding_LAM2_F2"] = (
        "tau(chi_4)^2 = -4 exactly: the derived seam phase at q=4 "
        "satisfies iota^2 = -1, the cut's defining relation. Recorded as "
        "an exact remark-level native hook, not a transport theorem."
    )
    cert["finding_LAM2_F3"] = (
        "The legacy winding twist has no conductor (explicit witness for "
        "every q <= 60), so no Gauss-seam derivation exists for it: the "
        "declared phase certified in LAM1-F1 was FORCED by the twist's "
        "own arithmetic, not an oversight. Together LAM1-F1 + LAM2-F1 "
        "state the N2 design constraint: a native functional equation "
        "must ride a twist whose phase is derivable as finite seam "
        "arithmetic."
    )
    cert["claim_boundary"] = {
        "classical_anchors": (
            "Gauss sums, separated-Gauss identity, twisted theta FE — "
            "pinned, finite instances verified, never claimed as new"
        ),
        "N1_native_continuation": "OPEN",
        "N2_native_functional_equation": "OPEN (template + impossibility "
                                         "half only)",
        "N3_identification": "OPEN",
        "K0_L0_RH": "OPEN",
        "YM_continuum_gates": "OPEN",
        "zeros_of_any_L_function": "no claims",
    }
    cert["arithmetic_discipline"] = (
        "exact integer cyclotomic-ring arithmetic (cyclic convolution in "
        "Z[x]/(x^L - 1), equality by reduction mod Phi_L) for all "
        "algebraic verdicts; directed rational intervals with 300-bit "
        "relative outward rounding for the analytic seam; no floating "
        "point in any verdict path"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "LAM2_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_LAM2.sha256"), "w") as f:
        f.write(digest + "\n")
    print("LAM2 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
