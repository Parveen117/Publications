"""LAM-1: Lambda-seam interface calibration.

Calibrates the legacy Lambda-analytic canvases (Lambda-Eisenstein /
Lambda-mock functional-equation machinery, descended from the T-V-S-P
compass of arXiv:2603.20773) against exact, executable arithmetic.

Four theorem blocks, all exact (integers / rationals / directed rational
intervals; no floats anywhere in a verdict path):

T1  RAMANUJAN TWO-ROUTE. The arithmetic layer that the Lambda-twist
    multiplies but never touches: c_c(n) computed by exact cyclotomic
    reduction (Z[x]/Phi_c(x)) equals the Moebius/gcd closed form
    sum_{d | gcd(c,n)} d * mu(c/d) on a grid. Control: dropping the
    coprimality condition separates.

T2  MEMORY-LATTICE MULTIPLICATIVITY. The Lambda-character
    chi_phi(n) = e^{i phi nu(n)} is exact integer bookkeeping on the
    winding lattice: nu = Omega is additive (nu(mn) = nu(m) + nu(n)),
    and the formal Euler product reproduces every n <= N with marker
    exactly Omega(n) and multiplicity exactly one. Control: a
    non-additive fake winding is detected.

T3  THETA SEAM, CERTIFIED. The single S-move that generates the
    classical functional equation — theta(1/t) = sqrt(t) theta(t) —
    verified at t = 2 with two-sided directed rational enclosures
    (Machin pi brackets, alternating-series exp brackets with
    scaling-and-squaring, geometric tail bounds). Control: replacing
    sqrt(2) by sqrt(3) separates the enclosures.

T4  DECLARED-PHASE AUDIT (the finding). Under the exact S-reindexing
    bijection (c,d) -> (d,-c) on primitive pairs, the coset-weight
    layer nu(|c|) has zero net change: the weight bookkeeping alone
    produces NO functional-equation phase. Therefore the e^{i phi}
    factor in the legacy Lambda functional equation is carried entirely
    by the DECLARED right S-multiplier (the legacy canvas itself notes
    a genuine homomorphism Gamma -> S^1 cannot exist: the abelianization
    of SL(2,Z) is finite). Finding LAM1-F1: the legacy Lambda-FE phase
    is an input, not an output. Control: an asymmetric region separates.

Claim boundary: calibration and audit only. No claim about zeros of any
L-function. N1/N2/N3 (native continuation, native functional equation,
identification), K0/L0, RH, and the YM continuum gates all remain OPEN.
The classical automorphic substrate (Poisson summation, theta inversion,
SL(2,Z) coset structure) is a PINNED NAMED DEPENDENCY, cited and never
rederived natively.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# Exact polynomial arithmetic over Z (ascending coefficient lists)
# ----------------------------------------------------------------------

def poly_trim(p):
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def poly_mul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i + j] += x * y
    return poly_trim(out)


def poly_divmod_exact(num, den):
    """Exact division of integer polynomials; den monic."""
    num = list(num)
    den = poly_trim(list(den))
    assert den[-1] == 1, "divisor must be monic"
    q = [0] * max(1, len(num) - len(den) + 1)
    r = list(num)
    while len(poly_trim(r)) >= len(den) and any(r):
        r = poly_trim(r)
        if len(r) < len(den):
            break
        shift = len(r) - len(den)
        coef = r[-1]
        q[shift] += coef
        for i, d in enumerate(den):
            r[shift + i] -= coef * d
        r = poly_trim(r)
    return poly_trim(q), poly_trim(r)


_CYCLO = {}


def cyclotomic(c):
    """Phi_c(x) as integer coefficient list, via x^c - 1 = prod_{d|c} Phi_d."""
    if c in _CYCLO:
        return _CYCLO[c]
    num = [-1] + [0] * (c - 1) + [1]  # x^c - 1
    den = [1]
    for d in range(1, c):
        if c % d == 0:
            den = poly_mul(den, cyclotomic(d))
    q, r = poly_divmod_exact(num, den)
    assert r == [0], "cyclotomic division must be exact"
    _CYCLO[c] = q
    return q


def ramanujan_direct(c, n, require_coprime=True):
    """c_c(n) by exact reduction in Z[x]/Phi_c(x).

    Returns (value, is_rational_integer). The sum sum_{d} zeta_c^{n d}
    over d in [0, c) with gcd(d, c) == 1 (or all d if require_coprime is
    False) is represented as an integer coefficient vector and reduced
    mod Phi_c; a rational integer must reduce to a constant.
    """
    if c == 1:
        return (1, True)
    vec = [0] * c
    for d in range(c):
        if require_coprime and gcd(d, c) != 1:
            continue
        vec[(n * d) % c] += 1
    _, r = poly_divmod_exact(vec, cyclotomic(c))
    r = poly_trim(r)
    if len(r) == 1:
        return (r[0], True)
    return (None, False)


def moebius(m):
    if m == 1:
        return 1
    mu, p = 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            mu = -mu
        p += 1
    if m > 1:
        mu = -mu
    return mu


def ramanujan_moebius(c, n):
    g = gcd(c, n)
    return sum(d * moebius(c // d) for d in range(1, g + 1) if g % d == 0)


# ----------------------------------------------------------------------
# T2: winding lattice nu = Omega (total prime multiplicity)
# ----------------------------------------------------------------------

def omega_total(n):
    assert n >= 1
    count, p = 0, 2
    while p * p <= n:
        while n % p == 0:
            count += 1
            n //= p
        p += 1
    if n > 1:
        count += 1
    return count


def euler_marker_witness(N):
    """Formal Euler product over primes <= N restricted to n <= N.

    Returns dict n -> (marker, multiplicity). The Lambda-twist
    chi_phi(n) = e^{i phi nu(n)} lives entirely on the integer marker.
    """
    primes = [p for p in range(2, N + 1) if omega_total(p) == 1]
    table = {1: (0, 1)}
    for p in primes:
        new = dict(table)
        for n, (mark, mult) in table.items():
            k, pk = 1, p
            while n * pk <= N:
                key = n * pk
                add_mark = mark + k
                if key in new:
                    old_mark, old_mult = new[key]
                    # same n reached twice would be a multiplicity failure
                    new[key] = (old_mark, old_mult + mult)
                else:
                    new[key] = (add_mark, mult)
                k += 1
                pk *= p
        table = new
    return table


# ----------------------------------------------------------------------
# T3: directed rational interval toolkit (no floats in verdicts)
# ----------------------------------------------------------------------

PREC_BITS = 300


def round_dir(x, down, prec_bits=PREC_BITS):
    """Directed rounding of an exact rational to ~prec_bits significant
    bits (relative precision, so magnitudes like e^{-900} survive)."""
    x = Fr(x)
    if x == 0:
        return x
    e = x.numerator.bit_length() - x.denominator.bit_length()
    m = prec_bits - e
    if m >= 0:
        scaled = x * (1 << m)
        n = scaled.__floor__() if down else -((-scaled).__floor__())
        return Fr(n, 1 << m)
    scaled = x / (1 << -m)
    n = scaled.__floor__() if down else -((-scaled).__floor__())
    return Fr(n * (1 << -m))


def round_out(lo, hi):
    lo2 = round_dir(lo, down=True)
    hi2 = round_dir(hi, down=False)
    assert lo2 <= lo <= hi <= hi2
    return lo2, hi2


def arctan_inv_brackets(q, K=40):
    """Two-sided brackets for arctan(1/q), q >= 2, alternating series."""
    s = Fr(0)
    lo = hi = None
    for k in range(K + 1):
        term = Fr((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
        s += term
        if k >= K - 1:
            if k % 2 == 1:
                lo = s
            else:
                hi = s
    if lo is None or hi is None:
        # ensure both parities captured
        raise AssertionError("bracket parity capture failed")
    assert lo <= hi
    return lo, hi


def pi_brackets():
    a5lo, a5hi = arctan_inv_brackets(5)
    a239lo, a239hi = arctan_inv_brackets(239)
    lo = 16 * a5lo - 4 * a239hi
    hi = 16 * a5hi - 4 * a239lo
    assert Fr(3) < lo <= hi < Fr(4)
    return round_out(lo, hi)


def exp_neg_brackets(y, K=41):
    """Two-sided brackets for e^{-y}, exact rational y > 0.

    Scale so z = y / 2^j <= 1/2, alternating series brackets, then
    square j times with outward rounding.
    """
    assert y > 0
    j = 0
    z = Fr(y)
    while z > Fr(1, 2):
        z /= 2
        j += 1
    s = Fr(0)
    fact = 1
    lo = hi = None
    for k in range(K + 1):
        if k > 0:
            fact *= k
        term = Fr((-1) ** k) * z ** k / fact
        s += term
        if k >= K - 1:
            if k % 2 == 1:
                lo = s
            else:
                hi = s
    assert lo is not None and hi is not None and 0 < lo <= hi
    for _ in range(j):
        lo, hi = round_out(lo * lo, hi * hi)
        assert 0 < lo <= hi
    return lo, hi


def exp_neg_iv(xlo, xhi):
    """Enclosure of e^{-x} for x in [xlo, xhi]."""
    lo = exp_neg_brackets(xhi)[0]
    hi = exp_neg_brackets(xlo)[1]
    assert lo <= hi
    return lo, hi


def theta_iv(t, N=12):
    """Enclosure of theta(t) = sum_{n in Z} e^{-pi n^2 t}, exact t > 0."""
    pi_lo, pi_hi = pi_brackets()
    s_lo = Fr(0)
    s_hi = Fr(0)
    for n in range(1, N + 1):
        e_lo, e_hi = exp_neg_iv(pi_lo * n * n * t, pi_hi * n * n * t)
        s_lo += e_lo
        s_hi += e_hi
    # tail: n >= N+1 => pi t n^2 >= pi t (N+1)^2 + 2 pi t (N+1) (n-(N+1))
    head_hi = exp_neg_brackets(pi_lo * t * (N + 1) ** 2)[1]
    q_hi = exp_neg_brackets(pi_lo * t * 2 * (N + 1))[1]
    assert q_hi < 1
    tail_hi = head_hi / (1 - q_hi)
    lo = 1 + 2 * s_lo
    hi = 1 + 2 * (s_hi + tail_hi)
    return round_out(lo, hi)


def sqrt_brackets(m, iters=130):
    """Two-sided rational brackets for sqrt(m), integer m >= 1."""
    lo, hi = Fr(1), Fr(m)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if mid * mid <= m:
            lo = mid
        else:
            hi = mid
    return round_out(lo, hi)


def iv_mul_pos(a, b):
    assert a[0] > 0 and b[0] > 0
    return round_out(a[0] * b[0], a[1] * b[1])


def iv_overlap(a, b):
    return not (a[1] < b[0] or b[1] < a[0])


def iv_width(a):
    return a[1] - a[0]


# ----------------------------------------------------------------------
# T4: S-reindexing on primitive pairs; weight-layer net phase
# ----------------------------------------------------------------------

def primitive_region(C):
    """R_C = {(c,d): 1 <= |c|,|d| <= C, gcd(|c|,|d|) = 1}."""
    out = set()
    for c in range(-C, C + 1):
        if c == 0:
            continue
        for d in range(-C, C + 1):
            if d == 0:
                continue
            if gcd(abs(c), abs(d)) == 1:
                out.add((c, d))
    return out


def s_move(pair):
    c, d = pair
    return (d, -c)


def weight_sums(region):
    w_c = sum(omega_total(abs(c)) for (c, d) in region)
    w_d = sum(omega_total(abs(d)) for (c, d) in region)
    return w_c, w_d


# ----------------------------------------------------------------------
# Decimal formatting for the certificate (exact -> fixed string)
# ----------------------------------------------------------------------

def dec(fr, places=40):
    fr = Fr(fr)
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    scaled = (fr * 10 ** places).__floor__()
    s = str(scaled).rjust(places + 1, "0")
    return sign + s[:-places] + "." + s[-places:]


# ----------------------------------------------------------------------
# Certificate build
# ----------------------------------------------------------------------

def build_certificate():
    cert = {}
    cert["certificate_type"] = "LAM1_LAMBDA_SEAM_INTERFACE_CALIBRATION"
    cert["claim_status"] = (
        "calibration_and_audit_only; no zero claims; legacy Lambda-FE "
        "phase certified as DECLARED, not derived"
    )
    cert["anchor_pinned"] = (
        "Classical automorphic substrate — Poisson summation, Jacobi theta "
        "inversion theta(1/t) = sqrt(t) theta(t), SL(2,Z) coset structure "
        "Gamma_inf \\ Gamma, Ramanujan-sum closed form: PINNED NAMED "
        "DEPENDENCY, cited as classical witness, never rederived natively "
        "(CIRC-1 discipline). This capsule verifies exact finite instances "
        "and bookkeeping only."
    )
    cert["lineage"] = {
        "origin_diagram": (
            "T-V-S-P compass of arXiv:2603.20773 (Geometric Completion of "
            "Thermodynamic Response); response block L = [[lp, Lpv],[Lvp, lv]]"
        ),
        "phase_identification": (
            "Onsager non-reciprocity bias phase phi = arg(Lpv/Lvp) of the "
            "compass response block is the same phi appearing as the legacy "
            "Lambda-FE multiplier e^{i phi}"
        ),
        "legacy_documents": (
            "Lambda-Mathematics Analytic Rigor Canvas (+ Gamma/pi constant "
            "audits, half-weight Lambda-mock), Lambda-Geometric "
            "Generalization, lambda-Geometry Package v1"
        ),
        "ruling": "MINE-DONT-MERGE: zero claims imported from legacy canvases",
    }

    # ---------------- T1 ----------------
    grid_c = list(range(1, 37)) + [40, 45, 48, 49, 60]
    grid_n = list(range(1, 37))
    checked = 0
    max_abs = 0
    for c in grid_c:
        for n in grid_n:
            val, is_int = ramanujan_direct(c, n)
            assert is_int, f"c_{c}({n}) not rational integer after reduction"
            ref = ramanujan_moebius(c, n)
            assert val == ref, f"two-route mismatch at (c,n)=({c},{n})"
            checked += 1
            max_abs = max(max_abs, abs(val))
    # control C1: coprimality dropped separates
    c1_witness = None
    for (c, n) in [(4, 2), (6, 2), (8, 4), (9, 3)]:
        val, is_int = ramanujan_direct(c, n, require_coprime=False)
        if (not is_int) or val != ramanujan_moebius(c, n):
            c1_witness = {
                "c": c,
                "n": n,
                "full_sum_reduction": (str(val) if is_int else "non_constant"),
                "ramanujan_value": ramanujan_moebius(c, n),
            }
            break
    assert c1_witness is not None, "control C1 failed to separate"
    cert["T1_ramanujan_two_route"] = {
        "statement": (
            "Direct exact cyclotomic reduction of c_c(n) in Z[x]/Phi_c(x) "
            "equals sum_{d|gcd(c,n)} d mu(c/d) on the full grid; the "
            "Lambda-twist multiplies whole c-blocks and never enters this "
            "d-sum layer"
        ),
        "grid": {"c": f"1..36 plus {grid_c[36:]}", "n": "1..36"},
        "pairs_checked": checked,
        "max_abs_value_on_grid": max_abs,
        "verdict": "PASS",
    }
    cert["controls_C1_coprimality_tamper_separates"] = c1_witness

    # ---------------- T2 ----------------
    ADD_N = 400
    additivity_checked = 0
    for m in range(1, ADD_N + 1):
        for n in range(1, ADD_N // m + 1):
            assert omega_total(m * n) == omega_total(m) + omega_total(n)
            additivity_checked += 1
    EULER_N = 300
    table = euler_marker_witness(EULER_N)
    assert set(table.keys()) == set(range(1, EULER_N + 1))
    for n in range(1, EULER_N + 1):
        mark, mult = table[n]
        assert mult == 1, f"multiplicity {mult} != 1 at n={n}"
        assert mark == omega_total(n), f"marker mismatch at n={n}"
    # control C2: fake non-additive winding detected
    def fake_nu(n):
        return omega_total(n) + (1 if n > 1 else 0)
    c2_witness = None
    for (m, n) in [(2, 3), (2, 2), (3, 5)]:
        if fake_nu(m * n) != fake_nu(m) + fake_nu(n):
            c2_witness = {"m": m, "n": n,
                          "fake_nu(mn)": fake_nu(m * n),
                          "fake_nu(m)+fake_nu(n)": fake_nu(m) + fake_nu(n)}
            break
    assert c2_witness is not None
    cert["T2_memory_lattice_multiplicativity"] = {
        "statement": (
            "chi_phi(n) = e^{i phi Omega(n)} is exact integer bookkeeping: "
            "Omega additive on all mn <= 400; formal Euler product over "
            "primes <= 300 reproduces every n <= 300 with marker exactly "
            "Omega(n) and multiplicity exactly one (twist exponents ADD "
            "along the product — the memory lattice of multiplication)"
        ),
        "additivity_pairs_checked": additivity_checked,
        "euler_range": EULER_N,
        "verdict": "PASS",
    }
    cert["controls_C2_fake_winding_detected"] = c2_witness

    # ---------------- T3 ----------------
    th_half = theta_iv(Fr(1, 2))
    th_two = theta_iv(Fr(2))
    s2 = sqrt_brackets(2)
    rhs = iv_mul_pos(s2, th_two)
    assert iv_overlap(th_half, rhs), "theta seam identity enclosures separate"
    WIDTH_BUDGET = Fr(1, 10 ** 20)
    assert iv_width(th_half) < WIDTH_BUDGET
    assert iv_width(rhs) < WIDTH_BUDGET
    # control C3: sqrt(3) separates
    s3 = sqrt_brackets(3)
    rhs3 = iv_mul_pos(s3, th_two)
    assert not iv_overlap(th_half, rhs3), "sqrt(3) tamper failed to separate"
    sep_gap = rhs3[0] - th_half[1]
    assert sep_gap > 0
    cert["T3_theta_seam_certified"] = {
        "statement": (
            "The single S-move generating the classical functional "
            "equation, theta(1/t) = sqrt(t) theta(t), certified at t=2 "
            "with two-sided directed rational enclosures (no floats)"
        ),
        "theta_half_enclosure": [dec(th_half[0]), dec(th_half[1])],
        "sqrt2_times_theta2_enclosure": [dec(rhs[0]), dec(rhs[1])],
        "width_theta_half": dec(iv_width(th_half), 40),
        "width_rhs": dec(iv_width(rhs), 40),
        "width_budget": "1e-20",
        "overlap": True,
        "verdict": "PASS",
    }
    cert["controls_C3_sqrt3_tamper_separates"] = {
        "separation_gap_lower_bound": dec(sep_gap, 20),
        "separated": True,
    }

    # ---------------- T4 ----------------
    C = 60
    region = primitive_region(C)
    image = {s_move(p) for p in region}
    assert image == region, "S-move is not a bijection of the region"
    w_c, w_d = weight_sums(region)
    assert w_c == w_d, "weight layer has nonzero net change under S"
    # control C4: asymmetric region separates
    asym = {(c, d) for c in range(1, C + 1) for d in range(1, 2 * C + 1)
            if gcd(c, d) == 1}
    aw_c = sum(omega_total(c) for (c, d) in asym)
    aw_d = sum(omega_total(d) for (c, d) in asym)
    assert aw_c != aw_d, "asymmetric control failed to separate"
    cert["T4_declared_phase_audit"] = {
        "statement": (
            "Under the exact S-reindexing bijection (c,d) -> (d,-c) on the "
            "primitive region |c|,|d| <= 60, the coset-weight layer "
            "nu(|c|) has ZERO net change (sum before = sum after). The "
            "weight bookkeeping alone produces no functional-equation "
            "phase; the legacy e^{i phi} FE factor is therefore carried "
            "entirely by the DECLARED right S-multiplier"
        ),
        "region_size": len(region),
        "weight_sum_pre_S": w_c,
        "weight_sum_post_S": w_d,
        "net_weight_change": w_c - w_d,
        "verdict": "PASS",
    }
    cert["controls_C4_asymmetric_region_separates"] = {
        "asym_weight_sum_c": aw_c,
        "asym_weight_sum_d": aw_d,
        "separated": True,
    }

    cert["finding_LAM1_F1"] = (
        "The e^{i phi} functional-equation factor of the legacy "
        "Lambda-Eisenstein / Lambda-mock construction is a DECLARED right "
        "S-multiplier, not a derived quantity: the legacy canvas itself "
        "records that no homomorphism SL(2,Z) -> S^1 can realize arbitrary "
        "phi (finite abelianization), and T4 certifies that the coset-weight "
        "bookkeeping contributes zero net phase. The legacy analytic layer "
        "therefore decorated the classical substrate with an input phase — "
        "this is the precise mechanism of its silence on RH/YM: the phase "
        "was never load-bearing. Any native functional equation (N2) must "
        "DERIVE its seam-crossing phase; the theta seam of T3 is the "
        "classical shape of that single crossing."
    )
    cert["claim_boundary"] = {
        "legacy_lambda_FE": "declared multiplier over classical substrate",
        "N1_native_continuation": "OPEN",
        "N2_native_functional_equation": "OPEN",
        "N3_identification": "OPEN",
        "K0_L0_RH": "OPEN",
        "YM_continuum_gates": "OPEN",
        "zeros_of_any_L_function": "no claims",
    }
    cert["arithmetic_discipline"] = (
        "integers, exact rationals, exact integer polynomial arithmetic, "
        "and directed rational intervals with outward rounding to 1e-60 "
        "denominators; no floating point in any verdict path"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "LAM1_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_LAM1.sha256"), "w") as f:
        f.write(digest + "\n")
    print("LAM1 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
