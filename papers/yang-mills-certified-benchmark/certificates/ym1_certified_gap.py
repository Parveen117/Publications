"""YM-1: Certified reduced cutoff gap for the SU(2) one-holonomy benchmark.

Claim boundary (declared, fail-closed):
  - This certifies the REDUCED FINITE-CUTOFF gap
        Delta_red(a, beta) = -(1/a) * log( I_2(beta) / I_1(beta) )
    of the normalized SU(2) Wilson convolution operator on class functions
    (MP adapters/yang_mills/cutoff_gap_benchmark.tex), as an exact two-sided
    rational enclosure. No floats participate in any verdict.
  - It is NOT a full-lattice gap, NOT a continuum gap, NOT the Clay predicate.
  - Source spectrum: lambda_j(beta) = I_{2j+1}(beta) / I_1(beta).

Method:
  - I_nu(beta) for integer nu >= 0 via the everywhere-positive series
        I_nu(x) = sum_k (x/2)^(nu+2k) / (k! (k+nu)!)
    with exact Fraction partial sums and a geometric-majorant tail bracket:
    term ratio t_{k+1}/t_k = (x/2)^2 / ((k+1)(k+nu+1)) <= r < 1 for k >= K,
    hence tail <= t_{K+1} / (1 - r).
  - log via atanh series: log(x) = 2*atanh((x-1)/(x+1)) for x > 0, with
    exact alternating-free positive series and geometric tail bound.
  - Directed rational rounding throughout; enclosure arithmetic only.

Controls:
  - C1 monotonicity: lambda_1(beta) strictly inside (0,1) => gap > 0.
  - C2 tamper: replacing I_2 by I_3 must strictly change the bracket.
"""

from fractions import Fraction as F
import hashlib
import json
import os
import sys

sys.set_int_max_str_digits(200000)  # exact rationals get large; verdicts stay exact

A_LATTICE = F(1)          # a = 1 (benchmark point)
BETA = F(2)               # beta = 2 (benchmark point)
TERMS = 60                # series depth (tail certified regardless)
LOG_TERMS = 80            # atanh series depth


# ---------------------------------------------------------------- intervals
class Iv:
    """Closed rational interval [lo, hi] with directed arithmetic."""

    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        lo = F(lo)
        hi = lo if hi is None else F(hi)
        if lo > hi:
            raise ValueError("empty interval")
        self.lo, self.hi = lo, hi

    def __add__(self, o):
        o = _iv(o)
        return Iv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = _iv(o)
        return Iv(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = _iv(o)
        c = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return Iv(min(c), max(c))

    def __truediv__(self, o):
        o = _iv(o)
        if o.lo <= 0 <= o.hi:
            raise ZeroDivisionError("interval contains zero")
        c = [self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi]
        return Iv(min(c), max(c))

    def __neg__(self):
        return Iv(-self.hi, -self.lo)

    def width(self):
        return self.hi - self.lo

    def strictly_inside(self, lo, hi):
        return F(lo) < self.lo and self.hi < F(hi)

    def separated_from(self, o):
        o = _iv(o)
        return self.hi < o.lo or o.hi < self.lo

    def __repr__(self):
        return f"Iv[{self.lo}, {self.hi}]"


def _iv(x):
    return x if isinstance(x, Iv) else Iv(x)


# ------------------------------------------------------- certified Bessel I
def bessel_I(nu: int, x: F, terms: int) -> Iv:
    """Exact enclosure of I_nu(x), integer nu >= 0, rational x > 0."""
    assert nu >= 0 and x > 0
    half = x / 2
    hs = half * half
    # t_0 = (x/2)^nu / nu!
    t = half ** nu
    fact = 1
    for i in range(2, nu + 1):
        fact *= i
    t = F(t, fact)
    s = t
    for k in range(1, terms + 1):
        t = t * hs / (k * (k + nu))
        s += t
    # tail bound: for k >= terms, ratio <= r := hs/((terms+1)(terms+nu+1)) < 1
    r = hs / ((terms + 1) * (terms + nu + 1))
    if r >= 1:
        raise ValueError("increase terms: tail ratio not < 1")
    t_next = t * hs / ((terms + 1) * (terms + nu + 1))
    tail = t_next / (1 - r)
    return Iv(s, s + tail)


# ---------------------------------------------------------- certified log
def log_iv(x: Iv, terms: int) -> Iv:
    """Enclosure of log over a positive rational interval via atanh series."""
    if x.lo <= 0:
        raise ValueError("log needs positive interval")
    return Iv(_log_lo(x.lo, terms), _log_hi(x.hi, terms))


def _log_point(x: F, terms: int):
    """Two-sided bracket of log(x) at a rational point, x > 0."""
    inv = False
    if x == 1:
        return F(0), F(0)
    if x < 1:
        x, inv = 1 / x, True
    u = (x - 1) / (x + 1)          # 0 < u < 1
    u2 = u * u
    term = u
    s = term
    for k in range(1, terms + 1):
        term = term * u2
        s += term / (2 * k + 1)
    # tail: sum_{k>terms} u^{2k+1}/(2k+1) <= u^{2*terms+3}/((2*terms+3)(1-u2))
    tail = (term * u2) / ((2 * terms + 3) * (1 - u2))
    lo, hi = 2 * s, 2 * (s + tail)
    if inv:
        lo, hi = -hi, -lo
    return lo, hi


def _log_lo(x, terms):
    return _log_point(F(x), terms)[0]


def _log_hi(x, terms):
    return _log_point(F(x), terms)[1]


# ------------------------------------------------------------------ theorem
def certified_reduced_gap(a: F, beta: F, terms=TERMS, log_terms=LOG_TERMS):
    I1 = bessel_I(1, beta, terms)
    I2 = bessel_I(2, beta, terms)
    lam1 = I2 / I1                              # lambda_1 = I_2/I_1
    gap = -(log_iv(lam1, log_terms)) / Iv(a)    # Delta = -(1/a) log lambda_1
    return I1, I2, lam1, gap


def run():
    I1, I2, lam1, gap = certified_reduced_gap(A_LATTICE, BETA)

    # C1: spectral sanity — lambda_1 strictly inside (0,1) => gap strictly > 0
    c1 = lam1.strictly_inside(0, 1) and gap.lo > 0

    # C2: tamper control — using I_3 (i.e. lambda_{3/2}-slot abuse) separates
    I3 = bessel_I(3, BETA, TERMS)
    lam_tamper = I3 / I1
    gap_tamper = -(log_iv(lam_tamper, LOG_TERMS)) / Iv(A_LATTICE)
    c2 = gap_tamper.separated_from(gap)

    # C3: width budget — enclosure tighter than 1e-30
    c3 = gap.width() < F(1, 10 ** 30)

    ok = c1 and c2 and c3
    cert = {
        "certificate_type": "YM1_CERTIFIED_REDUCED_CUTOFF_GAP",
        "claim_status": "reduced_finite_cutoff_only",
        "claim_boundary": {
            "certified": "Delta_red(a=1, beta=2) two-sided exact-rational enclosure",
            "not_certified": [
                "full-lattice gap", "continuum gap", "Clay mass-gap predicate",
            ],
            "source": "MP adapters/yang_mills/cutoff_gap_benchmark.tex",
        },
        "parameters": {"a": str(A_LATTICE), "beta": str(BETA),
                       "bessel_terms": TERMS, "log_terms": LOG_TERMS},
        "enclosures": {
            "I1_lo": str(I1.lo), "I1_hi": str(I1.hi),
            "I2_lo": str(I2.lo), "I2_hi": str(I2.hi),
            "lambda1_lo": str(lam1.lo), "lambda1_hi": str(lam1.hi),
            "gap_lo": str(gap.lo), "gap_hi": str(gap.hi),
            "gap_lo_decimal_40": _dec(gap.lo, 40),
            "gap_hi_decimal_40": _dec(gap.hi, 40),
            "width_lt_1e30": bool(c3),
        },
        "controls": {
            "C1_lambda1_in_unit_interval_gap_positive": bool(c1),
            "C2_I3_tamper_separates": bool(c2),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


def _dec(q: F, places: int) -> str:
    sign = "-" if q < 0 else ""
    q = abs(q)
    ip = q.numerator // q.denominator
    rem = q.numerator - ip * q.denominator
    digits = []
    for _ in range(places):
        rem *= 10
        d = rem // q.denominator
        digits.append(str(d))
        rem -= d * q.denominator
    return f"{sign}{ip}." + "".join(digits)


def canonical_sha(cert) -> str:
    blob = json.dumps(cert, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


if __name__ == "__main__":
    cert = run()
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "YM1_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(here, "EXPECTED_YM1.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"])
    print("gap in [", cert["enclosures"]["gap_lo_decimal_40"], ",",
          cert["enclosures"]["gap_hi_decimal_40"], "]")
    print("sha256:", sha)
