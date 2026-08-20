"""YM-13 INDEPENDENT WITNESSES — deliberately imports NOTHING from any YM
capsule. Only fractions/typing. Each function is an independently executed
obligation replacing a shared witness slot in the verification Gram; the
import isolation is machine-checked by ym13 via ast.

W1  Bessel ratio I_{nu+1}(x)/I_nu(x) by CONTINUED FRACTION with certified
    alternating truncation — a genuinely different algorithm from the
    power-series + geometric-tail route of the arithmetic engine.
    r_nu = 1 / ( 2(nu+1)/x + r_{nu+1} ),  r >= 0, and the CF map is
    DECREASING in the tail, so seeding the deepest tail with the two-sided
    a-priori bracket 0 <= r <= x/(2(nu+2)) and propagating both seeds up
    gives a certified two-sided enclosure.

W2  e^x two-sided bounds by the COMPOUND-INTEREST inequalities
    (1 + x/n)^n <= e^x <= (1 - x/n)^{-n}  for 0 <= x < n — pure rational
    powers, independent of the Taylor-series exp of the engine.

W3  Casimir recomputed from scratch: C(j) = j(j+1), twice-spin input.

W4  2x2 symmetric inertia by the DETERMINANT/TRACE sign rule (independent
    of the LDL elimination engine): for [[a,b],[b,c]],
    det > 0 -> both eigenvalues share the sign of a (2,0) or (0,2);
    det < 0 -> (1,1); det = 0 -> degenerate (refused).
"""

from fractions import Fraction as F


# ------------------------------------------------------------------- W1
def bessel_ratio_cf(nu: int, x: F, depth: int = 40):
    """Certified enclosure of I_{nu+1}(x)/I_nu(x) via continued fraction."""
    assert x > 0 and nu >= 0 and depth >= 1
    top = nu + depth
    lo_seed = F(0)
    hi_seed = x / (2 * (top + 2))
    vals = []
    for seed in (lo_seed, hi_seed):
        r = seed
        for k in range(top, nu - 1, -1):
            r = 1 / (2 * (k + 1) / x + r)
        vals.append(r)
    lo, hi = min(vals), max(vals)
    return lo, hi


# ------------------------------------------------------------------- W2
def exp_bounds_compound(x: F, n: int = 4096):
    """(1 + x/n)^n <= e^x <= (1 - x/n)^{-n} for 0 <= x < n; for x < 0 use
    reciprocals. Pure rational arithmetic."""
    if x == 0:
        return F(1), F(1)
    if x < 0:
        lo_p, hi_p = exp_bounds_compound(-x, n)
        return 1 / hi_p, 1 / lo_p
    assert x < n
    lo = (1 + x / n) ** n
    hi = (1 - x / n) ** (-n)
    return lo, hi


# ------------------------------------------------------------------- W3
def casimir_independent(twice_spin: int) -> F:
    j = F(twice_spin, 2)
    return j * (j + 1)


# ------------------------------------------------------------------- W4
def inertia_2x2_dettrace(a: F, b: F, c: F):
    """Inertia of [[a,b],[b,c]] by determinant/trace signs; None if
    singular (refused, matching the engine's fail-closed convention)."""
    det = a * c - b * b
    tr = a + c
    if det == 0:
        return None
    if det < 0:
        return (1, 1)
    return (2, 0) if tr > 0 else (0, 2)
