"""LAM-3: The crossing, derived — seam to functional equation, certified.

LAM-1 certified the theta seam theta(1/t) = sqrt(t) theta(t) and showed
the legacy Lambda phase was declared. LAM-2 showed what a derived phase
is. LAM-3 executes the DERIVATION OF THE s <-> 1-s CROSSING itself, at
machine level, with every value produced by the memory grammar
(recognized body + Bernoulli correction window + declared bracketed
tail) and every verdict an exact-rational interval statement.

The chain (classical shape, executed):

  theta seam (LAM-1 T3)
    ==> split-flip identity   [T2 here]
          int_0^1 t^{s/2-1} omega(t) dt
            = 1/(s-1) - 1/s + int_1^infty t^{(1-s)/2-1} omega(t) dt
    ==> manifestly s <-> 1-s symmetric completed kernel
    ==> equals the Dirichlet/Euler-product side on Re s > 1  [T3, T4]
          pi^{-s/2} Gamma(s/2) zeta(s)
            = 1/(s-1) - 1/s
              + int_1^infty (t^{s/2-1} + t^{(1-s)/2-1}) omega(t) dt

certified at s = 3 and s = 5 — points where the NATIVE zeta_Sigma
exists (Euler-product side). The functional-equation phase for zeta is
the DERIVED +1 of the untwisted seam: nothing is declared anywhere in
this chain.

Blocks:

T1  MEMORY-GRAMMAR VALUES. gamma (Euler-Mascheroni), ln 2, ln pi,
    zeta(3), zeta(5) as two-sided exact-rational enclosures via
    body + Bernoulli correction window + bracketed tail
    (Euler-Maclaurin, remainder bounded by first omitted term).
    Self-consistency: two independent truncation depths overlap.

T2  SPLIT-FLIP IDENTITY AT s=3, CERTIFIED. Left route: int_0^1 side
    via exact alternating incomplete-gamma series. Right route:
    polar terms 1/(s-1) - 1/s plus int_1^infty t^{-2} omega via the
    E_1 route (exact convergent series, exact cancellation). Overlap
    within budget. Controls: dropping the seam-generated polar terms
    separates (gap > 1/10); tampering the seam weight
    (omega(1/t) -> t*omega(t) + (t-1)/2) separates (gap > 2/5).

T3  COMPLETED IDENTITY AT s=3: zeta(3)/(2 pi) equals
    1/6 + int_1^infty (t^{1/2} + t^{-2}) omega(t) dt. Two independent
    enclosure routes overlap. Control: exponent tamper
    t^{(1-s)/2-1} -> t^{-s/2-1} (forgetting the seam's sqrt-t weight)
    separates.

T4  SECOND POINT s=5: 3 zeta(5)/(4 pi^2) equals
    1/20 + int_1^infty (t^{3/2} + t^{-3}) omega(t) dt. Kills
    point-luck.

Findings:
  LAM3-F1: the crossing is derived end-to-end and machine-certified on
  the Re s > 1 side; the FE phase +1 is an output of the seam, never
  an input.
  LAM3-F2: the memory grammar is load-bearing — body-only zeta
  separates from the certified value (public echo of "without memory
  terms there is no critical strip").

Honest boundary: continuation into the critical strip via the
symmetric kernel is DEFINITIONAL here, not a native theorem. The
classical substrate (theta seam, Mellin, Euler-Maclaurin remainder
envelopes) is pinned. N1/N2/N3 remain OPEN; what this capsule fixes is
the exact remaining native obligation: reproduce the split-flip step
(T2) inside the cut grammar, where the theta seam is a native theorem
rather than a pinned anchor. T3 is T2 composed with exact termwise
Mellin evaluation — recorded per TAUT-1 (independence is at the level
of enclosure routes, and the tails on both sides share the
Euler-Maclaurin utility; seam content lives in the leading digits,
where every tamper control separates by O(0.1)).
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

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
# Interval helpers (exact rationals throughout)
# ----------------------------------------------------------------------

def iv(lo, hi):
    lo, hi = Fr(lo), Fr(hi)
    assert lo <= hi
    return (lo, hi)


def iv_add(a, b):
    return lam1.round_out(a[0] + b[0], a[1] + b[1])


def iv_sub(a, b):
    return lam1.round_out(a[0] - b[1], a[1] - b[0])


def iv_neg(a):
    return (-a[1], -a[0])


def iv_scale(a, c):
    c = Fr(c)
    return (c * a[0], c * a[1]) if c >= 0 else (c * a[1], c * a[0])


def iv_inv_pos(a):
    assert a[0] > 0
    return lam1.round_out(1 / a[1], 1 / a[0])


BERNOULLI = {2: Fr(1, 6), 4: Fr(-1, 30), 6: Fr(1, 42), 8: Fr(-1, 30),
             10: Fr(5, 66), 12: Fr(-691, 2730), 14: Fr(7, 6)}


def factorial(n):
    out = 1
    for k in range(2, n + 1):
        out *= k
    return out


def rising(m, j):
    out = 1
    for i in range(j):
        out *= (m + i)
    return out


# ----------------------------------------------------------------------
# T1 machinery: memory-grammar values (body + correction + tail)
# ----------------------------------------------------------------------

def ln_atanh_at(x):
    """Two-sided enclosure of ln x for exact rational x >= 1 with
    (x-1)/(x+1) <= 1/2, via 2*atanh; all terms positive, geometric
    tail bound."""
    x = Fr(x)
    assert x >= 1
    r = (x - 1) / (x + 1)
    assert r < Fr(1, 2)
    r = lam1.round_dir(r, down=False, prec_bits=200)  # up: sound upper r
    s = Fr(0)
    K = 90
    rp = r
    r2 = r * r
    for k in range(K + 1):
        s += rp / (2 * k + 1)
        rp = lam1.round_dir(rp * r2, down=False, prec_bits=400)
    tail = 2 * rp / ((2 * K + 3) * (1 - r2))
    # r rounded up => s and tail are upper-biased; lower bound uses the
    # plain partial sum at rounded-up r minus nothing is not sound for
    # the lower side, so recompute with r rounded down.
    r_dn = lam1.round_dir((x - 1) / (x + 1), down=True, prec_bits=200)
    s_dn = Fr(0)
    rp = r_dn
    r2d = r_dn * r_dn
    for k in range(K + 1):
        s_dn += rp / (2 * k + 1)
        rp = lam1.round_dir(rp * r2d, down=True, prec_bits=400)
    return lam1.round_out(2 * s_dn, 2 * s + tail)


_LN_CACHE = {}


def ln2_iv():
    if "ln2" not in _LN_CACHE:
        _LN_CACHE["ln2"] = ln_atanh_at(Fr(2))
    return _LN_CACHE["ln2"]


def ln_pos(x):
    """Enclosure of ln x for exact rational x > 0, via range reduction
    x = 2^m * y with y in [1, 2)."""
    x = Fr(x)
    assert x > 0
    m = 0
    y = x
    while y >= 2:
        y /= 2
        m += 1
    while y < 1:
        y *= 2
        m -= 1
    base = ln_atanh_at(y)
    l2 = ln2_iv()
    lo = base[0] + m * (l2[0] if m >= 0 else l2[1])
    hi = base[1] + m * (l2[1] if m >= 0 else l2[0])
    return lam1.round_out(lo, hi)


def ln_pi_iv():
    if "lnpi" not in _LN_CACHE:
        pi_lo, pi_hi = lam1.pi_brackets()
        _LN_CACHE["lnpi"] = lam1.round_out(
            ln_pos(pi_lo)[0], ln_pos(pi_hi)[1])
    return _LN_CACHE["lnpi"]


def harmonic(N):
    s = Fr(0)
    for n in range(1, N + 1):
        s += Fr(1, n)
    return s


def gamma_euler_iv(N=128, K=6, bern=None):
    """Euler-Mascheroni gamma via Euler-Maclaurin at depth N:
    gamma = H_N - ln N - 1/(2N) + sum_{k<=K} B_2k/(2k) N^{-2k},
    remainder bounded by first omitted term (enveloping)."""
    bern = bern or BERNOULLI
    base = harmonic(N) - Fr(1, 2 * N)
    for k in range(1, K + 1):
        base += bern[2 * k] / (2 * k) / Fr(N) ** (2 * k)
    bound = abs(bern[2 * K + 2] / (2 * K + 2) / Fr(N) ** (2 * K + 2))
    lnN = ln_pos(N)
    return lam1.round_out(base - lnN[1] - bound, base - lnN[0] + bound)


def zeta_tail_iv(m, M, K=6):
    """Enclosure of sum_{n>=M} n^{-m}, integer m >= 2. Internally the
    Euler-Maclaurin start is pushed to max(M, 64) with an exact body in
    between, so the bracket is ~1e-28 regardless of M."""
    M0 = max(M, 64)
    body = Fr(0)
    for n in range(M, M0):
        body += Fr(1, n ** m)
    tail = Fr(1, (m - 1) * M0 ** (m - 1)) + Fr(1, 2 * M0 ** m)
    for k in range(1, K + 1):
        tail += (BERNOULLI[2 * k] / factorial(2 * k)
                 * rising(m, 2 * k - 1) / Fr(M0) ** (m + 2 * k - 1))
    bound = abs(BERNOULLI[2 * K + 2] / factorial(2 * K + 2)
                * rising(m, 2 * K + 1) / Fr(M0) ** (m + 2 * K + 1))
    return lam1.round_out(body + tail - bound, body + tail + bound)


def zeta_iv(m, N=64):
    """zeta(m) = recognized body + Bernoulli correction window +
    declared bracketed tail (memory grammar)."""
    body = Fr(0)
    for n in range(1, N):
        body += Fr(1, n ** m)
    t = zeta_tail_iv(m, N)
    return lam1.round_out(body + t[0], body + t[1])


# ----------------------------------------------------------------------
# Integral pieces on [0,1] and [1,inf) (exact series, exact cancellation)
# ----------------------------------------------------------------------

def _round_c(c_iv_):
    lo = lam1.round_dir(c_iv_[0], down=True, prec_bits=100)
    hi = lam1.round_dir(c_iv_[1], down=False, prec_bits=100)
    return lo, hi


def J_at(c, two_a_plus_1):
    """int_0^1 t^a e^{-ct} dt at exact rational c > 0, a = (two_a_plus_1-1)/2
    (half-integer): sum_k (-1)^k c^k / (k! (k+a+1)); alternating bracket."""
    c = Fr(c)
    K = 4 * (int(c) + 1) + 80
    s = Fr(0)
    term = Fr(2, two_a_plus_1 + 1)  # k=0: 1/(a+1) = 2/(2a+2)
    prev = None
    for k in range(K + 1):
        s += term if k % 2 == 0 else -term
        prev = term
        term = term * c / (k + 1) * Fr(2 * k + two_a_plus_1 + 1,
                                       2 * k + two_a_plus_1 + 3)
    assert term < prev, "alternating tail not yet decreasing"
    return lam1.round_out(s - term, s + term)


def J_iv(c_iv_, two_a_plus_1):
    lo, hi = _round_c(c_iv_)         # J decreasing in c
    a = J_at(hi, two_a_plus_1)
    b = J_at(lo, two_a_plus_1)
    return lam1.round_out(a[0], b[1])


def E1_at(c):
    """E_1(c) = -gamma - ln c + sum_{k>=1} (-1)^{k+1} c^k/(k k!),
    exact convergent series (cancellation is exact in rationals)."""
    c = Fr(c)
    K = 4 * (int(c) + 1) + 80
    s = Fr(0)
    term = c  # k=1: c/(1*1!)
    prev = None
    for k in range(1, K + 1):
        s += term if k % 2 == 1 else -term
        prev = term
        term = term * c * k / ((k + 1) * (k + 1))
    assert term < prev
    g = gamma_euler_iv()
    lnc = ln_pos(c)
    return lam1.round_out(-g[1] - lnc[1] + s - term,
                          -g[0] - lnc[0] + s + term)


def E1_iv(c_iv_):
    lo, hi = _round_c(c_iv_)         # E1 decreasing in c
    a = E1_at(hi)
    b = E1_at(lo)
    return lam1.round_out(a[0], b[1])


def I_m2_iv(c_iv_):
    """int_1^inf t^{-2} e^{-ct} dt = e^{-c} - c E_1(c)."""
    lo, hi = _round_c(c_iv_)
    e = lam1.exp_neg_iv(lo, hi)
    e1 = E1_iv((lo, hi))
    ce1 = lam1.round_out(lo * e1[0], hi * e1[1])
    out = iv_sub(e, ce1)
    assert out[0] > 0
    return out


def I_m3_iv(c_iv_):
    """int_1^inf t^{-3} e^{-ct} dt = e^{-c}/2 - (c/2) I_{-2}(c)."""
    lo, hi = _round_c(c_iv_)
    e = lam1.exp_neg_iv(lo, hi)
    im2 = I_m2_iv((lo, hi))
    half_c_im2 = lam1.round_out(lo * im2[0] / 2, hi * im2[1] / 2)
    out = iv_sub(iv_scale(e, Fr(1, 2)), half_c_im2)
    assert out[0] > 0
    return out


# ----------------------------------------------------------------------
# Assembly at s = 3 and s = 5
# ----------------------------------------------------------------------

def build_certificate():
    pi_iv = lam1.pi_brackets()
    inv_2pi = iv_inv_pos(iv_scale(pi_iv, 2))
    pi_sq = lam1.iv_mul_pos(pi_iv, pi_iv)

    def c_of(n):
        return (pi_iv[0] * n * n, pi_iv[1] * n * n)

    NMAX = 4

    # tail over n >= 5 of int_1^inf (t^{3/2}+t^{-3}) e^{-c t} dt
    # bounded by 3 e^{-c}/c per n (c >= 25 pi); summed with margin
    tail_hi = Fr(0)
    for n in range(5, 10):
        clo = lam1.round_dir(pi_iv[0] * n * n, down=True, prec_bits=100)
        tail_hi += 3 * lam1.exp_neg_brackets(clo)[1] / clo
    tail_hi *= 2  # margin absorbing the (geometric) remainder n >= 10
    tail_int = (Fr(0), tail_hi)

    # per-n pieces
    J12 = {n: J_iv(c_of(n), 2) for n in range(1, NMAX + 1)}   # a = 1/2
    J32 = {n: J_iv(c_of(n), 4) for n in range(1, NMAX + 1)}   # a = 3/2
    IM2 = {n: I_m2_iv(c_of(n)) for n in range(1, NMAX + 1)}
    IM3 = {n: I_m3_iv(c_of(n)) for n in range(1, NMAX + 1)}
    E1V = {n: E1_iv(c_of(n)) for n in range(1, NMAX + 1)}

    # Gamma(3/2)(pi n^2)^{-3/2} = n^{-3}/(2 pi); Gamma(5/2)(pi n^2)^{-5/2}
    # = 3 n^{-5}/(4 pi^2)  (half-integer Gamma cancels sqrt(pi) exactly)
    def mellin_full_12(n):
        return iv_scale(inv_2pi, Fr(1, n ** 3))

    def mellin_full_32(n):
        return iv_scale(iv_inv_pos(iv_scale(pi_sq, 4)), Fr(3, n ** 5))

    IP12 = {n: iv_sub(mellin_full_12(n), J12[n]) for n in range(1, NMAX + 1)}
    IP32 = {n: iv_sub(mellin_full_32(n), J32[n]) for n in range(1, NMAX + 1)}
    for n in range(1, NMAX + 1):
        assert IP12[n][0] > 0 and IP32[n][0] > 0

    cert = {}
    cert["certificate_type"] = "LAM3_CROSSING_DERIVATION_CERTIFIED"
    cert["claim_status"] = (
        "s<->1-s crossing derivation executed and certified at s=3 and "
        "s=5 on the Euler-product side; memory grammar load-bearing; "
        "no zero claims; continuation into the strip is definitional "
        "here, not a native theorem"
    )
    cert["anchor_pinned"] = (
        "Classical substrate — Jacobi theta seam (certified instance in "
        "LAM-1 T3), termwise Mellin evaluation, Euler-Maclaurin "
        "enveloping remainders for completely monotone integrands: "
        "PINNED NAMED DEPENDENCIES (CIRC-1), cited and never rederived "
        "natively. This capsule certifies exact finite instances."
    )

    # ---------------- T1 ----------------
    g128 = gamma_euler_iv(128)
    g64 = gamma_euler_iv(64)
    assert lam1.iv_overlap(g128, g64)
    z3 = zeta_iv(3, 64)
    z3b = zeta_iv(3, 128)
    assert lam1.iv_overlap(z3, z3b)
    z5 = zeta_iv(5, 64)
    z5b = zeta_iv(5, 128)
    assert lam1.iv_overlap(z5, z5b)
    ln2 = ln2_iv()
    lnpi = ln_pi_iv()
    for v in (g128, z3, z5, ln2, lnpi):
        assert lam1.iv_width(v) < Fr(1, 10 ** 24)
    # C3 memory control: body-only zeta(3) separates
    body_only = sum(Fr(1, n ** 3) for n in range(1, 64))
    assert body_only < z3[0] and z3[0] - body_only > Fr(1, 10 ** 5)
    # C4 gamma self-consistency + Bernoulli tamper
    tampered = dict(BERNOULLI)
    tampered[4] = Fr(1, 30)
    g64_bad = gamma_euler_iv(64, bern=tampered)
    assert not lam1.iv_overlap(g64_bad, g128)
    cert["T1_memory_grammar_values"] = {
        "statement": (
            "gamma, ln 2, ln pi, zeta(3), zeta(5) as two-sided "
            "exact-rational enclosures via recognized body + Bernoulli "
            "correction window + declared bracketed tail; two "
            "truncation depths overlap for gamma, zeta(3), zeta(5)"
        ),
        "gamma": [lam1.dec(g128[0]), lam1.dec(g128[1])],
        "zeta3": [lam1.dec(z3[0]), lam1.dec(z3[1])],
        "zeta5": [lam1.dec(z5[0]), lam1.dec(z5[1])],
        "ln2": [lam1.dec(ln2[0]), lam1.dec(ln2[1])],
        "ln_pi": [lam1.dec(lnpi[0]), lam1.dec(lnpi[1])],
        "width_budget": "1e-24",
        "verdict": "PASS",
    }
    cert["controls_C3_body_only_zeta_separates"] = {
        "body_only_n_lt_64": lam1.dec(body_only, 20),
        "certified_lower": lam1.dec(z3[0], 20),
        "gap_exceeds": "1e-5",
        "separated": True,
    }
    cert["controls_C4_bernoulli_tamper_separates"] = {
        "tamper": "B_4 sign flip in the correction window",
        "separated": True,
    }

    # ---------------- T2 ----------------
    # LHS: int_0^1 t^{1/2} omega dt = sum_n J(pi n^2; a=1/2)
    #      n<=4 exact series; n>=5 via n^{-3}/(2 pi) minus tiny I_p bound
    lhs = (Fr(0), Fr(0))
    for n in range(1, NMAX + 1):
        lhs = iv_add(lhs, J12[n])
    ztail5 = zeta_tail_iv(3, NMAX + 1)
    lhs = iv_add(lhs, lam1.iv_mul_pos(inv_2pi, ztail5))
    ip_tail_bound = Fr(0)
    for n in range(5, 9):
        clo = lam1.round_dir(pi_iv[0] * n * n, down=True, prec_bits=100)
        ip_tail_bound += 2 * lam1.exp_neg_brackets(clo)[1] / clo
    lhs = iv_sub(lhs, (Fr(0), 2 * ip_tail_bound))
    # RHS: 1/(s-1) - 1/s + sum_n I_{-2}(pi n^2)  (E_1 route)
    rhs_int = (Fr(0), Fr(0))
    for n in range(1, NMAX + 1):
        rhs_int = iv_add(rhs_int, IM2[n])
    rhs_int = iv_add(rhs_int, tail_int)
    rhs = iv_add((Fr(1, 6), Fr(1, 6)), rhs_int)
    budget = Fr(1, 10 ** 20)
    assert lam1.iv_overlap(lhs, rhs), "split-flip identity separates"
    assert lam1.iv_width(lhs) < budget and lam1.iv_width(rhs) < budget
    # C1 polar drop
    gap_c1 = lhs[0] - rhs_int[1]
    assert gap_c1 > Fr(1, 10)
    # C2 wrong seam weight: omega(1/t) -> t omega(t) + (t-1)/2 gives
    # 2/3 + sum_n I_{-3/2}; and I_{-2} <= I_{-3/2} <= E_1 on [1,inf)
    wrong_lo = Fr(2, 3)
    for n in range(1, NMAX + 1):
        wrong_lo += IM2[n][0]
    gap_c2 = wrong_lo - lhs[1]
    assert gap_c2 > Fr(2, 5)
    cert["T2_split_flip_identity_s3"] = {
        "statement": (
            "int_0^1 t^{1/2} omega(t) dt = 1/2 - 1/3 + int_1^inf "
            "t^{-2} omega(t) dt, certified by two independent enclosure "
            "routes (incomplete-gamma series vs E_1 route with exact "
            "rational cancellation); this identity holds precisely "
            "because of the theta seam — it is the crossing step"
        ),
        "lhs_enclosure": [lam1.dec(lhs[0]), lam1.dec(lhs[1])],
        "rhs_enclosure": [lam1.dec(rhs[0]), lam1.dec(rhs[1])],
        "width_budget": "1e-20",
        "overlap": True,
        "verdict": "PASS",
    }
    cert["controls_C1_polar_drop_separates"] = {
        "gap_lower_bound": lam1.dec(gap_c1, 20), "gap_exceeds": "1/10",
        "separated": True,
    }
    cert["controls_C2_seam_weight_tamper_separates"] = {
        "tamper": "omega(1/t) -> t*omega(t) + (t-1)/2 (wrong weight)",
        "gap_lower_bound": lam1.dec(gap_c2, 20), "gap_exceeds": "2/5",
        "separated": True,
    }

    # ---------------- T3 ----------------
    lhsB = lam1.iv_mul_pos(z3, inv_2pi)
    rhsB_int = (Fr(0), Fr(0))
    for n in range(1, NMAX + 1):
        rhsB_int = iv_add(rhsB_int, iv_add(IM2[n], IP12[n]))
    rhsB_int = iv_add(rhsB_int, tail_int)
    rhsB = iv_add((Fr(1, 6), Fr(1, 6)), rhsB_int)
    assert lam1.iv_overlap(lhsB, rhsB), "completed identity s=3 separates"
    assert lam1.iv_width(lhsB) < budget and lam1.iv_width(rhsB) < budget
    # C5 exponent tamper: t^{(1-s)/2-1} -> t^{-s/2-1}; then the
    # t^{1/2}-piece becomes I_{-5/2} with I_{-3} <= I_{-5/2} <= I_{-2}
    wrongB_hi = Fr(1, 6) + tail_hi
    for n in range(1, NMAX + 1):
        wrongB_hi += IM2[n][1] + IM2[n][1]
    gap_c5 = lhsB[0] - wrongB_hi
    assert gap_c5 > Fr(1, 200)
    cert["T3_completed_identity_s3"] = {
        "statement": (
            "zeta(3)/(2 pi) = 1/6 + int_1^inf (t^{1/2} + t^{-2}) "
            "omega(t) dt — the manifestly s<->1-s symmetric completed "
            "kernel equals the Dirichlet/Euler-product side at s=3; "
            "LHS via memory-grammar zeta(3), RHS via seam integrals"
        ),
        "lhs_enclosure": [lam1.dec(lhsB[0]), lam1.dec(lhsB[1])],
        "rhs_enclosure": [lam1.dec(rhsB[0]), lam1.dec(rhsB[1])],
        "width_budget": "1e-20",
        "overlap": True,
        "taut1_note": (
            "T3 = T2 composed with exact termwise Mellin evaluation; "
            "verified as independently assembled enclosures; tails on "
            "both sides share the Euler-Maclaurin utility"
        ),
        "verdict": "PASS",
    }
    cert["controls_C5_exponent_tamper_separates"] = {
        "tamper": "t^{(1-s)/2-1} -> t^{-s/2-1} (seam sqrt-t weight "
                  "forgotten in the reflected kernel)",
        "gap_lower_bound": lam1.dec(gap_c5, 20), "gap_exceeds": "1/200",
        "separated": True,
    }

    # ---------------- T4 ----------------
    lhsC = lam1.iv_mul_pos(z5, iv_inv_pos(iv_scale(pi_sq, Fr(4, 3))))
    rhsC_int = (Fr(0), Fr(0))
    for n in range(1, NMAX + 1):
        rhsC_int = iv_add(rhsC_int, iv_add(IM3[n], IP32[n]))
    rhsC_int = iv_add(rhsC_int, tail_int)
    rhsC = iv_add((Fr(1, 20), Fr(1, 20)), rhsC_int)
    assert lam1.iv_overlap(lhsC, rhsC), "completed identity s=5 separates"
    assert lam1.iv_width(lhsC) < budget and lam1.iv_width(rhsC) < budget
    cert["T4_completed_identity_s5"] = {
        "statement": (
            "3 zeta(5)/(4 pi^2) = 1/20 + int_1^inf (t^{3/2} + t^{-3}) "
            "omega(t) dt — second point, same derived crossing; kills "
            "point-luck"
        ),
        "lhs_enclosure": [lam1.dec(lhsC[0]), lam1.dec(lhsC[1])],
        "rhs_enclosure": [lam1.dec(rhsC[0]), lam1.dec(rhsC[1])],
        "width_budget": "1e-20",
        "overlap": True,
        "verdict": "PASS",
    }

    cert["finding_LAM3_F1"] = (
        "The s<->1-s crossing is DERIVED end-to-end at machine level on "
        "the Euler-product side: theta seam (LAM-1 T3, certified) => "
        "split-flip identity (T2, certified; seam tamper separates by "
        ">2/5, polar drop by >1/10) => manifestly symmetric completed "
        "kernel = Dirichlet side (T3 at s=3, T4 at s=5, certified "
        "two-route). The functional-equation phase for zeta is the "
        "derived +1 of the untwisted seam — an OUTPUT of the crossing, "
        "never an input. Contrast chain complete: LAM-1 declared phase "
        "(never load-bearing) / LAM-2 derived phase (load-bearing, "
        "finite seam arithmetic) / LAM-3 the derivation itself executed."
    )
    cert["finding_LAM3_F2"] = (
        "The memory grammar is load-bearing: body-only zeta(3) "
        "separates from the certified value by more than 1e-5 at "
        "enclosure widths below 1e-24; every LHS in T3/T4 is a "
        "body+correction+tail object. Public echo of the framework "
        "fact that without memory terms there is no critical strip."
    )
    cert["claim_boundary"] = {
        "certified_content": (
            "exact finite instances at s=3, s=5 on Re s > 1, where the "
            "native zeta_Sigma exists as Euler product"
        ),
        "continuation_into_strip": (
            "definitional via the symmetric kernel here — NOT a native "
            "theorem"
        ),
        "remaining_native_obligation_N2": (
            "reproduce the split-flip step (T2) inside the cut grammar, "
            "where the theta seam is a native theorem rather than a "
            "pinned classical anchor"
        ),
        "N1_native_continuation": "OPEN",
        "N2_native_functional_equation": "OPEN",
        "N3_identification": "OPEN",
        "K0_L0_RH": "OPEN",
        "YM_continuum_gates": "OPEN",
        "zeros_of_any_L_function": "no claims",
    }
    cert["arithmetic_discipline"] = (
        "exact rationals throughout; alternating and enveloping series "
        "with remainder <= first omitted term; exact cancellation in "
        "the E_1 route (no floating point anywhere, hence no "
        "catastrophic cancellation); directed relative rounding to "
        "100-300 bit significands; every verdict an interval statement"
    )
    return cert


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out_path = os.path.join(HERE, "LAM3_RESULT.json")
    with open(out_path, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(HERE, "EXPECTED_LAM3.sha256"), "w") as f:
        f.write(digest + "\n")
    print("LAM3 certificate written:", out_path)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
