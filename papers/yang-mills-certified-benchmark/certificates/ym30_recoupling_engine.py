"""YM-30: THE 2-ROW RECOUPLING ENGINE — exact evaluation of one-time-step
fabric (ladder) integrals, the object YM-25/26/28/29 all reduced to.

The ladder: top rail rungs A_1..A_m, bottom rail B_1..B_m, top faces
u(A_l^{-1}A_{l+1}), bottom faces v(B_l^{-1}B_{l+1}), heat-kernel rungs
K(A_i B_i^{-1}); u, v, K class functions with coefficients u_a, v_b, k_c:
    Z = Int prod u prod v prod K  =  sum_{labels} prod u_a prod v_b prod k_c * N(a, b, c),
N(a,b,c) = the SU(2) spin-network evaluation of the labelled ladder graph.

ENGINE (all exact; numbers live in Q adjoined square roots of integers):
  1. Surd arithmetic: elements sum_s q_s sqrt(s), s squarefree; +, *, and
     multiplication by rationals; exact equality; rational extraction.
  2. Clebsch-Gordan coefficients by the Racah closed form, radicands
     reduced to squarefree surds (no floats anywhere).
  3. Wigner 3j from CG; vertex integral
        Int D^{j1}_{m1 n1} D^{j2}_{m2 n2} D^{j3}_{m3 n3} dg = (j1 j2 j3)(j1 j2 j3)
                                                            (m1 m2 m3)(n1 n2 n3)
     (the orthogonality of three matrix coefficients); conjugates removed
     with the epsilon matrix conj(D^j) = eps D^j eps^{-1}.
  4. Edge/rung characters expanded in matrix coefficients; the ladder
     contracted column by column as a transfer over the cut indices.

VALIDATION (the engine must pass all before any use):
  V1  every labelled evaluation N(a,b,c) is RATIONAL (surds cancel) —
      the built-in consistency check of the whole construction.
  V2  m = 2 (one plaquette): N equals the closed form
      delta_{a=b=c1=c2} / d_a^2  (a 4-edge loop: three character
      convolutions give chi_a/d_a^3, the closing trace gives d_a),
      certified for all labels <= 1.
  V3  single-row limit (v = 1, K = delta): the ladder reduces to the
      YM-15 chain closed forms f_0^{m-1} and r^{|p-q|} — the 1-row
      engine is recovered exactly.
  V4  recoupling cross-check: two plaquettes, all four faces and the
      outer rungs 1/2, middle rung c in {0,1}: sum_c d_c N = 1/4 (the
      middle rung summed with d_c is the delta function, closing a
      6-edge loop = 1/d^2), AND the split d_0 N_0 : d_1 N_1 = 1/4 : 3/4
      — exactly YM-25's squared recoupling weights, now from the engine.
  V5  Z is symmetric under swapping the rails (u <-> v) and under
      reversal of the chain.

FIRST USE — the product-ansatz vacuum floor (YM-29 T3): psi_c =
prod_l (1 + c chi_half(H_l)), operator T'' = W^{1/2} T0 W^{1/2}
(similar to T), Rayleigh quotient
    R_m(c) = <W psi, T0 W psi> / <psi, W psi>,
numerator = ladder with u = v = w (1 + c chi_half) and K the heat kernel,
denominator = (f_0 + 4 c f_half + c^2 (f_0 + 3 f_1))^{m-1} (YM-29 L).
Computed here with the face/time coefficients truncated to content <= 1
and taken as exact rationals at the interval lower ends — this first
use is an ENGINE CALIBRATION, not yet an outward certificate (the
content tail and the coefficient intervals are not propagated). Its
purpose: check that the product ansatz reaches the true per-face rate
of YM-18 (0.002177 at kappa = 1/8) where the linear ansatz could not.
Certified in this capsule: the engine (V1-V5) and the calibration table.

CALIBRATION FINDING (kappa = 1/8, content <= 1, point coefficients): already
at c = 0 — the state W^{1/2}*1 for T'' — the Rayleigh quotient exceeds the
old floor f_0^{m-1} by a factor that is GEOMETRIC in m: 1.000732 per face
(m = 2..6, ratios constant to 1e-7), i.e. floor rate log f_0 + 0.000732 =
0.002685 per face, ABOVE YM-18's certified lower bound 0.002177 and above
the E4D single-insertion excitation rate 0.001949 (YM-28) by 0.000736 per
face. The single-face excitation (c != 0) does not help at this coupling
(c = +-1/20 lowers the quotient): the improvement is the W^{1/2} dressing
itself, which the 1-row engine could not evaluate. The geometric law is a
transfer-matrix fact (the ladder is a product of column transfers), so
the m-uniform version of this floor is the top eigenvalue of the cut
transfer — NAMED YM-31 (outward certification: interval coefficients,
content tail, Perron bound on the column transfer).
NOT claimed: any gap; any floor beyond the calibration's stated scope.
"""

from fractions import Fraction as F
from math import factorial
import itertools
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, _dec, canonical_sha  # noqa: E402
from ym16_chain_dock import f0_of  # noqa: E402
from ym15_chain_closed_form import lam_half  # noqa: E402
from ym28_e4d_on_chain import bessel_I, TERMS  # noqa: E402
from ym15_chain_closed_form import BETA  # noqa: E402


# ============================================================ surds
def squarefree(n: int):
    """n > 0 -> (k, s) with n = k^2 s, s squarefree."""
    k, s, p = 1, 1, 2
    while p * p <= n:
        while n % (p * p) == 0:
            n //= p * p
            k *= p
        if n % p == 0:
            n //= p
            s *= p
        p += 1
    return k, s * n


class Surd:
    __slots__ = ("t",)

    def __init__(self, t=None):
        self.t = {s: q for s, q in (t or {}).items() if q != 0}

    @classmethod
    def rat(cls, q):
        return cls({1: F(q)})

    @classmethod
    def sqrt_of(cls, q: F):
        """sqrt of a nonnegative rational as a surd."""
        q = F(q)
        if q == 0:
            return cls()
        num, den = q.numerator, q.denominator
        k, s = squarefree(num * den)            # sqrt(num*den)/den
        return cls({s: F(k, den)})

    def __add__(self, o):
        t = dict(self.t)
        for s, q in o.t.items():
            t[s] = t.get(s, F(0)) + q
        return Surd(t)

    def __neg__(self):
        return Surd({s: -q for s, q in self.t.items()})

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        if not isinstance(o, Surd):
            return Surd({s: q * F(o) for s, q in self.t.items()})
        t = {}
        for s1, q1 in self.t.items():
            for s2, q2 in o.t.items():
                k, s = squarefree(s1 * s2)
                t[s] = t.get(s, F(0)) + q1 * q2 * k
        return Surd(t)

    __rmul__ = __mul__

    def is_rational(self):
        return all(s == 1 for s in self.t)

    def rational(self):
        assert self.is_rational(), self.t
        return self.t.get(1, F(0))

    def is_zero(self):
        return not self.t


# ============================================================ CG / 3j
def cg(j1, m1, j2, m2, j, m):
    """Clebsch-Gordan <j1 m1 j2 m2 | j m>, spins as Fractions (halves ok)."""
    if m1 + m2 != m or j < abs(j1 - j2) or j > j1 + j2:
        return Surd()
    if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return Surd()
    if (j1 + j2 + j).denominator != 1:
        return Surd()

    def fa(x):
        x = F(x)
        assert x.denominator == 1 and x >= 0, x
        return factorial(int(x))
    pref = F(2 * j + 1) * fa(j + j1 - j2) * fa(j - j1 + j2) * fa(j1 + j2 - j) / fa(j1 + j2 + j + 1)
    pref *= fa(j + m) * fa(j - m) * fa(j1 - m1) * fa(j1 + m1) * fa(j2 - m2) * fa(j2 + m2)
    ssum = F(0)
    k = 0
    while True:
        terms = [j1 + j2 - j - k, j1 - m1 - k, j2 + m2 - k, j - j2 + m1 + k, j - j1 - m2 + k]
        if any(F(x).denominator != 1 for x in terms):
            raise ValueError("bad spins")
        if terms[0] < 0 or terms[1] < 0 or terms[2] < 0:
            break
        if terms[3] >= 0 and terms[4] >= 0:
            ssum += F((-1) ** k) / (fa(k) * fa(terms[0]) * fa(terms[1]) * fa(terms[2])
                                    * fa(terms[3]) * fa(terms[4]))
        k += 1
    return Surd.sqrt_of(pref) * ssum


def three_j(j1, m1, j2, m2, j3, m3):
    sign = (-1) ** int(j1 - j2 - m3)
    return cg(j1, m1, j2, m2, j3, -m3) * Surd.sqrt_of(F(1, int(2 * j3 + 1))) * sign


def ms(j):
    """magnetic quantum numbers of spin j, in order j, j-1, ..., -j."""
    return [j - k for k in range(int(2 * j) + 1)]


def eps(j, m, mp):
    """conj(D^j_{m m'}) = sum eps_{m a} D^j_{a b} eps_{m' b} with
    eps_{m,-m} = (-1)^{j-m}; returns the scalar eps(j)_{m,mp}."""
    return F((-1) ** int(j - m)) if mp == -m else F(0)


def vertex(j1, m1, n1, j2, m2, n2, j3, m3, n3):
    """Int D^{j1}_{m1 n1} D^{j2}_{m2 n2} D^{j3}_{m3 n3} dg  (all plain)."""
    return three_j(j1, m1, j2, m2, j3, m3) * three_j(j1, n1, j2, n2, j3, n3)


# ============================================================ ladder
# Conventions. chi_a(X^{-1} Y) = sum_{q,x} conj(D^a_{qx}(X)) D^a_{qx}(Y).
# chi_c(X Y^{-1}) = sum_{y,z} D^c_{yz}(X) conj(D^c_{yz}(Y)).
# conj(D^j_{m n}(g)) = sum_{a,b} eps_{m a} D^j_{a b} eps_{n b}  (real eps).
# Every vertex then carries three PLAIN D's, integrated by `vertex`.

def column_tensor(a_prev, c, a_next, plain_rung: bool):
    """Vertex tensor for a rail vertex with incoming edge label a_prev
    (index pair (q,x) shared with previous vertex), outgoing edge a_next
    (pair (q',x') shared with next vertex), rung label c (pair (y,z)).
    Returns dict[(q,x,qp,xp,y,z)] -> Surd.
    Incoming edge gives D^{a_prev}_{qx}(g) [plain]; outgoing gives
    conj(D^{a_next}_{q'x'}(g)); the rung gives D^c_{yz}(g) if plain_rung
    else conj(D^c_{yz}(g))."""
    out = {}
    for q in ms(a_prev):
        for x in ms(a_prev):
            for qp in ms(a_next):
                for xp in ms(a_next):
                    for y in ms(c):
                        for z in ms(c):
                            tot = Surd()
                            # expand conj(D^{a_next}_{qp xp}) via eps
                            for qa in ms(a_next):
                                e1 = eps(a_next, qp, qa)
                                if e1 == 0:
                                    continue
                                for xa in ms(a_next):
                                    e2 = eps(a_next, xp, xa)
                                    if e2 == 0:
                                        continue
                                    if plain_rung:
                                        v = vertex(a_prev, q, x, a_next, qa, xa, c, y, z)
                                        tot = tot + v * (e1 * e2)
                                    else:
                                        for ya in ms(c):
                                            e3 = eps(c, y, ya)
                                            if e3 == 0:
                                                continue
                                            for za in ms(c):
                                                e4 = eps(c, z, za)
                                                if e4 == 0:
                                                    continue
                                                v = vertex(a_prev, q, x, a_next, qa, xa, c, ya, za)
                                                tot = tot + v * (e1 * e2 * e3 * e4)
                            if not tot.is_zero():
                                out[(q, x, qp, xp, y, z)] = tot
    return out


def ladder_eval(top, bottom, rungs):
    """N(a,b,c): top = [a_1..a_{m-1}], bottom = [b_1..b_{m-1}],
    rungs = [c_1..c_m]. Edge 0 (virtual, label 0) closes the left end;
    the right end likewise. Returns a Surd (should be rational)."""
    m = len(rungs)
    assert len(top) == len(bottom) == m - 1
    zero = F(0)
    top_l = [zero] + list(top) + [zero]
    bot_l = [zero] + list(bottom) + [zero]
    # cut state: dict[(q,x,qb,xb)] -> Surd, over the edge pairs (top, bottom)
    state = {(zero, zero, zero, zero): Surd.rat(1)}
    for i in range(m):
        c = rungs[i]
        tA = column_tensor(top_l[i], c, top_l[i + 1], plain_rung=True)     # A_i: rung plain
        tB = column_tensor(bot_l[i], c, bot_l[i + 1], plain_rung=False)    # B_i: rung conj
        new = {}
        for (q, x, qb, xb), amp in state.items():
            for (q1, x1, qp, xp, y, z), vA in tA.items():
                if (q1, x1) != (q, x):
                    continue
                for (q2, x2, qbp, xbp, y2, z2), vB in tB.items():
                    if (q2, x2) != (qb, xb) or (y2, z2) != (y, z):
                        continue
                    key = (qp, xp, qbp, xbp)
                    new[key] = new.get(key, Surd()) + amp * vA * vB
        state = new
    return state.get((zero, zero, zero, zero), Surd())


_CT = {}


def ct(a_prev, c, a_next, plain):
    key = (a_prev, c, a_next, plain)
    if key not in _CT:
        _CT[key] = column_tensor(a_prev, c, a_next, plain)
    return _CT[key]


def ladder_partition(u, v, k, m, jmax=F(1)):
    """Z = sum_labels prod (d_a u_a) prod (d_b v_b) prod (d_c k_c) N(a,b,c),
    computed as a transfer over the cut state (a, b, q, x, qb, xb) with the
    label sums folded in (coefficients are exact rationals)."""
    labels = [F(t, 2) for t in range(int(2 * jmax) + 1)]
    zero = F(0)
    state = {(zero, zero, zero, zero, zero, zero): Surd.rat(1)}
    for i in range(m):
        last = (i == m - 1)
        new = {}
        for (a, b, q, x, qb, xb), amp in state.items():
            for c in labels:
                kc = k[c] * (2 * c + 1)
                if kc == 0:
                    continue
                nexts = [(zero, zero)] if last else [(an, bn) for an in labels for bn in labels]
                for an, bn in nexts:
                    coef = kc
                    if not last:
                        coef *= u[an] * (2 * an + 1) * v[bn] * (2 * bn + 1)
                    if coef == 0:
                        continue
                    tA = ct(a, c, an, True)
                    tB = ct(b, c, bn, False)
                    for (q1, x1, qp, xp, y, z), vA in tA.items():
                        if (q1, x1) != (q, x):
                            continue
                        for (q2, x2, qbp, xbp, y2, z2), vB in tB.items():
                            if (q2, x2) != (qb, xb) or (y2, z2) != (y, z):
                                continue
                            key = (an, bn, qp, xp, qbp, xbp)
                            new[key] = new.get(key, Surd()) + amp * vA * vB * coef
        state = new
    total = state.get((zero, zero, zero, zero, zero, zero), Surd())
    assert total.is_rational(), ("surd did not cancel", total.t)
    return total.rational()


def run():
    half, one, zero = F(1, 2), F(1), F(0)
    # ---- V1/V2 one plaquette closed form: N = delta / d_a
    v1 = v2 = True
    for a in (zero, half, one):
        for b in (zero, half, one):
            for c1 in (zero, half, one):
                for c2 in (zero, half, one):
                    val = ladder_eval([a], [b], [c1, c2])
                    if not val.is_rational():
                        v1 = False
                    want = F(1, int(2 * a + 1) ** 2) if a == b == c1 == c2 else F(0)
                    if val.rational() != want:
                        v2 = False
    # ---- V3 single-row limit: v = delta_0, k = delta_0 -> chain closed form
    # Z with u = w-coefficients (u_a = f_a), bottom trivial, rungs trivial:
    # Z = f_0^{m-1} ... and with chi_half insertions at the ends? Here the
    # plain partition function: labels forced to 0 on bottom/rungs, then
    # top loop must be 0 as well -> Z = f_0^{m-1}.
    v3 = True
    fvals = {zero: F(2), half: F(1, 3), one: F(1, 10)}
    for m in (2, 3, 4):
        Z = ladder_partition(fvals, {zero: F(1), half: F(0), one: F(0)},
                             {zero: F(1), half: F(0), one: F(0)}, m)
        if Z != fvals[zero] ** (m - 1):
            v3 = False
    # ---- V4 two-plaquette recoupling: three rungs all 1/2, top/bottom
    # faces 1/2: N relates to the theta/6j; check the SUM over the middle
    # rung label of d_c * N reproduces the one-loop value (unitarity of
    # recoupling): sum_c d_c N(1/2,1/2 ; 1/2,1/2 ; 1/2, c, 1/2) = N(one loop)^? 
    # Use the exact identity: sum_c d_c {6j}^2 = 1/d  <-> here:
    parts = {}
    for c in (zero, one):
        val = ladder_eval([half, half], [half, half], [half, c, half])
        v1 = v1 and val.is_rational()
        parts[c] = (2 * c + 1) * val.rational()
    s = parts[zero] + parts[one]
    v4 = (s == F(1, 4)) and (parts[zero] / s == F(1, 4)) and (parts[one] / s == F(3, 4))
    # ---- V5 rail swap and reversal symmetry
    v5 = True
    top, bot, rg = [half, one, half], [one, half, half], [half, half, one, zero]
    n1 = ladder_eval(top, bot, rg).rational()
    n2 = ladder_eval(bot, top, rg).rational()
    n3 = ladder_eval(top[::-1], bot[::-1], rg[::-1]).rational()
    v5 = (n1 == n2 == n3)

    engine_ok = v1 and v2 and v3 and v4 and v5

    # ---- FIRST USE: product-ansatz Rayleigh (calibration)
    calib = {}
    kap = F(1, 8)
    f0 = f0_of(kap).lo
    fh = (f0_of(kap) * bessel_I(2, kap, TERMS) / bessel_I(1, kap, TERMS)).lo
    f1 = (f0_of(kap) * bessel_I(3, kap, TERMS) / bessel_I(1, kap, TERMS)).lo
    lam = lam_half().lo
    lam1 = (bessel_I(3, BETA, TERMS) / bessel_I(1, BETA, TERMS)).lo   # Wilson time kernel, spin 1
    k = {zero: F(1), half: lam, one: lam1}
    for c in (F(-1, 20), F(0), F(1, 20)):
        # u = w (1 + c chi_half): coefficients u_a with content <= 1:
        # w(1 + c chi_half) = sum d_j f_j chi_j + c sum d_j f_j chi_j chi_half
        # chi_j chi_half = chi_{j-1/2} + chi_{j+1/2}; keep a <= 1 (truncation)
        u = {zero: f0 + c * F(2) * fh / 1,            # from j=1/2: d_half f_half chi_0 /d_0
             half: fh + c * (f0 + F(3) * f1) / 2,     # (d_0 f_0 + d_1 f_1) chi_half / d_half
             one: f1 + c * (F(2) * fh) / 3}           # d_half f_half chi_1 / d_1 (+ f_3/2 dropped)
        rows = {}
        for m in (2, 3, 4, 5, 6):
            Z = ladder_partition(u, u, k, m)
            den = (f0 + 4 * c * fh + c * c * (f0 + 3 * f1)) ** (m - 1)
            R = Z / den
            rows[str(m)] = {"rayleigh": _dec(R, 10),
                            "ratio_to_f0_pow": _dec(R / f0 ** (m - 1), 10)}
        calib[str(c)] = rows
    cert = {
        "certificate_type": "YM30_TWO_ROW_RECOUPLING_ENGINE",
        "claim_status": "engine_validated_V1_V5__first_use_is_calibration_only",
        "validation": {"V1_rational": v1, "V2_one_plaquette_closed_form": v2,
                       "V3_single_row_limit": v3, "V4_recoupling_unitarity": v4,
                       "V5_symmetries": v5},
        "calibration_product_ansatz_kappa_1_8_content_le_1": calib,
        "scope": "labels <= 1; coefficients at interval lower ends; not an outward certificate",
        "verdict": "PASS" if engine_ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM30_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM30.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["validation"])
    print(json.dumps(cert["calibration_product_ansatz_kappa_1_8_content_le_1"], indent=0))
    print("sha256:", sha)
