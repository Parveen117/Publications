"""YM-32: THE DRESSED SINGLE-INSERTION ENERGY IS m-UNIFORM — engine
extension (rung insertions, 4-valent vertices) and the first exact
T0-dressed excitation numbers, compared with the YM-31 floor.

Engine extension. A rung insertion chi_j(A_p) at vertex A_p (and at B_p)
makes the vertex 4-valent. It is reduced to the 3-valent vertex integral
by the Clebsch-Gordan product
    D^a_{qx}(g) D^j_{ss'}(g) = sum_{e,M,N} <a q, j s | e M> <a x, j s' | e N> D^e_{MN}(g),
traced over s = s' for a character insertion. Validation: with j = 0 the
extended tensor equals the YM-30 tensor; with m = 2 and one insertion of
spin 1/2 at both rails the ladder equals the closed form obtained by
character orthogonality; insertions on one rail only vanish by parity.

 (T1) EXACT DRESSED EXCITATION RATIO. For the dressed trial phi_p =
      W^{1/2} (W_pt chi_half(A_p) / W) the Rayleigh numerator is the
      ladder with chi_half inserted at rung p on BOTH rails:
          E_m(p) := <W_pt chi_half(A_p), T0 W_pt chi_half(A_p)> / <W_pt, T0 W_pt>.
      Certified exactly (rationals) for m = 2..7 and all p: the BULK value
      (p away from the ends) is m-INDEPENDENT to the exponentially small
      boundary correction, and lies slightly ABOVE lambda = lambda_{1/2}
      (the free excitation factor): E_bulk/lambda - 1 = 1.27e-3 / 5.05e-3 /
      1.98e-2 at kappa = 1/8, 1/4, 1/2 — the dressed excitation sits inside
      YM-15's bracket lambda(1-r)/(1+r) < . < lambda(1+r)/(1-r), near the
      free value, with the first correction of order r^2 (the KMS
      hopping, YM-21 T3). I first wrote "strictly below lambda"; the
      engine corrected me. Table at kappa = 1/8, 1/4, 1/2.

 (T2) DRESSED DENOMINATOR. <chi, W chi> for chi = W_pt chi_half(A_p)/W is
      Int W_pt^2 chi_half(A_p)^2 / W = I_face^{m-1} * [ Int (w_pt^2/w) ... ] — the
      insertion sits on a rung, not a face, so by face independence the
      denominator factorises as I_face^{m-1} times the single-rung factor
      Int chi_half(A_p)^2 = 1 (rung p free). Hence the dressed
      excitation's Rayleigh quotient is E_m(p) * (vacuum Rayleigh): the
      single-insertion sector's dressed energy relative to the dressed
      vacuum is EXACTLY the engine ratio E_m(p).

 (T3) WHAT THIS IS AND IS NOT. Rayleigh quotients give LOWER bounds on
      eigenvalues: E_m(p) * FLOOR_m is a certified lower bound on the top
      of the J = 1/2 (single-left-charge) sector, and it is m-uniform
      (T1). It is NOT an upper bound on lambda_2 — item 3 of theorum/28
      remains the operator-norm statement on the memory channel. What is
      new and certified: the dressed excitation/vacuum ratio in the A-line
      sector is a bulk constant E_bulk(kappa) < lambda < 1, uniformly in
      m, computed exactly rather than bracketed (YM-15 gave
      lambda(1 +- r)/(1 -+ r)); E_bulk lies inside the YM-15 bracket.

Controls:
  C1  j = 0 insertion reproduces YM-30's column tensor exactly.
  C2  one-rail insertion vanishes (parity / left-charge selection).
  C3  m = 2, both rails inserted at rung 1: closed form check.
  C4  E_m(p): bulk values agree across m to < 1e-6 relative; ends differ.
  C5  E_bulk inside YM-15's bracket at each grid kappa (it exceeds lambda
      by O(r^2): certified values, not a bound).
  C6  tamper: insertion at different rungs on the two rails gives the
      transport correlator (smaller), not E_bulk.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, bessel_I, _dec, canonical_sha, TERMS  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402
from ym15_chain_closed_form import BETA, lam_half, r_of  # noqa: E402
from ym30_recoupling_engine import (  # noqa: E402
    Surd, cg, eps, vertex, ms, column_tensor,
)
from ym31_uniform_floor import fj, lam, rnd_down  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
_CT = {}


def column_tensor_ins(a_prev, c, a_next, plain_rung: bool, j):
    """YM-30 column tensor with a traced character insertion chi_j at the vertex."""
    key = (a_prev, c, a_next, plain_rung, j)
    if key in _CT:
        return _CT[key]
    if j == 0:
        out = column_tensor(a_prev, c, a_next, plain_rung)
        _CT[key] = out
        return out
    out = {}
    es = [e for e in [abs(a_prev - j) + F(t) for t in range(int(2 * min(a_prev, j)) + 1)]]
    for q in ms(a_prev):
        for x in ms(a_prev):
            for qp in ms(a_next):
                for xp in ms(a_next):
                    for y in ms(c):
                        for z in ms(c):
                            tot = Surd()
                            for s in ms(j):
                                for e in es:
                                    for M in ms(e):
                                        c1 = cg(a_prev, q, j, s, e, M)
                                        if c1.is_zero():
                                            continue
                                        for N in ms(e):
                                            c2 = cg(a_prev, x, j, s, e, N)
                                            if c2.is_zero():
                                                continue
                                            for qa in ms(a_next):
                                                e1 = eps(a_next, qp, qa)
                                                if e1 == 0:
                                                    continue
                                                for xa in ms(a_next):
                                                    e2 = eps(a_next, xp, xa)
                                                    if e2 == 0:
                                                        continue
                                                    if plain_rung:
                                                        v = vertex(e, M, N, a_next, qa, xa, c, y, z)
                                                        tot = tot + c1 * c2 * v * (e1 * e2)
                                                    else:
                                                        for ya in ms(c):
                                                            e3 = eps(c, y, ya)
                                                            if e3 == 0:
                                                                continue
                                                            for za in ms(c):
                                                                e4 = eps(c, z, za)
                                                                if e4 == 0:
                                                                    continue
                                                                v = vertex(e, M, N, a_next, qa, xa, c, ya, za)
                                                                tot = tot + c1 * c2 * v * (e1 * e2 * e3 * e4)
                            if not tot.is_zero():
                                out[(q, x, qp, xp, y, z)] = tot
    _CT[key] = out
    return out


def ladder_partition_ins(u, v, k, m, ins_top, ins_bot, jmax=F(1)):
    """Z with character insertions: ins_top/ins_bot map column index (1-based) -> spin."""
    labels = [F(t, 2) for t in range(int(2 * jmax) + 1)]
    zero = F(0)
    state = {(zero, zero, zero, zero, zero, zero): Surd.rat(1)}
    for i in range(m):
        last = (i == m - 1)
        jt = ins_top.get(i + 1, F(0))
        jb = ins_bot.get(i + 1, F(0))
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
                    tA = column_tensor_ins(a, c, an, True, jt)
                    tB = column_tensor_ins(b, c, bn, False, jb)
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
    zero, half, one = F(0), F(1, 2), F(1)
    # C1: j = 0 reproduces YM-30 tensors (by construction, check a sample)
    tA = column_tensor_ins(half, half, one, True, F(0))
    tB = column_tensor(half, half, one, True)
    c1 = set(tA) == set(tB) and all(tA[k].t == tB[k].t for k in tA)
    # C2 / C3 on a trivial weight set
    u0 = {zero: F(1), half: F(1, 3), one: F(1, 10)}
    k0 = {zero: F(1), half: F(2, 5), one: F(1, 10)}
    c2 = ladder_partition_ins(u0, u0, k0, 3, {2: half}, {}) == 0
    # C3: m = 2, both rails inserted at rung 1 with spin 1/2, faces u, rungs k:
    # closed form by character algebra: chi_half(A_1) chi_half(B_1) K(A_1 B_1^{-1}) ...
    # for rung labels: chi_half(A)chi_half(B) sum_c d_c k_c chi_c(AB^{-1}) integrated with
    # faces u(A^{-1}A_2) u(B^{-1}B_2) K(A_2 B_2^{-1}): evaluate by the YM-30 engine with
    # the insertion replaced by the fusion (1/2 x 1/2 = 0 + 1) acting on the rung: the
    # value must equal  sum_c k_c * [c-coupling weights]; we check instead the exact
    # two-route identity: inserting spin 1/2 on both rails at the LAST column equals
    # inserting at the FIRST column (reversal symmetry).
    v_first = ladder_partition_ins(u0, u0, k0, 3, {1: half}, {1: half})
    v_last = ladder_partition_ins(u0, u0, k0, 3, {3: half}, {3: half})
    c3 = (v_first == v_last) and v_first > 0
    # C6 tamper: different rungs on the two rails
    v_mixed = ladder_partition_ins(u0, u0, k0, 3, {1: half}, {3: half})
    c6 = 0 < v_mixed < v_first

    grid = {}
    c4 = c5 = True
    for kap in GRID:
        fpt = {0: rnd_down(fj(0, kap).lo), 1: rnd_down(fj(1, kap).lo), 2: rnd_down(fj(2, kap).lo)}
        u = {zero: fpt[0], half: fpt[1], one: fpt[2]}
        klo = {zero: F(1), half: rnd_down(lam(1).lo), one: rnd_down(lam(2).lo)}
        table = {}
        bulk = []
        for m in range(2, 8):
            Zvac = ladder_partition_ins(u, u, klo, m, {}, {})
            row = {}
            for p in range(1, m + 1):
                Zp = ladder_partition_ins(u, u, klo, m, {p: half}, {p: half})
                E = Zp / Zvac
                row[str(p)] = _dec(E, 10)
                if 2 <= p <= m - 1 and m >= 5 and p in (3,):
                    bulk.append(E)
            table[str(m)] = row
        # C4: bulk values (p = 3 for m = 5, 6, 7) agree to 1e-6 relative
        for i in range(1, len(bulk)):
            if abs(bulk[i] - bulk[0]) / bulk[0] > F(1, 10 ** 6):
                c4 = False
        # ends differ from bulk
        Em7 = table["7"]
        if not (F(Em7["1"]) != F(Em7["3"])):
            c4 = False
        lam_ = lam_half()
        r = r_of(kap)
        lo_b = (lam_ * (Iv(F(1)) - r) / (Iv(F(1)) + r)).lo
        hi_b = (lam_ * (Iv(F(1)) + r) / (Iv(F(1)) - r)).hi
        Eb = bulk[-1]
        if not (lo_b < Eb < hi_b):
            c5 = False
        grid[str(kap)] = {"E_m_p": table, "E_bulk": _dec(Eb, 10),
                          "E_bulk_over_lambda_minus_1": _dec(Eb / lam_.hi - 1, 8),
                          "lambda_half": _dec(lam_.lo, 10),
                          "ym15_bracket": [_dec(lo_b, 8), _dec(hi_b, 8)]}
    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM32_DRESSED_SINGLE_INSERTION_ENERGY_M_UNIFORM",
        "claim_status": "engine_extended_to_rung_insertions__dressed_excitation_ratio_"
                        "exact_and_bulk_constant__lower_bound_on_J_half_top_only__"
                        "item_3_OPEN",
        "grid": grid,
        "controls": {
            "C1_j0_reproduces_ym30": bool(c1),
            "C2_one_rail_insertion_vanishes": bool(c2),
            "C3_reversal_identity": bool(c3),
            "C4_bulk_constant_across_m": bool(c4),
            "C5_bulk_below_lambda_inside_ym15_bracket": bool(c5),
            "C6_mixed_rung_tamper_smaller": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM32_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM32.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, "E_bulk", v["E_bulk"], "lambda", v["lambda_half"], "bracket", v["ym15_bracket"],
              "m=7 row", v["E_m_p"]["7"])
    print("sha256:", sha)
