"""YM-41: THE TOOLS DOCK — the operator-branch tools (RKF theorum/74, 75, 76)
consumed into the YM program, re-verified on the Publications carrier, and
E4D-C's remainder restated as what it has become.

Cross-repo consumption (YM-27 discipline: read, cite, re-verify on this
side; never rederive by narration): RKF branch theorem-49-local-to-uniform-
seam-gap, commits bfe40c2 (theorum/74 lopa-ledger contraction), d185e86 +
13b4315 (theorum/75 selective release, native regrading), 51afbcd
(theorum/76 seam-flow meter), audit 3534418/8679f0b. Publications cannot
fetch that repo in CI, so every number used here is RE-DERIVED below on
Publications' own pencil (ym37, pin ad598ee5...) with Publications' own
rounded Bessel enclosures; the RKF pins are recorded as the source
derivations.

 (T1) SILENCE-CHANNEL CONTRACTION on the YM pencil (theorum/74's tool,
      re-verified): with w = the silence vector e0 (the structural
      content-(0,0) line — NO Krylov iterate, defined on every carrier):
      defect exactly B-orthogonal, restricted doubled ceiling sigma_sil
      by one native elimination sign check, outward certificate
      beta_sil < 1 at every grid kappa. Non-vacuous refusal control
      (audit lesson): a heavy-complement pencil must be refused.

 (T2) LADDER LAW on Publications' own Bessel (theorum/75 T1,
      re-verified): f_{c+1/2}/f_c <= kappa/(2(2c+2)) for c = 0..3 on the
      enclosures, with the DISCRIMINATION check that the next sharper
      claim kappa/(2(2c+3)) FAILS — the law is tight, not loose. The
      full selective-release theorem (stationarity, seam-count
      retention at Lambda = 3/2, 2, outward certificate u + e < 1 for
      the untruncated column) is CONSUMED from theorum/75 by pin; its
      Publications-side re-instantiation is named YM-42 work, not
      claimed here.

 (T3) SEAM-FLOW METER on the YM pencil (theorum/76's instrument,
      re-verified): M is kappa-independent — the coupling path is lopa-
      ladder arithmetic on precomputed content-pair overlaps (native
      character projectors). Sector flips of k(kappa, mu) bracketed to
      2^-40 with recursive crossing separation; per bracket the count
      route and the native determinant route (T53 T6 product of
      weights) must AGREE (odd jump <=> sign flip). This is the
      instrument that will bracket the program's coupling thresholds
      exactly; here it is certified on the column pencil (honest label:
      the column pencil, NOT the chain carrier of YM-6/7).

 (T4) E4D-C RESTATED. Ledger update: with
        - m eliminated (YM-37: every quantity is one column),
        - the content tail released lawfully (theorum/75: outward
          certificate for the FULL untruncated column, ladder-justified
          declared tail beyond content 2),
        - the contraction structural (theorum/74: silence channel, no
          spectral objects, no Krylov),
      the alternating product K^{1/2} M_w K^{1/2}'s remaining open
      content is EXACTLY the t-direction: whether the time gap
      extracted from the t-rail column family {tau_t} persists
      t-uniformly. Named next: YM-42 = beta_sil(t = 3) on the 52-dim
      even sector (the t-ladder of silence contractions), and the
      t-uniform statement linking beta_sil(t) to the excitation/vacuum
      ratio limit. Nothing else remains of the m- or content-direction.
      Clay/OS/continuum gates: unchanged, OPEN.

Controls:
  C1  beta_sil < 1 at every kappa; r_sil'Bw = 0 exact; the heavy-
      complement refusal pencil is refused (control can fail).
  C2  ladder law holds AND the sharper claim fails (discrimination).
  C3  meter: sectors constant on certified intervals; every bracket's
      parity law (count vs det) holds; close pairs separated.
  C4  ladder identity: B rebuilt from content-pair overlaps equals the
      ym37 pencil B at 3 spot kappas, entry by entry, exactly.
  C5  tamper: moving (1/2,1/2) overlap mass into the (1,1) slot changes
      B (the (0,1/2) slot is EMPTY — centre superselection inside the
      meter, recorded as in theorum/76).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(1000000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import _dec, canonical_sha  # noqa: E402
from ym31_uniform_floor import fj, lam, rnd_down  # noqa: E402
import ym37_space_transfer as y37  # noqa: E402
import ym38_native_uniformity as y38  # noqa: E402

GRID = [F(1, 8), F(1, 4), F(1, 2)]
ZERO, HALF, ONE = F(0), F(1, 2), F(1)
CONTENTS = [ZERO, HALF, ONE]

RKF_SOURCES = {
    "theorum_74_lopa_ledger_contraction": "bfe40c2",
    "theorum_75_selective_release": "d185e86",
    "theorum_75_native_regrading": "13b4315",
    "theorum_76_seam_flow_meter": "51afbcd",
    "audit_74_75": "3534418+8679f0b",
}


def chi_of_u(c):
    polys = {F(0): [F(1)], HALF: [F(0), F(1)]}
    cur = HALF
    while cur < c:
        nxt = cur + HALF
        a, b = polys[cur], polys[cur - HALF]
        na = [F(0)] + a
        nb = b + [F(0)] * (len(na) - len(b))
        polys[nxt] = [x - y for x, y in zip(na, nb)]
        cur = nxt
    return polys[c]


def char_component(mono4, c, t=2):
    """content-c component of a rail monomial via the native character ladder
    (RKF theorum/75 regrading, re-implemented on Publications' engine)."""
    u = y37.padd({}, y37.dot(0, 1, t), F(2))
    coeffs = chi_of_u(c)
    chi = {}
    upow = y37.pconst(1, t)
    for k, a in enumerate(coeffs):
        if a != 0:
            chi = y37.padd(chi, upow, a)
        upow = y37.pmul(upow, u)
    my = {(0, 0, 0, 0) + mono4: F(1)}
    prod = y37.pmul(chi, my)
    out = {}
    for e, cf in prod.items():
        m = y37.moment4(e[4:8])
        if m == 0:
            continue
        out[e[0:4]] = out.get(e[0:4], F(0)) + cf * m * (2 * c + 1)
    return {e: v for e, v in out.items() if v != 0}


def pencil_and_overlaps(kap):
    fpt = {ZERO: rnd_down(fj(0, kap).lo), HALF: rnd_down(fj(1, kap).lo), ONE: rnd_down(fj(2, kap).lo)}
    klo = {ZERO: F(1), HALF: rnd_down(lam(1).lo), ONE: rnd_down(lam(2).lo)}
    basis, M, N, B, e0 = y37.build_pencil(2, fpt, klo)
    return basis, M, N, B, e0, fpt, klo


_S_CACHE = {}


def overlaps(basis):
    key = id(basis)
    if key in _S_CACHE:
        return _S_CACHE[key]
    n = len(basis)
    comps = []
    for b in basis:
        d = {}
        for e, cf in b.items():
            for c1 in CONTENTS:
                left = char_component(e[0:4], c1)
                if not left:
                    continue
                for c2 in CONTENTS:
                    right = char_component(e[4:8], c2)
                    if not right:
                        continue
                    tgt = d.setdefault((c1, c2), {})
                    for e1, v1 in left.items():
                        for e2, v2 in right.items():
                            k = e1 + e2
                            tgt[k] = tgt.get(k, F(0)) + cf * v1 * v2
        comps.append(d)
    S = {}
    for i in range(n):
        for j in range(n):
            for cc, poly in comps[j].items():
                v = y37.integrate(y37.pmul(basis[i], poly), 2)
                if v != 0:
                    S[(i, j, cc)] = v
    _S_CACHE[key] = S
    return S


def B_of(kap, n, S):
    fpt = {ZERO: rnd_down(fj(0, kap).lo), HALF: rnd_down(fj(1, kap).lo), ONE: rnd_down(fj(2, kap).lo)}
    B = [[F(0)] * n for _ in range(n)]
    for (i, j, (c1, c2)), v in S.items():
        B[i][j] += v / (fpt[c1] * fpt[c2])
    return B


def count(M, B, mu):
    n = len(M)
    return y37.inertia([[M[i][j] - mu * B[i][j] for j in range(n)] for i in range(n)])[0]


def det_weights(M, B, mu):
    n = len(M)
    A = [[M[i][j] - mu * B[i][j] for j in range(n)] for i in range(n)]
    detv = F(1)
    active = list(range(n))
    while active:
        p = next((i for i in active if A[i][i] != 0), None)
        if p is None:
            return F(0)
        d = A[p][p]
        detv *= d
        rest = [k for k in active if k != p]
        for r in rest:
            if A[r][p] != 0:
                fq = A[r][p] / d
                for c in rest:
                    A[r][c] -= fq * A[p][c]
        active = rest
    return detv


def run():
    grid = {}
    C = [True] * 5
    meter = {}
    for kap in GRID:
        basis, M, N, B, e0i, fpt, klo = pencil_and_overlaps(kap)
        n = len(basis)
        # ---- T1 silence contraction ----
        w = [F(0)] * n
        w[e0i] = F(1)
        W_w = B[e0i][e0i]
        tw = y37.tau_apply(M, B, w)
        theta = y38.bq(B, w, tw) / W_w
        r_sil = [tw[k] - theta * w[k] for k in range(n)]
        c1 = (y38.bq(B, r_sil, w) == 0)
        Rhat = y38.sqrt_up(y38.bq(B, r_sil, r_sil))
        U, Mp, Bp = y38.restricted(M, B, w, W_w)
        sig = y38.doubled_ceiling(Mp, Bp, hi=F(1))
        xs = [w]
        for _ in range(8):
            xs.append(y37.tau_apply(M, B, xs[-1]))
        t_up = []
        for x in xs:
            b = y38.bq(B, w, x) / W_w
            yv = [x[k] - b * w[k] for k in range(n)]
            t_up.append(y38.sqrt_up(y38.bq(B, yv, yv)) / b)
        tbar = max(t_up[2:]) * F(3, 2)
        th_lo = theta - Rhat * tbar / W_w
        beta = sig / th_lo + (theta + Rhat * tbar / W_w) * tbar * Rhat / (W_w * th_lo * th_lo)
        inv_ok = th_lo > 0 and (Rhat + sig * tbar) / th_lo <= tbar
        c1 = c1 and inv_ok and beta < 1
        # non-vacuous refusal (audit lesson)
        Mref = [[F(1), F(0)], [F(0), F(5)]]
        Bref = [[F(1), F(0)], [F(0), F(1)]]
        _, Mpr, Bpr = y38.restricted(Mref, Bref, [F(1), F(0)], F(1))
        c1 = c1 and (y38.doubled_ceiling(Mpr, Bpr, hi=F(8), steps=30) >= 1)
        C[0] = C[0] and c1
        # ---- T2 ladder law + discrimination ----
        c2 = True
        for c2i in range(0, 4):
            lo0 = fj(c2i, kap).lo
            hi1 = fj(c2i + 1, kap).hi
            if not (hi1 / lo0 <= kap / (2 * (c2i + 2))):
                c2 = False
            if hi1 / lo0 <= kap / (2 * (2 * c2i + 3)):
                c2 = False
        C[1] = C[1] and c2
        grid[str(kap)] = {"theta_sil": _dec(theta, 16), "sigma_sil": _dec(sig, 14),
                          "beta_sil": _dec(beta, 12), "beta_lt_1": bool(beta < 1)}
        # ---- C4/C5 ladder identity + tamper (once, at each kappa) ----
        S = overlaps(basis)
        c4 = (B_of(kap, n, S) == B)
        C[3] = C[3] and c4
        if kap == GRID[0]:
            assert not any(cc == (ZERO, HALF) for (_, _, cc) in S), "superselection"
            Sbad = {}
            for (i, j, cc), v in S.items():
                key = (i, j, (ONE, ONE)) if cc == (HALF, HALF) else (i, j, cc)
                Sbad[key] = Sbad.get(key, F(0)) + v
            C[4] = C[4] and (B_of(kap, n, Sbad) != B_of(kap, n, S))
            # ---- T3 meter scan (structural: M fixed, B(kappa) by ladder) ----
            mu = F(1, 2)
            KG = [F(k, 16) for k in range(1, 65)]
            counts = [count(M, B_of(k2, n, S), mu) for k2 in KG]
            flips = []
            c3 = True

            def bisect_all(lo, hi, kl, kh, depth=0):
                nonlocal c3
                if kh == kl:
                    return
                if hi - lo <= F(1, 2 ** 40) or depth >= 60:
                    dlo = det_weights(M, B_of(lo, n, S), mu)
                    dhi = det_weights(M, B_of(hi, n, S), mu)
                    flip = (dlo > 0) != (dhi > 0)
                    if (abs(kh - kl) % 2 == 1) != flip:
                        c3 = False
                    flips.append({"kappa_lo": str(lo), "jump": kh - kl, "det_flip": bool(flip)})
                    return
                mid = (lo + hi) / 2
                km = count(M, B_of(mid, n, S), mu)
                bisect_all(lo, mid, kl, km, depth + 1)
                bisect_all(mid, hi, km, kh, depth + 1)

            for t in range(len(KG) - 1):
                if counts[t] != counts[t + 1]:
                    bisect_all(KG[t], KG[t + 1], counts[t], counts[t + 1])
                else:
                    midk = (KG[t] + KG[t + 1]) / 2
                    if count(M, B_of(midk, n, S), mu) != counts[t]:
                        c3 = False
            C[2] = C[2] and c3
            meter = {"mu": str(mu), "kappa_range": [str(KG[0]), str(KG[-1])],
                     "sectors": sorted(set(counts)), "n_flips": len(flips), "flips": flips}
    ok = all(C)
    cert = {
        "certificate_type": "YM41_TOOLS_DOCK_E4DC_RESTATED",
        "claim_status": "silence_contraction_reverified_on_ym_pencil__ladder_law_reverified_with_discrimination__"
                        "seam_flow_meter_reverified__selective_release_CONSUMED_by_pin__"
                        "E4DC_remainder_is_the_t_direction_ONLY__ym42_named",
        "rkf_source_pins": RKF_SOURCES,
        "grid": grid,
        "meter_column_pencil": meter,
        "controls": {"C1_silence_contraction_and_nonvacuous_refusal": bool(C[0]),
                     "C2_ladder_law_and_discrimination": bool(C[1]),
                     "C3_meter_parity_per_bracket": bool(C[2]),
                     "C4_ladder_identity_exact": bool(C[3]),
                     "C5_superselection_and_tamper": bool(C[4])},
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM41_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    with open(os.path.join(HERE, "EXPECTED_YM41.sha256"), "w") as f:
        f.write(canonical_sha(cert) + "\n")
    print("verdict:", cert["verdict"])
    print(json.dumps(cert["controls"], indent=1))
    for k, v in cert["grid"].items():
        print(k, "beta_sil", v["beta_sil"], v["beta_lt_1"])
    print("meter:", cert["meter_column_pencil"]["sectors"], cert["meter_column_pencil"]["n_flips"])
