"""YM-5: CERTIFIED TWO-SIDED SPECTRAL GAP of the interacting theta-graph
transfer at finite coupling — complement control via the exact doubling
identity m_kappa^2 = m_{2 kappa}.

Lineage note (mine, don't merge): the old lambda-framework manuscripts
(Yang_Mill_1/2.tex, uncertified) reduce "positivity + Schur off-block
estimate => gap" (their Kill Lemma) conditionally on a conjectured vacuum.
This certificate realizes exactly that skeleton on the theta-graph carrier,
unconditionally and machine-checked, with every constant an exact enclosure.
No claim from those manuscripts is imported. See LINEAGE.md.

Claim boundary (declared, fail-closed):
  Carrier as in YM-2/3/4, S_kappa = T0^{1/2} M_kappa T0^{1/2} (isospectral
  to T_kappa), P = projection onto V = span{1, chi12(A), chi12(B)} — exact
  T0 eigenvectors (eigenvalues 1, lam, lam; lam = lambda_{1/2}) — so P and
  Q = I - P commute with T0.

  CERTIFIED, for each grid kappa where the verdict fires:
    (U1) block bound   |Q S Q| <= lam_next * e^{3 kappa},
         lam_next = lambda_{1/2}^2 (certified next free level in the
         Ad-invariant sector: lambda_1, lambda_{3/2} certified below it),
    (U2) off-block bound via EXACT DOUBLING m_kappa^2 = m_{2 kappa}:
             |P M Q|^2 <= lammax( P m_{2k} P - (P m_k P)^2 )
         with the 3x3 interval matrices from YM-4's certified compression
         and a Gershgorin lammax bound,
             |P S Q| <= sqrt(lam_next) * |P M Q|,
    (U3) Weyl:  lambda_2(T_kappa) <= max( lambda_2(P S P), |Q S Q| ) + |P S Q|,
    (L)  interlacing:  lambda_1(T_kappa) >= lambda_1(P S P),
    =>  certified gap ratio upper bound
             R(kappa) := lambda_2_upper / lambda_1_lower,
        and where R(kappa) < 1 the interacting transfer has a certified
        TWO-SIDED positive spectral gap:
             Delta(kappa) >= -log R(kappa) > 0.

  HEADLINE TARGET: certified positive gap at kappa BEYOND YM-2's sandwich
  threshold kappa_0 = 0.1394..., i.e. the Kill-Lemma route strictly
  outperforms the min-max sandwich.

  NOT certified: couplings where R >= 1 (certificate refuses, fail-closed),
  full 4d lattice, UV/IR limits, continuum OS reconstruction, the vacuum
  construction, the Clay predicate. Adapter verdict stays `hold`.

Controls:
  C1  doubling identity content check: P m_{2k} P - (P m_k P)^2 is a Gram
      matrix of (I-P) M phi_i, hence PSD; the certificate checks its crude
      Gershgorin lower bound stays above -1/25 (Gershgorin subtracts full
      off-diagonal rows, so a small genuine PSD matrix with comparable
      off-diagonals reads slightly negative — the slack is documented, and
      soundness never uses positivity of D, only its Gershgorin UPPER
      bound).
  C2  kappa -> 0: R reproduces the free ratio lambda_{1/2} (bracket overlap
      at tiny kappa within slack).
  C3  next-level certification: lambda_1, lambda_{3/2} strictly below
      lambda_{1/2}^2 (so lam_next is the true Q-top of the free transfer
      restricted to the invariant sector's low levels), plus the declared
      structural fact that T0 <= lam_next on ALL of Q (content pairs (j1,j2)
      != (0,0),(1/2,0),(0,1/2) have lambda_{j1} lambda_{j2} <= lam_next).
  C4  fail-closed demo: at large kappa (kappa = 1) the bound must NOT
      certify (R >= 1 there with these blocks) — refusal recorded.
"""

from fractions import Fraction as F
import hashlib
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS,
)
from ym2_theta_interacting_gap import exp_point  # noqa: E402
from ym4_symmetry_protected import (  # noqa: E402
    m_matrix as _m_matrix_raw, iv_sqrt as _iv_sqrt_raw, sym2_eigs, BETA,
)

_ROUND = 10 ** 30


def _r(x: Iv) -> Iv:
    """Outward rounding to denominator 1e30 — sound (bracket only widens)."""
    import math
    lo = F(math.floor(x.lo * _ROUND), _ROUND)
    hi = F(math.ceil(x.hi * _ROUND), _ROUND)
    return Iv(lo, hi)


_MCACHE = {}


def m_matrix(kappa, **kw):
    key = (kappa, tuple(sorted(kw.items())))
    if key not in _MCACHE:
        M, err = _m_matrix_raw(kappa, **kw)
        _MCACHE[key] = ([[_r(M[i][j]) for j in range(3)] for i in range(3)],
                        err)
    return _MCACHE[key]


def iv_sqrt(x: Iv, tol=F(1, 10 ** 20)) -> Iv:
    return _r(_iv_sqrt_raw(_r(x), tol=tol))

GRID5 = [F(1, 8), F(7, 50), F(1, 5), F(1, 4), F(3, 10)]
KAPPA_FAILDEMO = F(1)


def lam_levels():
    I1 = bessel_I(1, BETA, TERMS)
    lam = bessel_I(2, BETA, TERMS) / I1          # lambda_{1/2}
    lam1 = bessel_I(3, BETA, TERMS) / I1         # lambda_1
    lam32 = bessel_I(4, BETA, TERMS) / I1        # lambda_{3/2}
    return lam, lam1, lam32


def gershgorin_max(M3):
    """Certified upper bound on lammax of a symmetric 3x3 interval matrix."""
    bound = None
    for i in range(3):
        row = M3[i][i].hi
        for j in range(3):
            if j != i:
                row += max(abs(M3[i][j].lo), abs(M3[i][j].hi))
        bound = row if bound is None else max(bound, row)
    return bound


def gershgorin_min(M3):
    bound = None
    for i in range(3):
        row = M3[i][i].lo
        for j in range(3):
            if j != i:
                row -= max(abs(M3[i][j].lo), abs(M3[i][j].hi))
        bound = row if bound is None else min(bound, row)
    return bound


def psp_spectrum(kappa: F, lam: Iv):
    """Sector-resolved certified spectrum of P S_kappa P (YM-4 machinery)."""
    M, err = m_matrix(kappa)
    s_half = iv_sqrt(lam)
    a = M[0][0]
    b = ((M[0][1] + M[0][2]) * iv_sqrt(Iv(F(1, 2)))) * s_half
    c = (M[1][1] + M[1][2]) * lam
    ev_top, ev_second = sym2_eigs(a, b, c)
    odd = (M[1][1] - M[1][2]) * lam
    lam1_lower = ev_top.lo                       # interlacing lower bound
    lam2_psp_upper = max(ev_second.hi, odd.hi)
    return M, lam1_lower, lam2_psp_upper


def offblock_bound(kappa: F, lam_next_hi: F):
    """|P S Q| <= sqrt(lam_next) * sqrt( lammax(P m_{2k} P - (P m_k P)^2) )."""
    Mk, _ = m_matrix(kappa)
    M2k, _ = m_matrix(2 * kappa)
    # (P m_k P)^2 as interval matrix product
    sq = [[Iv(F(0)) for _ in range(3)] for _ in range(3)]
    for i in range(3):
        for j in range(3):
            acc = Iv(F(0))
            for k in range(3):
                acc = acc + Mk[i][k] * Mk[k][j]
            sq[i][j] = acc
    D = [[M2k[i][j] - sq[i][j] for j in range(3)] for i in range(3)]
    gmax = gershgorin_max(D)
    gmin = gershgorin_min(D)
    pmq = iv_sqrt(Iv(F(0), max(F(0), gmax))).hi
    psq = iv_sqrt(Iv(lam_next_hi)).hi * pmq
    return psq, gmin, gmax


def certified_two_sided(kappa: F):
    lam, lam1_lvl, lam32_lvl = lam_levels()
    lam_next = lam * lam                          # lambda_{1/2}^2
    # C3 level certification
    levels_ok = (lam1_lvl.hi < lam_next.lo and lam32_lvl.hi < lam_next.lo
                 and lam_next.hi < lam.lo)
    _, lam1_lower, lam2_psp_upper = psp_spectrum(kappa, lam)
    e3k_hi = exp_point(3 * kappa).hi
    qsq_upper = lam_next.hi * e3k_hi              # U1
    psq_upper, gmin, _ = offblock_bound(kappa, lam_next.hi)   # U2
    lam2_upper = max(lam2_psp_upper, qsq_upper) + psq_upper   # U3
    ratio_hi = lam2_upper / lam1_lower
    certifies = levels_ok and ratio_hi < 1 and lam1_lower > 0
    gap_lower = None
    if certifies:
        gap_lower = -(log_iv(Iv(ratio_hi), LOG_TERMS)).hi     # careful sign
        gap_lower = -(log_iv(Iv(ratio_hi, ratio_hi), LOG_TERMS).hi)
    return {
        "lam1_lower": lam1_lower, "lam2_psp_upper": lam2_psp_upper,
        "qsq_upper": qsq_upper, "psq_upper": psq_upper,
        "lam2_upper": lam2_upper, "ratio_hi": ratio_hi,
        "gram_gershgorin_min": gmin, "levels_ok": levels_ok,
        "certifies": certifies, "gap_lower": gap_lower,
    }


def run():
    lam, _, _ = lam_levels()
    grid_out = {}
    fired_beyond_threshold = False
    all_levels_ok = True
    c1_all = True
    for kap in GRID5:
        r = certified_two_sided(kap)
        all_levels_ok = all_levels_ok and r["levels_ok"]
        c1_all = c1_all and (r["gram_gershgorin_min"] > F(-1, 25))
        grid_out[str(kap)] = {
            "lambda1_lower": _dec(r["lam1_lower"], 25),
            "lambda2_upper": _dec(r["lam2_upper"], 25),
            "ratio_upper": _dec(r["ratio_hi"], 25),
            "two_sided_gap_certified": bool(r["certifies"]),
            "gap_lower": (_dec(r["gap_lower"], 25)
                          if r["gap_lower"] is not None else None),
        }
        if r["certifies"] and kap > F(1394538843, 10 ** 10):
            fired_beyond_threshold = True

    # C2 tiny-kappa reproduction
    tiny = certified_two_sided(F(1, 1000))
    c2 = (tiny["certifies"]
          and abs(tiny["ratio_hi"] - lam.hi) < F(1, 10))

    # C4 fail-closed at kappa = 1
    big = certified_two_sided(KAPPA_FAILDEMO)
    c4 = not big["certifies"]

    ok = all_levels_ok and c1_all and c2 and c4 and fired_beyond_threshold
    cert = {
        "certificate_type": "YM5_TWO_SIDED_GAP_KILL_LEMMA_REALIZATION",
        "claim_status": "certified_two_sided_gap_finite_coupling",
        "claim_boundary": {
            "certified": [
                "two-sided spectral gap of the interacting theta transfer "
                "at every grid kappa marked certified, via U1-U3 + L",
                "off-block control from the exact doubling m_k^2 = m_{2k}",
                "gap certified BEYOND YM-2's sandwich threshold kappa_0",
            ],
            "not_certified": [
                "couplings where the ratio bound >= 1 (refused)",
                "full 4d lattice", "UV/IR limits", "continuum OS",
                "vacuum construction", "Clay mass-gap predicate",
            ],
            "lineage": "lambda-framework Kill-Lemma skeleton realized; "
                       "no claim imported (see LINEAGE.md)",
            "sources": ["YM-1/2/3/4 (consumed, pinned)"],
        },
        "parameters": {"beta": str(BETA),
                       "grid": [str(x) for x in GRID5],
                       "fail_demo_kappa": str(KAPPA_FAILDEMO)},
        "grid_results": grid_out,
        "controls": {
            "C1_gram_difference_psd_within_tolerance": bool(c1_all),
            "C2_tiny_kappa_reproduces_free_ratio": bool(c2),
            "C3_next_level_certified": bool(all_levels_ok),
            "C4_large_kappa_fails_closed": bool(c4),
            "fired_beyond_sandwich_threshold": bool(fired_beyond_threshold),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


def run_all_and_pin():
    import ym1_certified_gap as ym1
    import ym2_theta_interacting_gap as ym2
    import ym3_crossing_direction as ym3
    import ym4_symmetry_protected as ym4mod
    out = {}
    for name, mod_run in (("YM1", ym1.run), ("YM2", ym2.run),
                          ("YM3", ym3.run), ("YM4", ym4mod.run),
                          ("YM5", run)):
        cert = mod_run()
        sha = canonical_sha(cert)
        with open(os.path.join(HERE, f"{name}_RESULT.json"), "w") as fjson:
            json.dump(cert, fjson, indent=2, sort_keys=True)
        with open(os.path.join(HERE, f"EXPECTED_{name}.sha256"), "w") as fsha:
            fsha.write(sha + "\n")
        out[name] = (cert["verdict"], sha)
        print(f"{name}: {cert['verdict']}  sha256:{sha[:16]}...")
    return out


if __name__ == "__main__":
    results = run_all_and_pin()
    assert all(v == "PASS" for v, _ in results.values())
    print("ALL CERTIFICATES PASS")
