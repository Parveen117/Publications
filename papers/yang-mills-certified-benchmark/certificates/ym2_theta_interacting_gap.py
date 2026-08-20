"""YM-2: Certified positive gap for the INTERACTING SU(2) theta-graph
transfer at certified small coupling.

Claim boundary (declared, fail-closed):
  - Carrier: theta graph, physical space H = L^2(SU(2)^2)^Ad, transfer
        T_kappa = M_kappa^{1/2} (K_beta (x) K_beta) M_kappa^{1/2},
    with K_beta the normalized Wilson convolution operator (YM-1 source) and
        m_kappa(A,B) = exp[(kappa/2)(Tr A + Tr B + Tr(A B^{-1}))]
    the three-face theta interaction (MP THETA_GRAPH_PROTOTYPE.md).
  - CERTIFIED here: for beta = 2 and every rational coupling
        0 <= kappa < kappa_0 := Delta_red(1,2) / 6,
    the interacting transfer has spectral-ratio gap
        lambda_2(T_kappa)/lambda_1(T_kappa) <= e^{6 kappa} * lambda_{1/2}(beta) < 1,
    equivalently reduced log-gap  Delta_theta(kappa) >= Delta_red - 6*kappa > 0.
    All quantities are exact two-sided rational enclosures.
  - NOT certified: large coupling, full lattice, UV/IR uniformity, continuum
    reconstruction, the Clay predicate. Adapter verdict stays `hold`.

Mathematical spine (each step elementary and stated in the certificate):
  S1  |Tr U| <= 2 on SU(2)  =>  s = TrA + TrB + Tr(AB^-1) in [-6, 6]
      =>  m_kappa in [e^{-3 kappa}, e^{3 kappa}] pointwise.
  S2  K_beta is PSD with character eigenvalues lambda_j = I_{2j+1}/I_1 > 0
      (positive Bessel series), so T_0 = K (x) K is PSD; its two top
      eigenvalues on the invariant sector are 1 and lambda_{1/2}.
  S3  Min-max sandwich: m_- T_0 <= T_kappa <= m_+ T_0 (operator order),
      hence lambda_k(T_kappa) in [m_- lambda_k(T_0), m_+ lambda_k(T_0)].
  S4  Ratio bound: lambda_2/lambda_1 (T_kappa) <= (m_+/m_-) lambda_{1/2}
      = e^{6 kappa} * I_2(beta)/I_1(beta).
  S5  Positive gap iff e^{6 kappa} < I_1/I_2, i.e. kappa < Delta_red/6.

Controls:
  C1  kappa = 0 recovers YM-1's bracket exactly (consumption check).
  C2  kappa above threshold: bound fails (ratio bracket not < 1) — the
      certificate refuses, fail-closed.
  C3  exp enclosure two-route: series vs squaring of half-argument overlap.
"""

from fractions import Fraction as F
import hashlib
import json
import os
import sys

sys.set_int_max_str_digits(200000)

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, log_iv, certified_reduced_gap, _dec, canonical_sha,
    TERMS, LOG_TERMS,
)

BETA = F(2)
EXP_TERMS = 60


# ------------------------------------------------------------ certified exp
def exp_point(x: F, terms=EXP_TERMS):
    """Two-sided bracket of e^x at rational x via Taylor + geometric tail."""
    neg = x < 0
    x = abs(x)
    # need tail ratio x/(k+1) < 1 eventually; terms chosen large enough
    t = F(1)
    s = F(1)
    for k in range(1, terms + 1):
        t = t * x / k
        s += t
    r = x / (terms + 1)
    if r >= 1:
        raise ValueError("increase terms")
    tail = t * r / (1 - r)
    lo, hi = s, s + tail
    if neg:
        lo, hi = 1 / hi, 1 / lo
    return Iv(lo, hi)


def exp_iv(x: Iv, terms=EXP_TERMS) -> Iv:
    return Iv(exp_point(x.lo, terms).lo, exp_point(x.hi, terms).hi)


# ---------------------------------------------------------------- YM-2 core
def certified_theta_gap(kappa: F, beta: F = BETA):
    """Certified bracket on the interacting reduced gap lower bound
    Delta_theta(kappa) >= Delta_red(beta) - 6*kappa, with the underlying
    ratio bound e^{6 kappa} * lambda_{1/2} bracketed exactly."""
    I1, I2, lam_half, gap_red = certified_reduced_gap(F(1), beta)
    ratio_bound = exp_iv(Iv(6 * kappa)) * lam_half         # S1+S4
    gap_lower = gap_red - Iv(6 * kappa)                    # S5 (log scale)
    positive = ratio_bound.hi < 1 and gap_lower.lo > 0
    return {
        "I1": I1, "I2": I2, "lam_half": lam_half,
        "gap_red": gap_red, "ratio_bound": ratio_bound,
        "gap_lower": gap_lower, "positive": positive,
    }


def kappa_threshold():
    """kappa_0 = Delta_red / 6, exact enclosure."""
    _, _, _, gap_red = certified_reduced_gap(F(1), BETA)
    return gap_red / Iv(6)


def run():
    k0 = kappa_threshold()

    # certified working coupling: a rational strictly below k0.lo
    kappa = F(1, 8)                                   # 0.125 < 0.13945...
    assert kappa < k0.lo
    main = certified_theta_gap(kappa)

    # C1: kappa = 0 consumption check — reduces exactly to YM-1
    free = certified_theta_gap(F(0))
    c1 = (free["ratio_bound"].lo == free["lam_half"].lo
          and free["ratio_bound"].hi == free["lam_half"].hi
          and free["gap_lower"].lo == free["gap_red"].lo)

    # C2: fail-closed above threshold — bound must NOT certify positivity
    above = certified_theta_gap(F(1, 4))              # 0.25 > kappa_0
    c2 = not above["positive"]

    # C3: exp two-route overlap at 6*kappa
    e1 = exp_point(6 * kappa)
    h = exp_point(3 * kappa)
    e2 = h * h
    c3 = not e1.separated_from(e2)

    ok = main["positive"] and c1 and c2 and c3
    cert = {
        "certificate_type": "YM2_CERTIFIED_INTERACTING_THETA_GAP_SMALL_COUPLING",
        "claim_status": "finite_cutoff_theta_graph_small_coupling",
        "claim_boundary": {
            "certified": ("theta-graph interacting transfer spectral-ratio gap "
                          "for 0 <= kappa < kappa_0 = Delta_red/6 at beta=2; "
                          "min-max sandwich m_- T_0 <= T_kappa <= m_+ T_0"),
            "not_certified": [
                "large coupling", "full 4d lattice", "UV/IR uniformity",
                "continuum OS reconstruction", "Clay mass-gap predicate",
            ],
            "sources": [
                "MP adapters/yang_mills/THETA_GRAPH_PROTOTYPE.md",
                "papers/yang-mills-certified-benchmark YM-1 (consumed)",
            ],
        },
        "parameters": {"beta": str(BETA), "kappa": str(F(1, 8)),
                       "exp_terms": EXP_TERMS, "bessel_terms": TERMS,
                       "log_terms": LOG_TERMS},
        "enclosures": {
            "kappa_threshold_lo": _dec(k0.lo, 40),
            "kappa_threshold_hi": _dec(k0.hi, 40),
            "ratio_bound_lo": _dec(main["ratio_bound"].lo, 40),
            "ratio_bound_hi": _dec(main["ratio_bound"].hi, 40),
            "gap_lower_lo": _dec(main["gap_lower"].lo, 40),
            "gap_lower_hi": _dec(main["gap_lower"].hi, 40),
        },
        "controls": {
            "C1_kappa0_recovers_YM1": bool(c1),
            "C2_above_threshold_fails_closed": bool(c2),
            "C3_exp_two_route_overlap": bool(c3),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM2_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM2.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"])
    print("kappa_0 in [", cert["enclosures"]["kappa_threshold_lo"], ",",
          cert["enclosures"]["kappa_threshold_hi"], "]")
    print("gap_lower(kappa=1/8) in [", cert["enclosures"]["gap_lower_lo"], ",",
          cert["enclosures"]["gap_lower_hi"], "]")
    print("sha256:", sha)
