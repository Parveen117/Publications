"""YM-15: THE 1D BLOCK-TRANSFER DOCK, OPENED — exact closed form of the
bounded-overlap chain compression for EVERY chain length m, and the
first m-UNIFORM statement in the overlap regime.

This is YM-14's named object (T4). Carrier: the chain S_m of m blocks
with holonomies (A_i, B_i), i = 1..m, nearest-neighbour bridges
    m_kappa = prod_{i=1}^{m-1} exp[ (kappa/2) Tr(A_i A_{i+1}^{-1}) ],
free transfer T0 with reduced level lambda = lambda_{1/2} on every
chi12(A_i) (YM-14 convention, beta = 2). Compression onto the A-LINE
carrier
    W_{m+1} = { 1, chi12(A_1), ..., chi12(A_m) }
(orthonormal, all exact T0 eigenvectors; B-contents exactly inert by
parity, as in YM-14 C2).

Write the single-bridge character expansion
    exp[(kappa/2) Tr U] = sum_j d_j f_j chi_j(U),   f_j = 2 I_{2j+1}(kappa)/kappa
and r := f_{1/2} / f_0 = I_2(kappa) / I_1(kappa)  (0 < r < 1 for kappa > 0).

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) CLOSED FORM FOR ALL m — THE CONVOLUTION LEMMA. The compressed
      bridge product on W_{m+1} is block-diagonal:
          vacuum entry      <1, m_kappa 1>                = f_0^{m-1}
          A-line block      <chi12(A_p), m_kappa chi12(A_q)> = f_0^{m-1} r^{|p-q|}
          vacuum/A coupling                                 = 0 (parity).
      The A-line block is f_0^{m-1} times the Kac-Murdock-Szego matrix
      KMS_m(r) = [ r^{|p-q|} ].  PROOF (induction along the chain, the
      only integral used is the SU(2) convolution
          Int dA chi_a(X A^{-1}) chi_b(A Y) = delta_{ab} chi_a(XY) / d_a
      plus Int dA chi_c(A) = delta_{c0}): integrating the sites left to
      right, a site carrying content c and one bridge to its right emits
      exactly f_c chi_c at the next site; the vacuum carries c = 0
      through every bridge (f_0 each); an insertion chi12(A_p) turns the
      carried content into 1/2 until the second insertion at A_q returns
      it to chi_0 + chi_1, whose chi_1 branch dies at the open end
      (Int chi_1 = 0) and whose chi_0 branch carries f_0 to the end.
      Hence f_0^{(p-1)} f_{1/2}^{(q-p)} f_0^{(m-q)} = f_0^{m-1} r^{q-p}.
      No truncation remainder exists on this carrier: every f_j with
      j > 1 is integrated to zero EXACTLY (the dead chi_1 branch is the
      only j = 1 appearance). The capsule MACHINE-CHECKS this identity
      as an exact polynomial identity in formal indeterminates
      (f_0, f_{1/2}, f_1, ...) for m = 2..8 and all matrix entries, via
      an independent symbolic chain-integration engine (the
      "convolution formula" realised as code, not as numbers).

 (T2) m-UNIFORM BRACKET ON THE A-TOP / VACUUM RATIO. With T0 weights,
      the compressed transfer on W_{m+1} is
          diag( f_0^{m-1} ,  lambda f_0^{m-1} KMS_m(r) ).
      Normalised to the vacuum (the YM-5/YM-6 ratio convention:
      rho = lambda_top-excited / lambda_vacuum, gap = -log rho), the
      A-line top satisfies EXACTLY, for every m >= 1,
          lambda * (1 - r)/(1 + r)  <  rho_m := lambda lammax(KMS_m(r))
                                    <  lambda * (1 + r)/(1 - r),
      both bounds m-independent (Gershgorin rows of KMS_m are
      1 + 2(r + ... ) < (1+r)/(1-r), strict for finite m; lower bound by
      the exact tridiagonal inverse, whose rows sum to >= (1-r)/(1+r)
      times the diagonal scale). Consequently the carrier's normalised
      A-line gap obeys
          gap_m  >=  -log lambda - log((1+r)/(1-r))   for ALL m,
      and this is certified positive (ratio < 1) on a kappa grid. Per-m
      exact brackets of lammax(KMS_m(r)) for m = 2..8 are also pinned
      (interval LDL inertia bisection, YM-7 pattern) and are seen to
      increase monotonically toward the uniform ceiling — the three
      spot checks of the working session are now a theorem plus a
      monotone ladder.

 (T3) FAIL-CLOSED EDGE. At kappa = 2 the uniform ceiling exceeds 1 and
      the capsule REFUSES the gap claim (the A-line would be degenerate
      with, or above, the vacuum within the carrier); the refusal is
      recorded, not hidden.

HONEST REMAINDER (unchanged in kind, narrowed in place):
  * This is the A-LINE CARRIER ONLY. The off-carrier complement of
    L^2(SU(2)^{2m})^Ad — B-contents, mixed channels (1/2,1/2)^k, all
    higher spins — GROWS with m and is NOT controlled here. YM-14's
    Haynsworth dock handled it for m = 2; extending the complement bound
    uniformly in m is the next dock (YM-16), named below. Until it is
    done, T2 is a statement about the compressed operator, not a
    certified gap of S_m.
  * The chain is still the toy carrier (2 holonomies per block, fixed
    beta = 2, Wilson T0). 2D lattice, AF trajectory, tightness, OS,
    non-triviality, metric universality, Clay: all OPEN.
  * NAMED NEXT (YM-16): complement top of S_m in the chain grammar —
    the complement's T0 top is lambda^2 for every m (two half-spins is
    the cheapest content off the carrier, m-independent), so the
    m-dependence sits entirely in |PMQ|; the doubling identity
    m_kappa^2 = m_{2kappa} survives product-of-bridges verbatim, which
    is why this dock is finite.

Controls:
  C1  symbolic engine == closed form, all entries, m = 2..8 (exact
      polynomial identity, zero tolerance).
  C2  pairing-normalisation tamper (1/d_j -> 1 in the convolution)
      breaks the identity already at m = 2 — the check bites.
  C3  free limit: kappa = 0 gives r = 0 and KMS_m = identity (exact).
  C4  exact tridiagonal inverse of KMS_m verified (interval product
      within enclosure) for m = 2..8.
  C5  per-m lammax brackets strictly increasing in m and strictly below
      the uniform ceiling at every grid kappa.
  C6  fail-closed at kappa = 2 (ceiling >= 1 -> claim refused).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import (  # noqa: E402
    Iv, bessel_I, log_iv, _dec, canonical_sha, TERMS, LOG_TERMS,
)
from ym4_symmetry_protected import chi_mul  # noqa: E402
from ym6_seam_integer_dock import ldl_inertia, _r  # noqa: E402

BETA = F(2)
GRID = [F(1, 8), F(1, 4), F(1, 2), F(1), F(2)]     # kappa = 2 must refuse
M_RANGE = list(range(2, 9))                         # symbolic / per-m checks
P_SYM = 4                                           # formal f_0..f_{P_SYM}
BISECT_TOL = F(1, 10 ** 12)

# ------------------------------------------------ formal polynomial ring
# monomial = tuple of exponents over (f_0, ..., f_{P_SYM}); poly = dict


def p_const(c):
    return {tuple([0] * (P_SYM + 1)): F(c)} if c != 0 else {}


def p_add(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def p_scale(a, c):
    return {k: v * c for k, v in a.items()} if c != 0 else {}


def p_mul_f(a, t):
    """multiply polynomial by the indeterminate f_t."""
    out = {}
    for k, v in a.items():
        kk = list(k)
        kk[t] += 1
        out[tuple(kk)] = out.get(tuple(kk), F(0)) + v
    return out


def p_mul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k = tuple(x + y for x, y in zip(ka, kb))
            out[k] = out.get(k, F(0)) + va * vb
    return {k: v for k, v in out.items() if v != 0}


def p_pow_f(t, e):
    out = p_const(1)
    for _ in range(e):
        out = p_mul_f(out, t)
    return out


# ---------------------------------------- symbolic chain-integration engine
def chain_entry(m: int, ins: dict, tamper=False):
    """Exact formal value of Int prod_i chi_{ins[i]}(A_i) * prod_bridges
    sum_j d_j f_j chi_j(A_i A_{i+1}^{-1})  over SU(2)^m (Haar), as a
    polynomial in the indeterminates f_j.  `ins` maps site -> twice-spin
    of the inserted character (missing = chi_0).

    State after integrating sites 1..i-1: dict content_at_site_i -> poly.
    Transfer across bridge i: chi_c(A_i) * sum_j d_j f_j chi_j(A_i A^{-1})
    integrates to  d_c f_c chi_c(A_{i+1}) / d_c  = f_c chi_c(A_{i+1})
    (tamper: keeps the d_c factor -> d_c f_c).
    """
    state = {0: p_const(1)}
    for i in range(1, m + 1):
        # multiply in the insertion at site i
        t = ins.get(i, 0)
        nstate = {}
        for c, poly in state.items():
            for c2 in chi_mul(c, t):
                nstate[c2] = p_add(nstate.get(c2, {}), poly)
        state = nstate
        if i == m:
            # last site: Int chi_c dA = delta_{c0}
            return state.get(0, {})
        # bridge i -> i+1
        nstate = {}
        for c, poly in state.items():
            if c > P_SYM:
                raise ValueError("raise P_SYM")
            q = p_mul_f(poly, c)
            if tamper:
                q = p_scale(q, c + 1)          # keeps d_c: wrong pairing
            nstate[c] = p_add(nstate.get(c, {}), q)
        state = nstate
    return state.get(0, {})


def closed_form_entry(m: int, p: int, q: int):
    """Formal closed form: p = q = 0 vacuum -> f_0^{m-1};
    p,q in 1..m -> f_0^{m-1-|p-q|} f_{1/2}^{|p-q|}; mixed -> 0."""
    if p == 0 and q == 0:
        return p_pow_f(0, m - 1)
    if p == 0 or q == 0:
        return {}
    d = abs(p - q)
    return p_mul(p_pow_f(0, m - 1 - d), p_pow_f(1, d))


def t1_symbolic_check(m: int, tamper=False) -> bool:
    for p in range(0, m + 1):
        for q in range(p, m + 1):
            ins = {}
            if p:
                ins[p] = 1
            if q:
                ins[q] = ins.get(q, 0)
            if p and q and p == q:
                ins = {p: 1}
                # two insertions at the same site: chi12*chi12 = chi0+chi1
                eng = {}
                for c2 in chi_mul(1, 1):
                    eng = p_add(eng, chain_entry_with_site_content(m, p, c2,
                                                                   tamper))
            else:
                if q:
                    ins[q] = 1
                eng = chain_entry(m, ins, tamper)
            if eng != closed_form_entry(m, p, q):
                return False
    return True


def chain_entry_with_site_content(m, site, content, tamper=False):
    return chain_entry(m, {site: content}, tamper) if content else \
        chain_entry(m, {}, tamper)


# NOTE on the diagonal: chi12(A_p)^2 = chi_0(A_p) + chi_1(A_p); the
# chi_1 branch is propagated by the engine with f_1 factors and killed
# at the open end — the closed form f_0^{m-1} must come out with NO f_1
# monomial surviving.  That cancellation is the content of "no
# truncation remainder on this carrier".


# ------------------------------------------------------- interval pieces
def lam_half() -> Iv:
    return bessel_I(2, BETA, TERMS) / bessel_I(1, BETA, TERMS)


def r_of(kappa: F) -> Iv:
    if kappa == 0:
        return Iv(F(0))
    return bessel_I(2, kappa, TERMS) / bessel_I(1, kappa, TERMS)


def kms(m: int, r: Iv):
    out = [[None] * m for _ in range(m)]
    for i in range(m):
        for j in range(m):
            d = abs(i - j)
            v = Iv(F(1))
            for _ in range(d):
                v = v * r
            out[i][j] = _r(v)
    return out


def kms_inverse(m: int, r: Iv):
    """(1 - r^2) KMS^{-1} = tridiag(-r ; 1, 1+r^2, ..., 1+r^2, 1)."""
    rr = r * r
    z = Iv(F(0))
    out = [[z] * m for _ in range(m)]
    for i in range(m):
        out[i][i] = Iv(F(1)) if i in (0, m - 1) else Iv(F(1)) + rr
        if i + 1 < m:
            out[i][i + 1] = out[i + 1][i] = -r
    if m == 1:
        out[0][0] = Iv(F(1)) - rr
    return out


def mat_mul(A, B):
    n = len(A)
    return [[_r(sum((A[i][k] * B[k][j] for k in range(n)), Iv(F(0))))
             for j in range(n)] for i in range(n)]


def count_above(M, mu: F):
    """certified #eigenvalues > mu, or None on refusal."""
    n = len(M)
    S = [[_r(M[i][j] - Iv(mu if i == j else F(0))) for j in range(n)]
         for i in range(n)]
    res = ldl_inertia(S)
    return None if res is None else res[0]


def lammax_bracket(m: int, r: Iv):
    """exact-LDL inertia bisection for the top eigenvalue of KMS_m(r)."""
    M = kms(m, r)
    lo, hi = F(1), ((Iv(F(1)) + r) / (Iv(F(1)) - r)).hi   # ceiling
    # invariant: count_above(lo) >= 1, count_above(hi) == 0
    assert count_above(M, hi) == 0
    while hi - lo > BISECT_TOL:
        mid = (lo + hi) / 2
        c = count_above(M, mid)
        if c is None:
            # pivot straddle: shrink from the side that stays certified
            lo2 = lo + (mid - lo) / 3
            c2 = count_above(M, lo2)
            if c2 is None:
                break
            if c2 >= 1:
                lo = lo2
            else:
                hi = lo2
            continue
        if c >= 1:
            lo = mid
        else:
            hi = mid
    return Iv(lo, hi)


def uniform_bounds(kappa: F):
    r = r_of(kappa)
    lam = lam_half()
    one = Iv(F(1))
    ceil = lam * (one + r) / (one - r)
    floor = lam * (one - r) / (one + r)
    return r, lam, _r(floor), _r(ceil)


# ------------------------------------------------------------------ run
def run():
    # ---- T1 / C1 / C2
    c1 = all(t1_symbolic_check(m) for m in M_RANGE)
    c2 = not t1_symbolic_check(2, tamper=True)
    # diagonal has no f_1 monomial (no-remainder statement)
    diag = {}
    for c2_ in chi_mul(1, 1):
        diag = p_add(diag, chain_entry_with_site_content(4, 2, c2_))
    no_f1 = all(k[2] == 0 for k in diag)

    # ---- T2 per kappa
    lam = lam_half()
    rows = {}
    all_cert = True
    refused_at_2 = False
    c5 = True
    c4 = True
    for kap in GRID:
        r, _, floor, ceil = uniform_bounds(kap)
        row = {"r_lo": _dec(r.lo, 25), "r_hi": _dec(r.hi, 25),
               "rho_uniform_floor_lo": _dec(floor.lo, 20),
               "rho_uniform_ceiling_hi": _dec(ceil.hi, 20)}
        if ceil.hi < 1:
            g = -log_iv(Iv(ceil.hi), LOG_TERMS)
            row["status"] = "CERTIFIED_UNIFORM_IN_m"
            row["gap_lower_bound_all_m"] = _dec(g.lo, 20)
        else:
            row["status"] = "REFUSED_CEILING_NOT_BELOW_1"
            if kap == F(2):
                refused_at_2 = True
            elif kap < F(2):
                all_cert = False
        per_m = {}
        prev = None
        for m in M_RANGE:
            b = lammax_bracket(m, r)
            per_m[str(m)] = [_dec(b.lo, 12), _dec(b.hi, 12)]
            # C5: strictly increasing, strictly below the ceiling
            if prev is not None and not (b.lo > prev.hi):
                c5 = False
            if not ((lam * b).hi < ceil.hi + F(1, 10 ** 9)):
                c5 = False
            if not (b.hi < ((Iv(F(1)) + r) / (Iv(F(1)) - r)).hi):
                c5 = False
            prev = b
            # C4 inverse identity
            P = mat_mul(kms(m, r), kms_inverse(m, r))
            scale = Iv(F(1)) - r * r
            for i in range(m):
                for j in range(m):
                    target = scale if i == j else Iv(F(0))
                    if not (P[i][j].lo <= target.hi + F(1, 10 ** 25)
                            and P[i][j].hi >= target.lo - F(1, 10 ** 25)):
                        c4 = False
        row["lammax_KMS_m_brackets"] = per_m
        rows[str(kap)] = row

    # ---- C3 free limit
    r0 = r_of(F(0))
    K0 = kms(4, r0)
    c3 = all(K0[i][j].lo == (1 if i == j else 0) and
             K0[i][j].hi == (1 if i == j else 0) for i in range(4)
             for j in range(4))

    ok = c1 and c2 and no_f1 and c3 and c4 and c5 and all_cert and refused_at_2
    cert = {
        "certificate_type": "YM15_CHAIN_CLOSED_FORM_BLOCK_TRANSFER_DOCK",
        "claim_status": "first_m_uniform_statement_in_overlap_regime_"
                        "A_line_carrier_only",
        "theorems": {
            "T1_closed_form_all_m":
                "compressed bridge product on W_{m+1} = diag(f0^{m-1}, "
                "f0^{m-1} KMS_m(r)), r = f_{1/2}/f0 = I2(kappa)/I1(kappa); "
                "proof = SU(2) convolution lemma by induction along the "
                "chain; machine-checked as an exact polynomial identity "
                "in formal f_j for m = 2..8; NO truncation remainder on "
                "this carrier",
            "T2_m_uniform_ratio_bracket":
                "lambda(1-r)/(1+r) < rho_m < lambda(1+r)/(1-r) for every m; "
                "normalised A-line gap >= -log lambda - log((1+r)/(1-r)) "
                "uniformly in m; certified positive on the grid kappa in "
                "{1/8,1/4,1/2,1}; per-m exact brackets pinned and monotone",
            "T3_fail_closed_edge":
                "kappa = 2: ceiling >= 1, claim refused",
        },
        "convention": {
            "normalisation": "rho = (A-line top)/(vacuum) of the COMPRESSED "
                             "transfer; gap = -log rho (YM-5/6 convention)",
            "carrier": "W_{m+1} = {1, chi12(A_1..A_m)}, orthonormal, exact "
                       "T0 eigenvectors; B-contents inert by parity",
            "T0": "Wilson beta = 2, lambda = I2(2)/I1(2) (YM-14 convention)",
        },
        "grid": rows,
        "honest_remainder": {
            "off_carrier": ("complement of W_{m+1} grows with m and is NOT "
                            "bounded here; T2 is a compressed-operator "
                            "statement, not a certified gap of S_m"),
            "named_next": ("YM-16: m-uniform complement bound — complement "
                           "T0 top is lambda^2 for all m; doubling identity "
                           "survives the bridge product verbatim"),
            "scope": "toy carrier; 2D, AF trajectory, tightness, OS, "
                     "non-triviality, metric universality, Clay OPEN",
        },
        "controls": {
            "C1_symbolic_engine_equals_closed_form_m2_to_8": bool(c1),
            "C2_pairing_normalisation_tamper_breaks_identity": bool(c2),
            "C2b_no_f1_monomial_on_diagonal": bool(no_f1),
            "C3_free_limit_identity": bool(c3),
            "C4_tridiagonal_inverse_identity": bool(c4),
            "C5_per_m_monotone_below_ceiling": bool(c5),
            "C6_fail_closed_kappa_2": bool(refused_at_2),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM15_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM15.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"])
    for k, v in cert["controls"].items():
        print(" ", k, v)
    for k, v in cert["grid"].items():
        print(" kappa", k, v["status"], v.get("gap_lower_bound_all_m"),
              "ceiling", v["rho_uniform_ceiling_hi"][:12])
        print("   lammax m=2..8:",
              [x[0][:8] for x in v["lammax_KMS_m_brackets"].values()])
    print("sha256:", sha)
