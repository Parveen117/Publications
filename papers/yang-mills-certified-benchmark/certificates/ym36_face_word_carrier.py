"""YM-36: T54' ON THE FABRIC AT CONTENT 1/2 — the alternating product
T = K^{1/2} M_w K^{1/2} compressed to the FACE-WORD carrier
    W_m = span{ B_S = prod_{l in S} chi_half(H_l) : S subset of the m-1 bridge faces }
(dimension 2^{m-1}), i.e. EVERY superposition of content-1/2 face insertions —
the object YM-28 T5 / theorum/28 item 3 left open ("E4D-C on superpositions").

Engine: YM-33's t-row fabric with per-column face labels on the outer rails
(the word S on rail 0, S' on rail 2), the truncated face weight
w_pt/f_0 = 1 + 2 mu chi_half on the middle rail, and K^{1/2} as the two rung
kernels.  K^{1/2} needs sqrt(lambda_c): we use a DECLARED RATIONAL-SQUARE
KERNEL  lambda_half = (11/16)^2, lambda_1 = (3/8)^2  (close to a = 1:
0.4727 vs 0.4724, 0.1406 vs 0.1353) so that every matrix element is an exact
rational and no surd survives (the engine asserts it).  B_S are orthonormal
(gauge change of variables: faces independent, int chi_half^k = 0,1,0,2 for
k=0..3 — no Gram needed).

 T1  Exact rational 2^{m-1} x 2^{m-1} matrices T_m(mu) for m = 2..6 at
     mu = 1/32, 1/16, 1/8 (mu = r = f_half/f_0; 1/32 ~ kappa 1/8, 1/8 ~ kappa 1/2).
     Symmetric; vacuum column T_{00} = 1 exactly (mean law, YM-28 T1).
 T2  rho_m := (second eigenvalue)/(vacuum eigenvalue) — the excitation/vacuum
     ratio over ALL superpositions.  Successive increments rho_{m+1} - rho_m
     are positive and decay geometrically: every ratio of consecutive
     increments lies in [2/5, 3/5] at all three mu (certified), so under the
     DECLARED continuation q <= 2/3 the limit is bracketed
         rho_6 <= rho_inf <= rho_6 + (rho_6 - rho_5) * q / (1 - q).
     The m-uniform limit sits at lambda^2 (1 + O(mu^2)) for small mu and
     BELOW lambda^2 for mu = 1/8 (the vacuum grows faster than the top).
 T3  Superposition test (the E4D-C question at this content): the top
     eigenvector's overlap with the single-insertion sector is computed; the
     top is NOT a single insertion (single-insertion weight 0.98 / 0.92 / 0.73
     at m = 6), yet rho_inf exceeds the single-insertion ratio (m = 2, one
     face) by at most 0.24 % (mu = 1/32) / 0.93 % (1/16) / 3.5 % (mu = 1/8)
     — superpositions do not open a gap-closing channel at content 1/2.
 T4  HONEST BOUNDARY.  (i) Compression: by interlacing every compressed
     eigenvalue is a LOWER bound of the true one, so rho_m is a statement
     about the compressed operator, not an upper bound on the true
     excitation; (ii) contents j >= 1 of faces and of fused sites are
     outside W_m (the YM-31 rung truncation) — their first effect is the
     YM-35 adjacent commutator (energy lambda^2 (3/4)(1-sqrt lambda_1)^2,
     relative size ~mu^2); (iii) rational surrogate kernel; (iv) m <= 6 in
     the certificate (m = 7 reproduced once off-line: rho_7 - rho_6 =
     1.5e-5 / 5.8e-5 / 2.1e-4, ratios 0.56, consistent with T2).
     E4D-C remains OPEN; at content 1/2 on the face-word carrier it HOLDS
     in compression with a geometric tail.
"""

from fractions import Fraction as F
import itertools
import json
import os
import sys

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ym1_certified_gap import canonical_sha  # noqa: E402
import ym33_t_row_engine as E  # noqa: E402

mp.mp.dps = 40
half, one, zero = F(1, 2), F(1), F(0)

# per-column outer-rail labels: same engine, rails_u[s] may be a function of the column
_src = open(os.path.join(HERE, "ym33_t_row_engine.py")).read()
_code = "def fabric_partition_cols" + _src.split("def fabric_partition")[1].split("\n_VT = {}")[0]
_code = _code.replace("rails_u[s][an]", "(rails_u[s](i) if callable(rails_u[s]) else rails_u[s]).get(an, F(0))")
_ns = dict(E.__dict__)
exec(_code, _ns)  # noqa: S102 — engine variant, same code path
fabric_partition_cols = _ns["fabric_partition_cols"]

S_HALF, S_ONE = F(11, 16), F(3, 8)          # sqrt(lambda_half), sqrt(lambda_1): declared rational-square kernel
K_HALF = {zero: F(1), half: S_HALF, one: S_ONE}
MMAX = 6


def rail(S):
    return lambda i: ({half: F(1, 2)} if i in S else {zero: F(1)})


def word_matrix(mu: F, m: int):
    faces = list(range(m - 1))
    words = [frozenset(c) for r in range(m) for c in itertools.combinations(faces, r)]
    mid = {zero: F(1), half: mu}
    T = {}
    for S in words:
        for Sp in words:
            if (Sp, S) in T:
                T[(S, Sp)] = T[(Sp, S)]
                continue
            T[(S, Sp)] = fabric_partition_cols([rail(S), mid, rail(Sp)], K_HALF, m, 3, {}, jmax=one)
    return words, T


def spectrum(words, T):
    n = len(words)
    M = mp.matrix(n, n)
    for i, S in enumerate(words):
        for j, Sp in enumerate(words):
            M[i, j] = mp.mpf(T[(S, Sp)].numerator) / T[(S, Sp)].denominator
    ev, Q = mp.eigsy(M)
    order = sorted(range(n), key=lambda i: -ev[i])
    vac, top = ev[order[0]], ev[order[1]]
    v = Q[:, order[1]]
    single = sum(v[i] ** 2 for i, S in enumerate(words) if len(S) == 1)
    return vac, top, single


def run():
    out = {}
    ok = True
    for mu in (F(1, 32), F(1, 16), F(1, 8)):
        rows = {}
        rho = {}
        sym_ok = True
        for m in range(2, MMAX + 1):
            words, T = word_matrix(mu, m)
            if T[(frozenset(), frozenset())] != 1:
                ok = False
            if any(T[(S, Sp)] != T[(Sp, S)] for S in words for Sp in words):
                sym_ok = False
            vac, top, single = spectrum(words, T)
            rho[m] = top / vac
            rows[m] = {"dim": len(words), "vacuum": mp.nstr(vac, 15), "top": mp.nstr(top, 15), "rho": mp.nstr(rho[m], 15),
                       "top_vector_single_insertion_weight": mp.nstr(single, 8)}
        inc = {m: rho[m + 1] - rho[m] for m in range(2, MMAX)}
        pos = all(v > 0 for v in inc.values())
        qs = {m: inc[m + 1] / inc[m] for m in range(2, MMAX - 1)}
        q_ok = all(mp.mpf(2) / 5 <= q <= mp.mpf(3) / 5 for q in qs.values())
        q_decl = mp.mpf(2) / 3
        lo, hi = rho[MMAX], rho[MMAX] + inc[MMAX - 1] * q_decl / (1 - q_decl)
        lam2 = (S_HALF ** 2) ** 2
        excess = hi / rho[2] - 1
        out[str(mu)] = {"rows": {str(k): v for k, v in rows.items()},
                        "increments": {str(k): mp.nstr(v, 8) for k, v in inc.items()},
                        "increment_ratios": {str(k): mp.nstr(v, 6) for k, v in qs.items()},
                        "rho_inf_bracket_under_q_le_2_3": [mp.nstr(lo, 12), mp.nstr(hi, 12)],
                        "lambda_squared_reference": str(lam2),
                        "excess_over_single_insertion_ratio_rho2": mp.nstr(excess, 6),
                        "symmetric": sym_ok, "increments_positive": pos, "increment_ratios_in_2_5_to_3_5": q_ok}
        ok = ok and sym_ok and pos and q_ok and excess < mp.mpf(1) / 25
    return {
        "certificate_type": "YM36_T54PRIME_ON_FABRIC_CONTENT_HALF_FACE_WORD_CARRIER",
        "declared": {"rational_square_kernel": {"sqrt_lambda_half": str(S_HALF), "sqrt_lambda_1": str(S_ONE)},
                     "face_weight": "w_pt/f0 = 1 + 2 mu chi_half (content-1/2 truncation, YM-31)",
                     "carrier": "all 2^(m-1) face-insertion words, orthonormal", "tail_continuation": "q <= 2/3 (declared, observed in [2/5,3/5])"},
        "grid": out,
        "honest_boundary": ["compression (interlacing): compressed eigenvalues are lower bounds of the true ones; rho_m is a compressed-operator statement",
                            "contents j>=1 outside the carrier; first effect = YM-35 adjacent commutator, relative size ~mu^2",
                            "m <= 6 certified; m = 7 reproduced off-line once (increment ratios 0.56)",
                            "E4D-C OPEN; holds at content 1/2 on the face-word carrier in compression with geometric tail"],
        "verdict": "PASS" if ok else "FAIL",
    }


if __name__ == "__main__":
    cert = run()
    json.dump(cert, open(os.path.join(HERE, "YM36_RESULT.json"), "w"), indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    open(os.path.join(HERE, "EXPECTED_YM36.sha256"), "w").write(sha + "\n")
    print(cert["verdict"])
    for mu, g in cert["grid"].items():
        print(mu, g["rho_inf_bracket_under_q_le_2_3"], g["increment_ratios"], "excess", g["excess_over_single_insertion_ratio_rho2"], g["rows"]["6"]["top_vector_single_insertion_weight"])
    print("sha256:", sha)
