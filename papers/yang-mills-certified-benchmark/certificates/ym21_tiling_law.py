"""YM-21: THE TILING LAW — what Recognition-Stokes actually gives the
gap, derived inside the fabric (YM-F1), with one ill-posed question of
my own killed by an exact witness first.

STANDING CORRECTION (YM-F1 "native form of YM-21"). I posed: does the
boundary-holonomy residue of a fabric with exact interior cancellation
control the face-residue sum uniformly in the number of faces? NO, by
type and by witness: faces H_1 = g, H_2 = g^{-1} have boundary I
(residue 0) and face-residue sum 2 rho(g), arbitrary in [0,4]. Boundary
closure never controls interior residue — that is the fabric appendix's
own point (visible return != recognition return), and I wrote the
question against it. Certified here as a control, then replaced.

THE RIGHT NATIVE OBJECT IS THE TILING, NOT THE BOUNDARY.
Stokes telescoping (YM-F1 T3) says a rung-to-rung transport A_p^{-1}A_q
is the product of exactly the faces between p and q. Under the
face-weighted measure each face contributes its coefficient ladder
f_j/f_0 =: r_j (r_{1/2} = r). So a correlation between two rungs is, at
leading order, the product of r over the faces a minimal tiling must
cross. This is the native reading of YM-15's closed form.

CERTIFIED (exact rational / interval; no floats in verdicts):

 (T1) NEGATIVE WITNESS (the correction): for every rational unit g,
      the two-face fabric (g, g^{-1}) has boundary residue exactly 0
      and face-residue sum exactly 2 rho(g) > 0 for g != I. Boundary
      residue does not bound face residue: no constant works.

 (T2) STOKES SUB-TELESCOPING => SPATIAL TILING LAW, EXACT AND m-UNIFORM.
      For 1 <= p < q <= m the product of the faces l_p ... l_{q-1}
      equals A_p^{-1} A_q exactly (sub-telescoping, certified m=2..8),
      and YM-15's carrier entry <chi12(A_p), m_kappa chi12(A_q)> /
      f_0^{m-1} equals r^{q-p} = product over exactly those faces of
      the face coefficient r (pinned YM-15 closed form, re-read here as
      a face product). Spatial decay along the chain is therefore
      EXACTLY geometric in the number of faces crossed, with rate
      -Log_Sigma(r) = A_Sigma-free native quantity, independent of m.
      This is a theorem about the compressed bridge product and holds
      for every m (YM-15 T1); it is the 1D "area law" of the fabric.

 (T3) TIME DECAY AT LEADING TILING ORDER IS m-UNIFORM AND EQUALS
      YM-15'S BRACKET. Build the 2D fabric m x t (t time steps). A
      time face carries the free coefficient lambda (= lambda_{1/2} of
      T0), a space face carries r. The leading-order two-point function
      of chi12 at rungs p, q separated by t time steps is the tiling
      sum over monotone face paths = [(lambda KMS_m(r))^t]_{pq}, i.e.
      the t-th power of YM-15's compressed A-block (vacuum-normalised).
      Its growth is bracketed for EVERY m by
          lambda^t ((1-r)/(1+r))^t  <=  row growth  <=  lambda^t ((1+r)/(1-r))^t
      (KMS row sums), so the leading-order time decay rate is at least
      -Log(lambda) - A_Sigma(r) uniformly in m — exactly YM-15 T2 /
      YM-20 T2, now DERIVED as a tiling count rather than observed as
      a matrix bound. Certified: exact powers t = 1..6 for m = 2..8
      lie inside the bracket; the bracket is m-independent.

 (T4) WHAT REMAINS, NAMED IN FABRIC LANGUAGE. The full two-point
      function adds tilings with higher face content (coefficients r_j,
      j >= 1, geometric in kappa: r_j = f_j/f_0 certified decreasing)
      and branching tilings (faces shared by several paths). The
      expansion parameter is the coefficient ladder r_j; its
      convergence uniformly in m is the native cluster statement
      (YM-18 T3 design, now with polymers = connected face sets of the
      fabric and activities = face coefficients). NOT derived here.
      Also certified: the ladder r_1/r_{1/2} and r_{3/2}/r_1 are < 1
      at every grid kappa, so higher contents are subleading face by
      face — necessary for, not sufficient for, convergence.

HONEST STATUS. Spatial decay: exact, native, m-uniform (T2). Time decay:
leading order native and m-uniform (T3); full order OPEN (T4). No
volume-uniform gap claimed. Cutoff a -> 0: untouched.

Controls:
  C1  negative witness: boundary residue 0, face sum 2 rho(g) for 20 g.
  C2  sub-telescoping exact for all 1<=p<q<=m, m=2..8.
  C3  YM-15 entry == r^{#faces between} (pinned closed form).
  C4  tiling powers inside the m-independent bracket, t=1..6.
  C5  coefficient ladder strictly decreasing at every grid kappa.
  C6  tamper: a fabric with one unpaired interior edge breaks C2.
"""

from fractions import Fraction as F
import json
import os
import random
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, log_iv, _dec, canonical_sha, LOG_TERMS  # noqa: E402
from ym4_symmetry_protected import face_coeffs  # noqa: E402
from ym6_seam_integer_dock import _r  # noqa: E402
from ym15_chain_closed_form import (  # noqa: E402
    closed_form_entry, p_pow_f, p_mul, kms, r_of, lam_half,
)
from ymf1_chain_fabric import (  # noqa: E402
    rational_unit, ladder_faces, QONE, Quat,
)

random.seed(20260822)
GRID = [F(1, 8), F(1, 4), F(1, 2), F(1)]
M_RANGE = list(range(2, 9))
T_MAX = 6


def rnd_unit():
    return rational_unit(F(random.randint(-9, 9), random.randint(1, 7)),
                         F(random.randint(-9, 9), random.randint(1, 7)),
                         F(random.randint(-9, 9), random.randint(1, 7)))


def mat_pow_rows(M, t):
    n = len(M)
    P = [[Iv(F(1) if i == j else F(0)) for j in range(n)] for i in range(n)]
    for _ in range(t):
        P = [[_r(sum((P[i][k] * M[k][j] for k in range(n)), Iv(F(0))))
              for j in range(n)] for i in range(n)]
    return P


def run():
    random.seed(20260822)   # reseed: pins must not depend on test order
    # ---- T1 / C1 negative witness
    c1 = True
    for _ in range(20):
        g = rnd_unit()
        faces = [g, g.inv()]
        boundary = faces[0] * faces[1]
        if boundary != QONE or boundary.residue() != 0:
            c1 = False
        if g != QONE and not (faces[0].residue() + faces[1].residue() > 0):
            c1 = False
        if faces[0].residue() + faces[1].residue() != 2 * g.residue():
            c1 = False

    # ---- T2 / C2 / C3 / C6 sub-telescoping and face product
    c2 = c3 = True
    c6 = False
    for m in M_RANGE:
        rungs = [rnd_unit() for _ in range(m)]
        faces = ladder_faces(rungs)
        for p in range(1, m + 1):
            for q in range(p + 1, m + 1):
                prod = QONE
                for i in range(p - 1, q - 1):
                    prod = prod * faces[i]
                if prod != rungs[p - 1].inv() * rungs[q - 1]:
                    c2 = False
                # closed-form entry is f0^{m-1-(q-p)} f_{1/2}^{q-p}
                want = p_mul(p_pow_f(0, m - 1 - (q - p)), p_pow_f(1, q - p))
                if closed_form_entry(m, p, q) != want:
                    c3 = False
        if m >= 4:
            bad = list(faces)
            bad[1] = rungs[1].inv() * rungs[2].inv()     # unpaired edge
            prod = QONE
            for i in range(0, 3):
                prod = prod * bad[i]
            if prod != rungs[0].inv() * rungs[3]:
                c6 = True

    # ---- T3 / C4 tiling powers inside the m-independent bracket
    lam = lam_half()
    c4 = True
    rows = {}
    for kap in GRID:
        r = r_of(kap)
        one = Iv(F(1))
        up = lam * (one + r) / (one - r)
        lo = lam * (one - r) / (one + r)
        per_m = {}
        for m in M_RANGE:
            M = [[_r(lam * x) for x in row] for row in kms(m, r)]
            for t in range(1, T_MAX + 1):
                P = mat_pow_rows(M, t)
                for i in range(m):
                    rs = sum((P[i][j] for j in range(m)), Iv(F(0)))
                    up_t = one
                    lo_t = one
                    for _ in range(t):
                        up_t = up_t * up
                        lo_t = lo_t * lo
                    if not (rs.hi <= up_t.hi + F(1, 10 ** 20)
                            and rs.lo >= lo_t.lo - F(1, 10 ** 20)):
                        c4 = False
            per_m[str(m)] = "inside"
        rate = (-log_iv(lam, LOG_TERMS)) - log_iv((one + r) / (one - r),
                                                   LOG_TERMS)
        rows[str(kap)] = {
            "r": _dec(r.lo, 15),
            "leading_order_time_decay_rate_uniform_in_m": _dec(rate.lo, 15),
            "spatial_rate_-log_r": _dec((-log_iv(r, LOG_TERMS)).lo, 15),
            "tiling_powers_checked_m": list(per_m),
        }
    # ---- T4 / C5 coefficient ladder
    c5 = True
    ladder = {}
    for kap in GRID:
        f = face_coeffs(kap, 6)
        r1 = f[2] / f[0]
        r12 = f[1] / f[0]
        r32 = f[3] / f[0]
        if not (r1.hi < r12.lo and r32.hi < r1.lo):
            c5 = False
        ladder[str(kap)] = [_dec(r12.lo, 10), _dec(r1.lo, 10),
                            _dec(r32.lo, 10)]

    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM21_TILING_LAW_NATIVE_DECAY_ON_THE_FABRIC",
        "claim_status": "spatial_decay_exact_m_uniform__time_decay_leading_"
                        "order_m_uniform__full_order_OPEN",
        "standing_correction": {
            "target": "YM-F1 'native form of YM-21'",
            "content": "boundary residue cannot control face residue "
                       "(witness g, g^{-1}); replaced by the tiling law",
        },
        "theorems": {
            "T1_negative_witness": bool(c1),
            "T2_spatial_tiling_law":
                "sub-telescoping exact; YM-15 entry = r^{faces crossed}; "
                "spatial decay exactly geometric, rate -Log r, all m",
            "T3_time_decay_leading_order":
                "leading tiling two-point = (lambda KMS_m(r))^t; bracketed "
                "by lambda^t((1-/+r)/(1+/-r))^t for every m => rate >= "
                "-Log(lambda) - A_Sigma(r) uniformly — YM-15/20 bound DERIVED "
                "as a tiling count",
            "T4_named_remainder":
                "higher-content and branching tilings; expansion parameter "
                "= face-coefficient ladder r_j (certified decreasing); "
                "uniform convergence = native cluster statement, OPEN",
        },
        "grid": rows,
        "coefficient_ladder_r_half_r_1_r_3half": ladder,
        "controls": {
            "C1_negative_witness": bool(c1),
            "C2_sub_telescoping_exact": bool(c2),
            "C3_entry_is_face_product": bool(c3),
            "C4_tiling_powers_in_bracket": bool(c4),
            "C5_ladder_decreasing": bool(c5),
            "C6_unpaired_edge_breaks": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM21_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM21.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for k, v in cert["grid"].items():
        print(" kappa", k, v)
    print("ladder:", cert["coefficient_ladder_r_half_r_1_r_3half"])
    print("sha256:", sha)
