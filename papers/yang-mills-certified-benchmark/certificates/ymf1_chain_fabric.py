"""YM-F1: THE CHAIN AS A RECOGNITION FABRIC — the native carrier for the
Yang-Mills program. Source: vault APPENDIX_DISCRETE_HOLONOMY_FABRICS_AND_
RECOGNITION_STOKES (oriented complex (V,E,F), edge transports T_e, face
holonomies H_l, Recognition-Stokes closure) + EMK-1 determinant seam
ladder + F00-E native exponential. After this capsule the carrier layer
of YM-1..18 is NATIVE in the framework's own vocabulary; not one number
of the program changes.

THE FOUR MOVES (all exact, rational; no floats in verdicts):

 (T1) THE GAUGE GROUP IS THE EMK BLOCK OVER THE CUT-COMPLEX FIELD.
      EMK-1's primitive block M = aI + bK + cR + dRK (K^2 = I, R^2 = -I,
      RK = -KR) has real matrix representative [[a-d, b-c],[b+c, a+d]]
      and determinant (a^2 - b^2) + (c^2 - d^2) = Delta_par + Delta_perp.
      Give the SEAM-ODD coefficients b, d the cut-iota weight
      iota_Sigma (F00-E: iota^2 = -1). Then the three elements
          R,   iota K,   iota RK
      pairwise anticommute and each squares to -I, so
          span_{R_Sigma}{ I, R, iota K, iota RK }
      is Hamilton's quaternion algebra, and the determinant identity
      becomes  det = (a^2 + b^2) + (c^2 + d^2) = N_Sigma(M):
      the seam channel and the rotational channel are the two halves of
      the quaternion norm. SU(2) = { det = 1 } = { Delta_par + Delta_perp = 1 }
      is the unit sphere of the EMK block with iota-twisted seam sectors.
      Certified: (a) the three anticommutation/square relations as exact
      2x2 identities over the Gaussian rationals Q(iota); (b) quaternion
      product == matrix product (map is a homomorphism) on a rational
      grid; (c) det == N_Sigma; (d) M^dagger M = I exactly on rational
      unit quaternions (stereographic Pythagorean parametrisation — no
      trigonometry, no float). The SU(2) character chi_{1/2}(U) = Tr U
      = 2a is twice the IDENTITY-sector coefficient: the Wilson face
      reads the I-sector of the EMK block.

 (T2) THE CHAIN IS A LADDER FABRIC; BRIDGE FACES ARE ITS PLAQUETTES.
      Fabric: vertices b_1..b_m (bottom), t_1..t_m (top); rails
      b_i -> b_{i+1}, t_i -> t_{i+1} carry the identity transport (the
      gauge-fixed tree of YM-11 gate 1); rungs b_i -> t_i carry A_i.
      Faces l_i = (b_i, b_{i+1}, t_{i+1}, t_i). Exact:
          H_{l_i} = A_i^{-1} A_{i+1},   Tr H_{l_i} = Tr(A_i A_{i+1}^{-1}),
      so the bridge weight exp[(kappa/2) Tr(A_i A_{i+1}^{-1})] IS the
      face weight of the ladder fabric. Certified on random rational
      SU(2) data for m = 2..8.

 (T3) RECOGNITION-STOKES CLOSURE IS EXACT ON THE CHAIN. The product of
      face holonomies (based along the bottom rail) telescopes:
          H_{l_1} H_{l_2} ... H_{l_{m-1}} = A_1^{-1} A_m = H_{boundary},
      every interior rung A_i (1 < i < m) appearing once as A_i and once
      as A_i^{-1} and cancelling EXACTLY — the appendix's interior
      residue E_int is identically zero on the chain, and the fabric's
      boundary holonomy is the end-to-end rung ratio. Certified for
      m = 2..8; CONTROL: reversing one interior rung's orientation
      (an unpaired interior edge) breaks the identity.

 (T4) THE ACTION IS A SUM OF FACE RESIDUES; THE WEIGHT IS A NATIVE
      EXPONENTIAL. With lawful face holonomy I, define the face residue
          rho_l = 1 - a(H_l) = 1 - (1/2) Tr H_l  in [0, 2]  (exact, rational),
      the Cayley-free seam distance of the plaquette from lawful closure
      (rho = 0 iff H = I; rho = 2 iff H = -I, the center). Then
          prod_l exp[(kappa/2) Tr H_l] = Exp_Sigma(kappa (m-1)) *
                                         Exp_Sigma(- kappa sum_l rho_l),
      i.e. the chain weight is the native exponential of minus the
      fabric's total face residue (the Wilson action as a residue sum),
      times the constant Exp_Sigma(kappa) per face. YM-16's "sup price
      e^kappa per bridge" is exactly the rho = 0 face. Certified
      two-route (direct vs residue form) with enclosures.

WHAT THIS CHANGES. YM-1..18 used "L^2(SU(2)^m), Haar, characters". In
fabric language: the state space is functions of rung transports in the
EMK block over C_Sigma; Haar integration is the stationary-path
extraction Phi_Sigma of RH T01 (the pairing delta_{a=b=c}/d_a that the
program has used since YM-3 is its value on characters — that
identification is the ONE remaining shadow, recorded); the transfer
matrix is the fabric's face-weighted transport; seam counts (YM-6/7)
are fabric residues counted above threshold (gold/05's spectral-flow
law typed for the chain). The volume question becomes: does the
boundary-holonomy residue of a fabric of m-1 faces with EXACT interior
cancellation (T3) control the face-residue sum uniformly in m? That is
the native form of YM-21, now posed on a native object.

NOT CLAIMED: any gap; any decay; the Haar <-> Phi_Sigma identification
(declared, shadow); 2D fabrics; continuum.

Controls:
  C1  anticommutation/square relations exact over Q(iota).
  C2  quaternion-vs-matrix product agreement on rational grid.
  C3  det == Delta_par + Delta_perp == N_Sigma; unit sphere = SU(2).
  C4  Stokes telescoping exact m = 2..8; orientation tamper breaks it.
  C5  face residue rho in [0,2] with rho = 0 iff H = I exactly.
  C6  two-route weight enclosures overlap.
"""

from fractions import Fraction as F
import json
import os
import random
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, _dec, canonical_sha  # noqa: E402
from ym2_theta_interacting_gap import exp_point  # noqa: E402

random.seed(20260822)
M_RANGE = list(range(2, 9))
KAPPA = F(1, 4)


# --------------------------------------------------- Gaussian rationals
class GQ:
    """x + iota*y with x, y in Q ; iota^2 = -1 (the cut-iota element)."""
    __slots__ = ("x", "y")

    def __init__(self, x, y=0):
        self.x, self.y = F(x), F(y)

    def __add__(self, o):
        return GQ(self.x + o.x, self.y + o.y)

    def __sub__(self, o):
        return GQ(self.x - o.x, self.y - o.y)

    def __mul__(self, o):
        return GQ(self.x * o.x - self.y * o.y, self.x * o.y + self.y * o.x)

    def __neg__(self):
        return GQ(-self.x, -self.y)

    def conj(self):
        return GQ(self.x, -self.y)

    def __eq__(self, o):
        return self.x == o.x and self.y == o.y


ZERO, ONE, IOTA = GQ(0), GQ(1), GQ(0, 1)


def m2(a, b, c, d):
    return [[a, b], [c, d]]


def mmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0],
             A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0],
             A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(2)] for i in range(2)]


def mscale(A, s: GQ):
    return [[A[i][j] * s for j in range(2)] for i in range(2)]


def meq(A, B):
    return all(A[i][j] == B[i][j] for i in range(2) for j in range(2))


def mdagger(A):
    return [[A[0][0].conj(), A[1][0].conj()], [A[0][1].conj(), A[1][1].conj()]]


def mdet(A):
    return A[0][0] * A[1][1] - A[0][1] * A[1][0]


# EMK-1 primitive generators in the real 2x2 representation
I2 = m2(ONE, ZERO, ZERO, ONE)
K_ = m2(ZERO, ONE, ONE, ZERO)           # seam reflection, K^2 = I
R_ = m2(ZERO, -ONE, ONE, ZERO)          # quarter turn, R^2 = -I
RK = mmul(R_, K_)
NEG_I2 = mscale(I2, -ONE)


def emk_block(a, b, c, d, twist=True):
    """M = aI + bK + cR + dRK with seam-odd coefficients b, d carrying
    iota (twist=True). Returns the 2x2 matrix over Q(iota)."""
    sb = GQ(b) * IOTA if twist else GQ(b)
    sd = GQ(d) * IOTA if twist else GQ(d)
    return madd(madd(mscale(I2, GQ(a)), mscale(K_, sb)),
                madd(mscale(R_, GQ(c)), mscale(RK, sd)))


# ------------------------------------------------ rational quaternions
class Quat:
    __slots__ = ("a", "b", "c", "d")       # a + b i + c j + d k

    def __init__(self, a, b, c, d):
        self.a, self.b, self.c, self.d = F(a), F(b), F(c), F(d)

    def __mul__(self, o):
        a1, b1, c1, d1 = self.a, self.b, self.c, self.d
        a2, b2, c2, d2 = o.a, o.b, o.c, o.d
        return Quat(a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2,
                    a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                    a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2,
                    a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2)

    def inv(self):
        n = self.norm()
        return Quat(self.a / n, -self.b / n, -self.c / n, -self.d / n)

    def norm(self):
        return self.a ** 2 + self.b ** 2 + self.c ** 2 + self.d ** 2

    def __eq__(self, o):
        return (self.a, self.b, self.c, self.d) == (o.a, o.b, o.c, o.d)

    def trace(self):
        return 2 * self.a

    def residue(self):
        """rho = 1 - a = 1 - Tr/2, the face residue from lawful closure."""
        return 1 - self.a


def to_matrix(q: Quat):
    """quaternion a + b i + c j + d k  ->  EMK block with i = R,
    j = iota K, k = iota RK  (T1 identification)."""
    return emk_block(q.a, q.c, q.b, q.d)   # b-slot: j coeff, d-slot: k coeff


def rational_unit(p: F, q: F, r: F) -> Quat:
    """stereographic unit quaternion: no sqrt, no trig."""
    s = p * p + q * q + r * r
    return Quat((1 - s) / (1 + s), 2 * p / (1 + s), 2 * q / (1 + s),
                2 * r / (1 + s))


def random_unit() -> Quat:
    return rational_unit(F(random.randint(-9, 9), random.randint(1, 7)),
                         F(random.randint(-9, 9), random.randint(1, 7)),
                         F(random.randint(-9, 9), random.randint(1, 7)))


QI = Quat(0, 1, 0, 0)
QJ = Quat(0, 0, 1, 0)
QK = Quat(0, 0, 0, 1)
QONE = Quat(1, 0, 0, 0)


# ------------------------------------------------------------- fabric
def ladder_faces(rungs):
    """H_{l_i} = A_i^{-1} A_{i+1} for the ladder with identity rails."""
    return [rungs[i].inv() * rungs[i + 1] for i in range(len(rungs) - 1)]


def stokes_product(faces):
    out = QONE
    for h in faces:
        out = out * h
    return out


def run():
    # ---- T1 / C1: R, iota K, iota RK anticommute, square to -I
    iK = mscale(K_, IOTA)
    iRK = mscale(RK, IOTA)
    gens = [R_, iK, iRK]
    c1 = all(meq(mmul(g, g), NEG_I2) for g in gens) and all(
        meq(madd(mmul(gens[i], gens[j]), mmul(gens[j], gens[i])),
            mscale(I2, ZERO)) for i in range(3) for j in range(3) if i != j)
    # quaternion relations i j = k etc. under the identification
    c1 = c1 and meq(mmul(R_, iK), iRK) and meq(mmul(iK, iRK), R_) \
        and meq(mmul(iRK, R_), iK)

    # ---- C2: homomorphism on a rational grid
    c2 = True
    for _ in range(40):
        q1, q2 = random_unit(), random_unit()
        if not meq(to_matrix(q1 * q2), mmul(to_matrix(q1), to_matrix(q2))):
            c2 = False
    # ---- C3: det = Delta_par + Delta_perp = N ; unit sphere is SU(2)
    c3 = True
    for _ in range(40):
        a, b, c, d = (F(random.randint(-5, 5), random.randint(1, 4))
                      for _ in range(4))
        M = emk_block(a, b, c, d)
        det = mdet(M)
        if not (det == GQ(a * a + b * b + c * c + d * d)):
            c3 = False
        # untwisted block has the EMK-1 split-signature determinant
        Mu = emk_block(a, b, c, d, twist=False)
        if not (mdet(Mu) == GQ((a * a - b * b) + (c * c - d * d))):
            c3 = False
    for _ in range(20):
        q = random_unit()
        M = to_matrix(q)
        if not (meq(mmul(mdagger(M), M), I2) and mdet(M) == ONE):
            c3 = False
        # character = twice the identity-sector coefficient
        if not (M[0][0] + M[1][1] == GQ(q.trace())):
            c3 = False

    # ---- T2 / T3 / C4 fabric
    c4 = True
    t2 = True
    for m in M_RANGE:
        rungs = [random_unit() for _ in range(m)]
        faces = ladder_faces(rungs)
        # T2: face trace equals bridge trace Tr(A_i A_{i+1}^{-1})
        for i, h in enumerate(faces):
            if h.trace() != (rungs[i] * rungs[i + 1].inv()).trace():
                t2 = False
        # T3: telescoping to boundary holonomy
        if stokes_product(faces) != rungs[0].inv() * rungs[-1]:
            c4 = False
        # tamper: flip orientation of one interior rung
        if m >= 3:
            bad = list(rungs)
            bad_faces = list(faces)
            j = m // 2
            # replace face j-1 by its version with rung j reversed
            bad_faces[j - 1] = rungs[j - 1].inv() * rungs[j].inv()
            if stokes_product(bad_faces) == rungs[0].inv() * rungs[-1]:
                # generic data: equality would mean the tamper did not bite
                c4 = False
    # ---- T4 / C5 / C6 residue form of the weight
    c5 = True
    c6 = True
    residues_demo = None
    for m in (3, 5):
        rungs = [random_unit() for _ in range(m)]
        faces = ladder_faces(rungs)
        rhos = [h.residue() for h in faces]
        if not all(0 <= r <= 2 for r in rhos):
            c5 = False
        if QONE.residue() != 0 or Quat(-1, 0, 0, 0).residue() != 2:
            c5 = False
        # direct: prod exp[(kappa/2) Tr H] ; residue: e^{kappa(m-1)} e^{-kappa sum rho}
        direct = Iv(F(1))
        for h in faces:
            e = exp_point(KAPPA / 2 * h.trace())
            direct = direct * Iv(e.lo, e.hi)
        e1 = exp_point(KAPPA * (m - 1))
        e2 = exp_point(-KAPPA * sum(rhos))
        res_form = Iv(e1.lo, e1.hi) * Iv(e2.lo, e2.hi)
        if direct.hi < res_form.lo or res_form.hi < direct.lo:
            c6 = False
        if residues_demo is None:
            residues_demo = {"m": m, "face_residues": [str(r) for r in rhos],
                             "total_residue": str(sum(rhos)),
                             "weight_direct": [_dec(direct.lo, 15),
                                               _dec(direct.hi, 15)],
                             "weight_residue_form": [_dec(res_form.lo, 15),
                                                     _dec(res_form.hi, 15)]}

    ok = c1 and c2 and c3 and c4 and t2 and c5 and c6
    cert = {
        "certificate_type": "YMF1_CHAIN_AS_RECOGNITION_FABRIC_NATIVE_CARRIER",
        "claim_status": "native_carrier_certified__no_gap_claim",
        "source": ["vault APPENDIX_DISCRETE_HOLONOMY_FABRICS_AND_RECOGNITION_"
                   "STOKES (fabric, face holonomy, Stokes closure)",
                   "EMK-1 determinant seam ladder (pinned, Publications)",
                   "F00-E native exponential / iota_Sigma (pinned, RKF)",
                   "YM-11 gate 1 (tree gauge fixing) — consumed"],
        "theorems": {
            "T1_gauge_group_is_EMK_block_over_C_Sigma":
                "span{I, R, iota K, iota RK} = Hamilton quaternions; "
                "det = Delta_par + Delta_perp = N_Sigma; SU(2) = unit sphere; "
                "chi_1/2 = 2 x identity-sector coefficient",
            "T2_chain_is_ladder_fabric":
                "bridge face weight = ladder plaquette weight, H_l = A_i^{-1} A_{i+1}",
            "T3_recognition_stokes_exact_on_chain":
                "prod H_l = A_1^{-1} A_m; interior rungs cancel exactly "
                "(E_int = 0); orientation tamper breaks it",
            "T4_action_is_face_residue_sum":
                "weight = Exp(kappa(m-1)) Exp(-kappa sum rho_l), "
                "rho_l = 1 - Tr H_l / 2 in [0,2]",
        },
        "remaining_shadow": "Haar integration <-> Phi_Sigma stationary-path "
                            "extraction (pairing delta_{a=b=c}/d_a): DECLARED",
        "native_form_of_YM21": "does the boundary-holonomy residue of a fabric "
                               "with exact interior cancellation control the "
                               "face-residue sum uniformly in the number of "
                               "faces?",
        "demo": residues_demo,
        "controls": {
            "C1_quaternion_relations_over_Q_iota": bool(c1),
            "C2_homomorphism_quat_vs_matrix": bool(c2),
            "C3_det_is_norm_unit_sphere_is_SU2": bool(c3),
            "C4_stokes_exact_and_tamper_bites": bool(c4),
            "T2_face_trace_is_bridge_trace": bool(t2),
            "C5_residue_range_and_endpoints": bool(c5),
            "C6_two_route_weight_overlap": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YMF1_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YMF1.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    print("demo:", cert["demo"])
    print("sha256:", sha)
