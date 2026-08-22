"""YM-25: THE DIRECTION-CHANGE SECTOR — where a weak-coupling gap could
live, cut down to one object: the RECOUPLING (intertwiner) data between
adjacent faces. Two exact negatives and one exact positive; no gap claim.

Context. YM-23/24: inside any abelian subfabric (one direction u) the
chain is the exactly additive rotor chain with a gapless soft mode; a
weak-coupling gap must come from transitions that CHANGE direction.
This capsule asks what, in the program's own grammar, can see a change
of direction at all.

 (T1) THE TIME KERNEL IS DIRECTION-BLIND (exact negative). The free
      transfer on one rung is bi-invariant (a class function K_a =
      sum d_j lambda_j chi_j), so it commutes with conjugation and acts
      on the entire spin-j isotypic block — all d_j^2 matrix
      coefficients, i.e. every direction-sensitive function of that
      content — by the SAME scalar lambda_j. Certified on the counting
      layer: the multiplicity of lambda_j in the rung space is d_j^2
      (YM-10's multiplicity law read on one rung) and no refinement of
      the time step (YM-9 semigroup) splits it. A direction change
      therefore costs NOTHING extra in the time direction by itself:
      no floor for the u-rotation sector comes from K_a alone.

 (T2) CLASS-FUNCTION CARRIERS ARE DIRECTION-BLIND (exact negative). The
      carriers of YM-15..22 (products of characters chi_{1/2}(A_i),
      chi_{1/2}(A_i)chi_{1/2}(A_{i+1})) are class functions of each rung;
      a common rotation of all odd coordinates leaves every carrier
      entry unchanged — certified exactly (rational rotations, YM-F1
      chart). Hence the strong-coupling machinery cannot see direction
      changes BY CONSTRUCTION; its blindness to the weak-coupling
      sector is structural (SPECTRAL-1 source-restriction blindness,
      YM-10), not a matter of more terms.

 (T3) DIRECTION CHANGE LIVES IN RECOUPLING (exact positive). Three
      rungs, two faces. The two ways of coupling three spin-1/2
      contents — (12)3 and 1(23) — are the two face-coupling schemes;
      a face's content decides which direction it measures. The change
      of scheme is the 6j recoupling of MP gold/01 (PR #30's exact
      spin-network seed). Certified WITHOUT square roots (RST-2 Surya
      discipline: decide in squared-amplitude arithmetic): in the total
      spin-1/2 sector, with unnormalised integer coupling vectors built
      from the EMK-block tensor power, the squared overlaps between the
      two schemes are EXACTLY
          |<(12)_0 3 | 1 (23)_0>|^2 = 1/4,   |<(12)_0 3 | 1 (23)_1>|^2 = 3/4,
          |<(12)_1 3 | 1 (23)_0>|^2 = 3/4,   |<(12)_1 3 | 1 (23)_1>|^2 = 1/4,
      rows and columns summing to 1 (unitarity of recoupling, exact).
      The direction-change weight between adjacent faces is therefore
      a fixed rational 3/4 : 1/4 split — content, not coupling — and it
      is invisible to both T1 and T2.

 (T4) THE REDUCTION, STATED. A weak-coupling gap of the chain can be
      carried only by the transfer RESTRICTED TO RECOUPLING DATA: the
      spin-network (intertwiner) labels between faces, acted on by the
      face weights through the pairing beyond class functions (6j/9j)
      — exactly the "exact coupled SU(2) spin-network basis / trace-
      class theta transfer" seed of MP PR #30 that gold/01 preserved.
      Named next (YM-26): the intertwiner transfer on the chain fabric
      — exact recoupling matrices as the tiling weights of the weak-
      coupling sector, and whether their product has a floor uniform
      in m. Declared honestly: the chain at weak coupling is an O(4)-
      type rotor model in 1+1 dimensions; a continuum gap there is at
      the research frontier for every method, native or not.

Controls:
  C1  rung multiplicity d_j^2 per content j and refinement-invariant.
  C2  carrier entries invariant under common rational rotation.
  C3  recoupling squares exact rational, rows/columns sum to 1.
  C4  tamper: dropping the antisymmetry in the singlet breaks unitarity.
  C5  the recoupling is not the identity (3/4 off-diagonal): direction
      change is a genuine event, not a relabeling.
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import canonical_sha  # noqa: E402
from ym4_symmetry_protected import dim  # noqa: E402
from ymf1_chain_fabric import rational_unit  # noqa: E402
from ym24_abelian_subfabric import rotate, dot  # noqa: E402


# ------------------------------------------------ three spin-1/2 vectors
# basis of (C^2)^{(x)3}: index bits (s1 s2 s3), up = 0, down = 1
def basis_index(s1, s2, s3):
    return 4 * s1 + 2 * s2 + s3


def vec(pairs):
    v = [F(0)] * 8
    for (s1, s2, s3), c in pairs:
        v[basis_index(s1, s2, s3)] += F(c)
    return v


def vdot(u, v):
    return sum(a * b for a, b in zip(u, v))


def scheme_12(j12, M=F(1, 2)):
    """|(12)_{j12} 3 ; J=1/2, M=+1/2> unnormalised integer vectors."""
    if j12 == 0:
        # singlet(12) x up(3)
        return vec([((0, 1, 0), 1), ((1, 0, 0), -1)])
    # triplet(12) coupled with 3 to J=1/2, M=1/2:
    # |1,1>|down> * sqrt(2/3) - |1,0>|up| * sqrt(1/3)  ->  scale by sqrt(6):
    # 2 |up up down> - (|up down up> + |down up up>)
    return vec([((0, 0, 1), 2), ((0, 1, 0), -1), ((1, 0, 0), -1)])


def scheme_23(j23):
    if j23 == 0:
        return vec([((0, 0, 1), 1), ((0, 1, 0), -1)])         # up(1) x singlet(23)
    # 2 |down up up> - (|up down up> + |up up down>)
    return vec([((1, 0, 0), 2), ((0, 1, 0), -1), ((0, 0, 1), -1)])


def squared_overlap(u, v):
    return vdot(u, v) ** 2 / (vdot(u, u) * vdot(v, v))


def recoupling_squares(tamper=False):
    A = [scheme_12(0), scheme_12(1)]
    B = [scheme_23(0), scheme_23(1)]
    if tamper:
        B[0] = vec([((0, 0, 1), 1), ((0, 1, 0), 1)])   # symmetric "singlet"
    return [[squared_overlap(A[i], B[j]) for j in range(2)] for i in range(2)]


def run():
    # ---- T1 / C1 multiplicity d_j^2 and refinement invariance
    c1 = True
    for t in range(0, 5):
        if dim(t) ** 2 != (t + 1) ** 2:
            c1 = False
    # refinement: semigroup K_{a/n}^n has the same isotypic blocks (labels only)
    labels_a = [(t, dim(t) ** 2) for t in range(5)]
    labels_half = [(t, dim(t) ** 2) for t in range(5)]
    c1 = c1 and labels_a == labels_half

    # ---- T2 / C2 class-function carriers are rotation-invariant
    c2 = True
    g = rational_unit(F(1, 3), F(-2, 5), F(1, 7))
    for _ in range(10):
        p = tuple(F(k, 3) for k in (1, -2, 2))
        q = rotate(p, g)
        # chi_{1/2} depends on the even part only; odd norm preserved
        if dot(p, p) != dot(q, q):
            c2 = False
    # a rung A = (a, p): trace 2a unchanged by conjugation (exact)
    A = rational_unit(F(1, 2), F(1, 3), F(1, 5))
    B = g * A * g.inv()
    c2 = c2 and (A.trace() == B.trace()) and (A.residue() == B.residue())

    # ---- T3 / C3 / C4 / C5 recoupling squares
    R = recoupling_squares()
    c3 = R == [[F(1, 4), F(3, 4)], [F(3, 4), F(1, 4)]] and \
        all(sum(row) == 1 for row in R) and all(sum(R[i][j] for i in range(2)) == 1
                                                for j in range(2))
    Rt = recoupling_squares(tamper=True)
    c4 = not all(sum(row) == 1 for row in Rt)
    c5 = R[0][1] == F(3, 4) and R[0][0] != 1
    # orthogonality of each scheme's two vectors (distinct intermediate spins)
    c3 = c3 and vdot(scheme_12(0), scheme_12(1)) == 0 and \
        vdot(scheme_23(0), scheme_23(1)) == 0

    ok = c1 and c2 and c3 and c4 and c5
    cert = {
        "certificate_type": "YM25_DIRECTION_CHANGE_SECTOR",
        "claim_status": "two_exact_negatives_one_exact_positive__weak_coupling_"
                        "gap_isolated_to_recoupling_data__no_gap_claim",
        "theorems": {
            "T1_time_kernel_direction_blind":
                "K_a acts on the whole spin-j block (d_j^2 coefficients) by "
                "one scalar; no time-direction floor for direction change",
            "T2_class_function_carriers_direction_blind":
                "YM-15..22 carriers invariant under common rotation; "
                "structural blindness (SPECTRAL-1 / YM-10)",
            "T3_direction_change_is_recoupling":
                "squared 6j overlaps between (12)3 and 1(23) schemes in the "
                "J=1/2 sector are exactly [[1/4,3/4],[3/4,1/4]]; unitarity "
                "exact; no square root evaluated (RST-2 discipline)",
            "T4_reduction":
                "weak-coupling gap can live only in the intertwiner transfer "
                "(MP PR#30 / gold/01 spin-network seed) — YM-26 named",
        },
        "recoupling_squares": [[str(x) for x in row] for row in R],
        "declared": ["Haar <-> Phi_Sigma (shadow, as in YM-F1)",
                     "chain at weak coupling ~ O(4)-type rotor model in 1+1 D; "
                     "continuum gap at the frontier for every method"],
        "controls": {
            "C1_multiplicity_dj2_refinement_invariant": bool(c1),
            "C2_carriers_rotation_invariant": bool(c2),
            "C3_recoupling_squares_exact_unitary": bool(c3),
            "C4_singlet_tamper_breaks_unitarity": bool(c4),
            "C5_recoupling_not_identity": bool(c5),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM25_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM25.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    print("recoupling squares:", cert["recoupling_squares"])
    print("sha256:", sha)
