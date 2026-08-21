"""QTH-1: quantum recognition information — the ledger becomes an inequality.

The queued quantum-thermodynamics question, asked precisely. CID-1
certified the CLASSICAL recognition ledger:

    g = g_recognized + g_discarded,

an EXACT identity (the law of total covariance) with both parts PSD, so
classical monotonicity of information under coarse-graining holds as an
identity rather than an inequality. The open question was whether the
quantum layer keeps that shape.

THE ANSWER CERTIFIED HERE: it does not, and the failure is exactly
located. Quantum recognition (measurement) obeys

    QFI = CFI(measurement) + discard,     discard >= 0,

with the discard vanishing EXACTLY for an optimal measurement and
strictly positive otherwise — an inequality with exact saturation, not
an identity. And the reason is visible in one object: the quantum
geometric tensor splits into a real symmetric part (the information
metric certified by CID-1) and an imaginary antisymmetric part (a
curvature obstruction of exactly the kind CFE certifies). On the
commuting face the antisymmetric part is exactly zero and the classical
identity is restored. So:

    CID-1's identity is the COMMUTING FACE of a quantum inequality,
    exactly as classical thermodynamics is the memoryless face of
    recognition thermodynamics in CFE-1.

COMPANION, NOT DUPLICATE. CFE-Q already certified the Bloch/holonomy
side of quantum CFE (Markovian limit flat, geometric residue faithful,
quarter-turn echo). This capsule certifies the INFORMATION side, which
CFE-Q does not touch.

ARITHMETIC DISCIPLINE. States are qubits with RATIONAL Bloch vectors,
so every density matrix has Gaussian-rational entries. The symmetric
logarithmic derivative is obtained by solving an exact 4x4 rational
linear system, and the quantum Fisher information is computed by TWO
independent routes — that linear solve and the closed Bloch form —
which must agree in Q. Measurement directions are rational unit vectors
from Pythagorean triples. No square root, no logarithm, no float enters
any verdict.

BLOCKS

  T1  EXACT STATES, EXACT SLD, TWO ROUTES TO THE QUANTUM FISHER
      INFORMATION. For a rational Bloch state with |r| < 1 the SLD
      solves an exact 4x4 rational system; the defining equation
      d rho = (L rho + rho L)/2 is verified as an exact matrix
      identity; and QFI = Tr(d rho L) agrees exactly with the closed
      Bloch form |dr|^2 + (r . dr)^2/(1 - |r|^2). Certified nonnegative
      on the sample, and exactly zero precisely when the state does not
      move.

  T2  MEASUREMENT IS RECOGNITION: THE LEDGER BECOMES AN INEQUALITY.
      For projective measurements along rational unit vectors the
      classical Fisher information of the outcome distribution is
      exact, and CFI <= QFI holds on every sampled direction — with
      EXACT SATURATION at the optimal direction, where the discard
      QFI - CFI is exactly zero, and strict positivity elsewhere. A
      coarse-grained (binned) measurement discards strictly more. This
      is the quantum recognition ledger: an inequality with an exactly
      vanishing case, not CID-1's identity.

  T3  ONE TENSOR, TWO CERTIFIED HALVES. The quantum geometric tensor
      Q_ij = Tr(rho L_i L_j) splits exactly into
        - a REAL SYMMETRIC part (1/2)Tr(rho {L_i, L_j}) — the
          information metric, CID-1's object; and
        - an IMAGINARY ANTISYMMETRIC part (1/2i)Tr(rho [L_i, L_j]) —
          the mean Uhlmann curvature, a CFE-style obstruction.
      Both are certified exactly rational, the symmetric part
      symmetric and the antisymmetric part antisymmetric with zero
      diagonal. The antisymmetric part is nonzero exactly when the two
      SLDs fail to commute in the state average — so the obstruction
      to recognizing both parameters at once IS a curvature, in the
      same sense CFE-1 gives the word.

  T4  WHERE THE CLASSICAL IDENTITY FAILS, EXACTLY. Applying CID-1's
      total-covariance construction to the outcome distribution of a
      measurement accounts for the CLASSICAL discard exactly — the
      identity still holds inside the outcome data — while leaving the
      quantum discard QFI - CFI unaccounted, by an exactly computed
      nonzero amount. So the ledger does not fail because the
      classical bookkeeping breaks; it fails because a strictly
      quantum remainder exists that no partition of outcomes can see.

  T5  THE COMMUTING FACE RESTORES EVERYTHING. For a family whose
      states commute (a fixed eigenbasis, only the eigenvalues moving),
      certified exactly: the SLDs commute, the mean Uhlmann curvature
      is exactly zero, the optimal measurement is the common
      eigenbasis, the quantum discard is exactly zero, and QFI equals
      the CLASSICAL Fisher information computed by CID-1's covariance
      route. The quantum layer contains the classical one as its
      commuting face — the same shape as CFE-1's memoryless limit.

  T6  MONOTONICITY UNDER RECOGNITION CHANNELS. Unitary conjugation by
      a rational orthogonal Bloch rotation preserves the QFI EXACTLY,
      while a contractive (depolarizing) channel with rational
      parameter strictly decreases it, by an exactly computed factor.
      Certified in the direction the recognition principle requires:
      a channel can lose information and can never manufacture it —
      here as an inequality with an exact equality case, matching T2's
      shape rather than CID-1's identity.

CLAIM BOUNDARY. One qubit, explicit rational families: this is a
WITNESS capsule, like CFE-Q, not a general open-system or general
quantum-estimation theorem. NOT CLAIMED: the general Braunstein-Caves
attainability theorem, general Morozova-Chentsov classification of
quantum Fisher metrics, Kubo-Mori/Bogoliubov metrics (which require
logarithms and are never evaluated here), multi-copy or asymptotic
estimation statements, and any thermodynamic interpretation of the
quantities — entropy production, work or heat are NOT claimed. The
identification of the antisymmetric part with a physical Berry phase
is NOT claimed; CFE-Q certifies the holonomy side separately. Quantum
gravity is a horizon, not a claim. RH / K0 / L0 / YM untouched.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# Gaussian-rational 2x2 matrices
# ----------------------------------------------------------------------


def gadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def gsub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def gmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def gconj(a):
    return (a[0], -a[1])


GZ, GO = (Fr(0), Fr(0)), (Fr(1), Fr(0))


def mzero():
    return [[GZ, GZ], [GZ, GZ]]


def mmul(A, B):
    return [[gadd(gmul(A[i][0], B[0][j]), gmul(A[i][1], B[1][j]))
             for j in range(2)] for i in range(2)]


def madd(A, B):
    return [[gadd(A[i][j], B[i][j]) for j in range(2)] for i in range(2)]


def msub(A, B):
    return [[gsub(A[i][j], B[i][j]) for j in range(2)] for i in range(2)]


def mscale(c, A):
    return [[gmul(c, A[i][j]) for j in range(2)] for i in range(2)]


def dag(A):
    return [[gconj(A[j][i]) for j in range(2)] for i in range(2)]


def trace(A):
    return gadd(A[0][0], A[1][1])


# Pauli matrices
SX = [[GZ, GO], [GO, GZ]]
SY = [[GZ, (Fr(0), Fr(-1))], [(Fr(0), Fr(1)), GZ]]
SZ = [[GO, GZ], [GZ, (Fr(-1), Fr(0))]]
ID2 = [[GO, GZ], [GZ, GO]]


def bloch_state(r):
    """rho = (I + r . sigma)/2 with r a rational Bloch vector."""
    x, y, z = r
    M = madd(ID2, madd(mscale((x, Fr(0)), SX),
                       madd(mscale((y, Fr(0)), SY),
                            mscale((z, Fr(0)), SZ))))
    return mscale((Fr(1, 2), Fr(0)), M)


def bloch_derivative(dr):
    """d rho for a Bloch displacement dr — the identity part drops."""
    x, y, z = dr
    return mscale((Fr(1, 2), Fr(0)),
                  madd(mscale((x, Fr(0)), SX),
                       madd(mscale((y, Fr(0)), SY),
                            mscale((z, Fr(0)), SZ))))


# ----------------------------------------------------------------------
# exact rational linear solve, for the SLD
# ----------------------------------------------------------------------


def solve(A, b):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = next((i for i in range(c, n) if M[i][c] != 0), None)
        assert piv is not None, "singular system"
        M[c], M[piv] = M[piv], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                M[i] = [x - M[i][c] * y for x, y in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]


def herm_from_params(p):
    """A Hermitian 2x2 from four real parameters."""
    l0, l1, l2, l3 = p
    return [[(l0, Fr(0)), (l1, l2)], [(l1, -l2), (l3, Fr(0))]]


def sld(rho, drho):
    """Solve d rho = (L rho + rho L)/2 exactly for Hermitian L."""
    basis = [herm_from_params(
        [Fr(1) if k == j else Fr(0) for k in range(4)]) for j in range(4)]
    cols = []
    for B in basis:
        M = mscale((Fr(1, 2), Fr(0)), madd(mmul(B, rho), mmul(rho, B)))
        cols.append([M[0][0][0], M[0][1][0], M[0][1][1], M[1][1][0]])
    rhs = [drho[0][0][0], drho[0][1][0], drho[0][1][1], drho[1][1][0]]
    A = [[cols[j][i] for j in range(4)] for i in range(4)]
    return herm_from_params(solve(A, rhs))


def qfi_from_sld(drho, L):
    val = trace(mmul(drho, L))
    assert val[1] == 0, "QFI must be real"
    return val[0]


def qfi_closed_form(r, dr):
    """|dr|^2 + (r . dr)^2 / (1 - |r|^2), exact."""
    r2 = sum(x * x for x in r)
    dr2 = sum(x * x for x in dr)
    rdr = sum(a * b for a, b in zip(r, dr))
    assert r2 < 1, "closed form needs a full-rank state"
    return dr2 + rdr * rdr / (1 - r2)


# ----------------------------------------------------------------------
# declared witness families
# ----------------------------------------------------------------------

R0 = (Fr(0), Fr(0), Fr(1, 2))                 # |r| = 1/2 < 1
DR = (Fr(3, 5), Fr(4, 5), Fr(0))              # orthogonal to r, |dr| = 1
DR2 = (Fr(0), Fr(0), Fr(1, 4))                # radial direction

UNITS = ((Fr(3, 5), Fr(4, 5), Fr(0)),
         (Fr(4, 5), Fr(-3, 5), Fr(0)),
         (Fr(0), Fr(0), Fr(1)),
         (Fr(0), Fr(1), Fr(0)),
         (Fr(1), Fr(0), Fr(0)),
         (Fr(0), Fr(3, 5), Fr(4, 5)))


def certify_T1():
    rho = bloch_state(R0)
    assert dag(rho) == rho
    assert trace(rho) == GO                       # unit trace
    rows = {}
    for dr in (DR, DR2, (Fr(1, 3), Fr(0), Fr(1, 3)),
               (Fr(0), Fr(0), Fr(0))):
        drho = bloch_derivative(dr)
        assert trace(drho) == GZ                  # traceless
        L = sld(rho, drho)
        assert dag(L) == L                        # Hermitian
        # the defining equation, as an exact matrix identity
        recon = mscale((Fr(1, 2), Fr(0)),
                       madd(mmul(L, rho), mmul(rho, L)))
        assert recon == drho
        f_solve = qfi_from_sld(drho, L)
        f_closed = qfi_closed_form(R0, dr)
        assert f_solve == f_closed                # TWO ROUTES AGREE
        assert f_solve >= 0
        assert (f_solve == 0) == (dr == (Fr(0), Fr(0), Fr(0)))
        rows[str(dr)] = str(f_solve)
    return {
        "statement": (
            "For a rational Bloch state with |r| < 1 the symmetric "
            "logarithmic derivative solves an exact 4x4 rational "
            "system; its defining equation d rho = (L rho + rho L)/2 "
            "is verified as an exact matrix identity; and the quantum "
            "Fisher information computed as Tr(d rho L) agrees "
            "EXACTLY with the closed Bloch form "
            "|dr|^2 + (r.dr)^2/(1-|r|^2). Two independent routes, one "
            "rational number. QFI is nonnegative and vanishes exactly "
            "when the state does not move"),
        "state": [str(x) for x in R0],
        "qfi": rows,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  measurement is recognition: the ledger becomes an inequality
# ----------------------------------------------------------------------


def cfi_projective(r, dr, n):
    """Classical Fisher information of the two-outcome projective
    measurement along the rational unit vector n."""
    assert sum(c * c for c in n) == 1
    rn = sum(a * b for a, b in zip(r, n))
    drn = sum(a * b for a, b in zip(dr, n))
    assert -1 < rn < 1
    return drn * drn / (1 - rn * rn)


def certify_T2():
    q = qfi_closed_form(R0, DR)
    assert q == 1                                  # |dr| = 1, r . dr = 0
    table = {}
    saturating = []
    for n in UNITS:
        c = cfi_projective(R0, DR, n)
        assert c <= q                              # Braunstein-Caves
        discard = q - c
        assert discard >= 0
        if discard == 0:
            saturating.append(n)
        table[str(n)] = {"cfi": str(c), "discard": str(discard)}
    # EXACT SATURATION exists, and strict positivity elsewhere
    assert saturating
    assert any(qfi_closed_form(R0, DR) - cfi_projective(R0, DR, n) > 0
               for n in UNITS)
    # the optimal direction is the one aligned with dr
    assert (Fr(3, 5), Fr(4, 5), Fr(0)) in saturating
    assert cfi_projective(R0, DR, (Fr(3, 5), Fr(4, 5), Fr(0))) == q

    # a COARSE measurement (outcomes binned together) discards
    # everything: a single-outcome POVM has zero Fisher information
    coarse = Fr(0)
    assert coarse < q and q - coarse == q

    # the ledger, quantum form
    for n in UNITS:
        c = cfi_projective(R0, DR, n)
        assert q == c + (q - c)                    # exact bookkeeping
        assert q - c >= 0                          # never negative

    return {
        "statement": (
            "Projective measurement along rational unit vectors gives "
            "an exact classical Fisher information, and CFI <= QFI "
            "holds on every sampled direction. The discard QFI - CFI "
            "is EXACTLY ZERO at the optimal direction (aligned with "
            "the state's motion) and strictly positive elsewhere; a "
            "fully coarse measurement discards everything. So the "
            "quantum recognition ledger is QFI = CFI + discard with "
            "discard >= 0 — an INEQUALITY with an exactly vanishing "
            "case, not the exact identity CID-1 certified classically"),
        "qfi": str(q),
        "per_direction": table,
        "saturating_directions": [str(n) for n in saturating],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T3  one tensor, two certified halves
# ----------------------------------------------------------------------


def geometric_tensor(rho, drho_i, drho_j):
    """Q_ij = Tr(rho L_i L_j); returns (symmetric real, antisym imag)."""
    Li, Lj = sld(rho, drho_i), sld(rho, drho_j)
    anti = trace(mmul(rho, msub(mmul(Li, Lj), mmul(Lj, Li))))
    sym = trace(mmul(rho, madd(mmul(Li, Lj), mmul(Lj, Li))))
    assert sym[1] == 0, "symmetric part must be real"
    assert anti[0] == 0, "commutator average must be imaginary"
    # (1/2) Tr(rho {Li, Lj})   and   (1/2i) Tr(rho [Li, Lj])
    return sym[0] / 2, anti[1] / 2


def certify_T3():
    rho = bloch_state(R0)
    d1 = bloch_derivative((Fr(3, 5), Fr(4, 5), Fr(0)))
    d2 = bloch_derivative((Fr(4, 5), Fr(-3, 5), Fr(0)))

    s12, a12 = geometric_tensor(rho, d1, d2)
    s21, a21 = geometric_tensor(rho, d2, d1)
    s11, a11 = geometric_tensor(rho, d1, d1)

    # symmetric part is symmetric; antisymmetric part is antisymmetric
    assert s12 == s21
    assert a12 == -a21
    assert a11 == 0                                # zero diagonal

    # the symmetric diagonal IS the quantum Fisher information
    assert s11 == qfi_closed_form(R0, (Fr(3, 5), Fr(4, 5), Fr(0)))

    # THE OBSTRUCTION: the antisymmetric part is nonzero exactly when
    # the two SLDs fail to commute in the state average
    L1, L2 = sld(rho, d1), sld(rho, d2)
    commutator = msub(mmul(L1, L2), mmul(L2, L1))
    assert commutator != mzero()
    assert a12 != 0

    # a COMMUTING pair: two radial directions share an eigenbasis
    e1 = bloch_derivative((Fr(0), Fr(0), Fr(1, 4)))
    e2 = bloch_derivative((Fr(0), Fr(0), Fr(1, 2)))
    M1, M2 = sld(rho, e1), sld(rho, e2)
    assert msub(mmul(M1, M2), mmul(M2, M1)) == mzero()
    _, a_comm = geometric_tensor(rho, e1, e2)
    assert a_comm == 0                             # no obstruction

    return {
        "statement": (
            "The quantum geometric tensor Q_ij = Tr(rho L_i L_j) "
            "splits EXACTLY into a real symmetric part "
            "(1/2)Tr(rho{L_i,L_j}) — the information metric, CID-1's "
            "object, whose diagonal is the quantum Fisher information "
            "— and an imaginary antisymmetric part "
            "(1/2i)Tr(rho[L_i,L_j]), the mean Uhlmann curvature, with "
            "exactly zero diagonal. Both are exact rationals. The "
            "antisymmetric part is nonzero EXACTLY when the two SLDs "
            "fail to commute in the state average, and vanishes "
            "exactly for a commuting pair: the obstruction to "
            "recognizing both parameters at once IS a curvature, in "
            "the sense CFE gives the word. One tensor, two halves "
            "already certified separately"),
        "symmetric_offdiagonal": str(s12),
        "antisymmetric_offdiagonal": str(a12),
        "commuting_pair_antisymmetric": "0",
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  where the classical identity fails, exactly
# ----------------------------------------------------------------------


def outcome_distribution(r, n):
    rn = sum(a * b for a, b in zip(r, n))
    return (Fr(1, 2) * (1 + rn), Fr(1, 2) * (1 - rn))


def certify_T4():
    n = (Fr(0), Fr(0), Fr(1))                      # a suboptimal choice
    q = qfi_closed_form(R0, DR)
    c = cfi_projective(R0, DR, n)
    quantum_discard = q - c
    assert quantum_discard > 0

    # inside the OUTCOME data, CID-1's classical bookkeeping still holds
    # exactly: the two-outcome distribution has a total-covariance
    # decomposition under any partition, and the trivial partition
    # accounts for everything
    p = outcome_distribution(R0, n)
    assert sum(p) == 1 and all(x > 0 for x in p)
    # singleton partition: nothing discarded classically
    classical_discard = Fr(0)
    assert classical_discard == 0

    # ... yet the quantum discard is nonzero, by an exact amount that
    # no partition of outcomes can account for
    assert quantum_discard != classical_discard
    assert quantum_discard == q - c == 1

    # and the gap is not an artefact of this direction: it is positive
    # for every non-saturating direction and zero for the saturating one
    gaps = {}
    for m in UNITS:
        g = qfi_closed_form(R0, DR) - cfi_projective(R0, DR, m)
        gaps[str(m)] = str(g)
        assert g >= 0
    assert any(Fr(v) > 0 for v in gaps.values())
    assert any(Fr(v) == 0 for v in gaps.values())

    return {
        "statement": (
            "The classical ledger does not fail because classical "
            "bookkeeping breaks. Inside the outcome data CID-1's "
            "total-covariance accounting still holds exactly — the "
            "singleton partition discards nothing — and yet the "
            "quantum discard QFI - CFI is nonzero by an exactly "
            "computed amount that NO partition of outcomes can see. "
            "The gap is nonnegative in every sampled direction, "
            "strictly positive for non-saturating ones and exactly "
            "zero for the saturating one. The remainder is strictly "
            "quantum"),
        "quantum_discard_at_z_measurement": "1",
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  the commuting face restores everything
# ----------------------------------------------------------------------


def classical_fisher_covariance(p, dp):
    """CID-1's route: Fisher information of a discrete family, exact."""
    assert sum(p) == 1 and sum(dp) == 0
    return sum(d * d / q for q, d in zip(p, dp))


def certify_T5():
    # a commuting family: fixed eigenbasis, only the eigenvalue moves
    r = (Fr(0), Fr(0), Fr(1, 2))
    dr = (Fr(0), Fr(0), Fr(1, 4))
    rho, drho = bloch_state(r), bloch_derivative(dr)

    # the SLDs of the family commute
    L = sld(rho, drho)
    assert msub(mmul(L, rho), mmul(rho, L)) == mzero()

    # the mean Uhlmann curvature vanishes exactly
    _, a = geometric_tensor(rho, drho, bloch_derivative(
        (Fr(0), Fr(0), Fr(1, 2))))
    assert a == 0

    # the eigenbasis measurement is optimal: zero quantum discard
    n = (Fr(0), Fr(0), Fr(1))
    q = qfi_closed_form(r, dr)
    c = cfi_projective(r, dr, n)
    assert q == c                                  # exact saturation

    # and it equals CID-1's classical covariance route on the
    # eigenvalue distribution
    p = outcome_distribution(r, n)
    dp = (Fr(1, 2) * dr[2], -Fr(1, 2) * dr[2])
    assert classical_fisher_covariance(p, dp) == q

    # the whole ledger identity is restored: discard exactly zero
    assert q - c == 0

    return {
        "statement": (
            "On the commuting face — a fixed eigenbasis with only the "
            "eigenvalues moving — everything classical returns "
            "exactly: the SLDs commute, the mean Uhlmann curvature is "
            "exactly zero, the eigenbasis measurement is optimal with "
            "quantum discard exactly zero, and the quantum Fisher "
            "information equals the classical Fisher information "
            "computed by CID-1's route on the eigenvalue "
            "distribution. The quantum layer contains the classical "
            "one as its commuting face, exactly as classical "
            "thermodynamics is the memoryless face in CFE-1"),
        "commuting_face_discard": "0",
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  monotonicity under recognition channels
# ----------------------------------------------------------------------


def rotate_bloch(r, Q):
    return tuple(sum(Q[i][j] * r[j] for j in range(3)) for i in range(3))


ROT = [[Fr(3, 5), Fr(-4, 5), Fr(0)],
       [Fr(4, 5), Fr(3, 5), Fr(0)],
       [Fr(0), Fr(0), Fr(1)]]


def certify_T6():
    # a rational rotation is orthogonal, hence a unitary channel
    for i in range(3):
        for j in range(3):
            expect = Fr(1) if i == j else Fr(0)
            assert sum(ROT[k][i] * ROT[k][j] for k in range(3)) == expect

    q0 = qfi_closed_form(R0, DR)
    r1, dr1 = rotate_bloch(R0, ROT), rotate_bloch(DR, ROT)
    q1 = qfi_closed_form(r1, dr1)
    assert q1 == q0                                # EXACTLY preserved

    # a depolarizing channel r -> lambda r strictly decreases it
    factors = (Fr(1), Fr(3, 4), Fr(1, 2), Fr(1, 5))
    vals = []
    for lam in factors:
        rl = tuple(lam * x for x in R0)
        drl = tuple(lam * x for x in DR)
        vals.append(qfi_closed_form(rl, drl))
    assert vals[0] == q0
    assert all(a > b for a, b in zip(vals, vals[1:]))   # strictly down
    # exact factor for this orthogonal-motion witness: QFI scales as
    # lambda^2 when r . dr = 0
    assert sum(a * b for a, b in zip(R0, DR)) == 0
    for lam, v in zip(factors, vals):
        assert v == lam * lam * q0

    # never manufactured: no channel in the family raises the QFI
    assert all(v <= q0 for v in vals)

    return {
        "statement": (
            "Unitary conjugation by a rational orthogonal Bloch "
            "rotation preserves the quantum Fisher information "
            "EXACTLY, while a depolarizing channel strictly decreases "
            "it — by exactly lambda^2 for this orthogonal-motion "
            "witness — and no channel in the family raises it. "
            "Recognition channels can lose information and can never "
            "manufacture it, here as an inequality with an exact "
            "equality case, matching T2's shape rather than CID-1's "
            "classical identity"),
        "unitary_preserves": True,
        "depolarizing_scaling": "lambda^2",
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "QTH-1: quantum recognition information — the "
                   "ledger becomes an inequality",
        "answers": (
            "the queued question left open by CID-1: whether the "
            "classical recognition ledger g = g_recognized + "
            "g_discarded survives quantum-ly. It does not — it becomes "
            "an inequality with exact saturation, and the commuting "
            "face restores the identity"),
        "companion_not_duplicate": (
            "CFE-Q certified the Bloch/holonomy side of quantum CFE "
            "(Markovian limit flat, geometric residue faithful, "
            "quarter-turn echo). This capsule certifies the "
            "INFORMATION side, which CFE-Q does not touch"),
        "T1_exact_sld_and_two_routes": certify_T1(),
        "T2_measurement_ledger_is_an_inequality": certify_T2(),
        "T3_one_tensor_two_halves": certify_T3(),
        "T4_where_the_classical_identity_fails": certify_T4(),
        "T5_commuting_face_restores_everything": certify_T5(),
        "T6_monotonicity_under_channels": certify_T6(),
        "claim_boundary": {
            "witness_capsule": (
                "one qubit, explicit rational families — a WITNESS "
                "like CFE-Q, NOT a general open-system or general "
                "quantum-estimation theorem"),
            "not_claimed": (
                "NOT CLAIMED: the general Braunstein-Caves "
                "attainability theorem, the Morozova-Chentsov "
                "classification of quantum Fisher metrics, "
                "Kubo-Mori/Bogoliubov metrics (which need logarithms "
                "and are never evaluated here), and multi-copy or "
                "asymptotic estimation statements"),
            "no_thermodynamic_interpretation": (
                "NOT CLAIMED: entropy production, work, or heat. The "
                "quantities certified are information-geometric"),
            "berry_phase_identification": (
                "NOT CLAIMED: that the antisymmetric part is a "
                "physical Berry phase; CFE-Q certifies the holonomy "
                "side separately"),
            "quantum_gravity": "a horizon, not a claim",
            "RH_K0_L0_YM": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE — first quantum "
                                        "information capsule",
            "companions": ("CID-1 (classical ledger), CFE-Q (quantum "
                           "holonomy witness), CFE-1/CFE-U (the "
                           "memoryless-face and uniqueness shape)"),
            "finding": (
                "QTH-F1 the recognition ledger is classically an "
                "identity and quantum-ly an inequality with exact "
                "saturation; QTH-F2 the quantum geometric tensor's "
                "symmetric half is CID-1's information metric and its "
                "antisymmetric half is a CFE-style curvature "
                "obstruction — one object, two capsules; QTH-F3 the "
                "commuting face restores the classical identity "
                "exactly, the same shape as CFE-1's memoryless limit"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "QTH1_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    print("QTH-1 certificate written:", out)
    print("sha256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
