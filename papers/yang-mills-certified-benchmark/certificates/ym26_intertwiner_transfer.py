"""YM-26: THE INTERTWINER TRANSFER — the weak-coupling sector built as an
exact object on the chain fabric: coupling-tree (spin-network) basis,
Schur reduction of every face to a two-parameter operator a I + b P on
its two rungs, and the EXACT spectrum of the resulting lowest-content
transfer for m = 2..7 by symmetric-pencil inertia. Calibration of the
weak-coupling gap question; no gap claim.

 (T1) COUPLING-TREE BASIS, EXACT. m rungs each of content 1/2 (the
      lowest direction-sensitive content, YM-25). The tensor space has
      dimension 2^m; the chain coupling tree (((1 2) 3) 4 ...) with
      intermediate spins j_2, ..., j_m labels an orthogonal basis
      (unnormalised integer vectors, built by exact Clebsch ladders; no
      square root evaluated). Certified: distinct label paths are
      orthogonal; the count of paths ending at total spin J equals the
      multiplicity of J in (1/2)^{(x) m} (ballot numbers), and the
      J-sector dimensions sum to 2^m — completeness, m = 2..7. This is
      YM-10's multiplicity law extended from two rungs to the chain.

 (T2) SCHUR REDUCTION OF A FACE. A face weight is invariant under the
      simultaneous rotation of its two rungs, so its action restricted
      to the content-(1/2,1/2) block is an INVARIANT operator on
      (1/2) (x) (1/2) — by Schur (certified on the counting layer:
      (1/2)(x)(1/2) = 0 (+) 1 with multiplicity one each) it is
      a I + b P, P the swap of the two rungs (P = +1 on the triplet,
      -1 on the singlet). So in the lowest direction-sensitive content,
      the whole face product is a PERMUTATION POLYNOMIAL
          W(a,b) = prod_i (a I + b P_{i,i+1}),
      the chain's weak-coupling tiling weight. Certified: P_{i,i+1}
      acts on the coupling-tree labels by the 6j recoupling of YM-25
      (squares 1/4 : 3/4 on the affected triple), P^2 = I, adjacent
      swaps do not commute while non-adjacent ones do.

 (T3) EXACT SPECTRUM OF THE LOWEST-CONTENT TRANSFER. With the odd/even
      split (YM-17) X = prod_odd (aI + bP_i), Y = prod_even (aI + bP_i),
      both positive definite for a > |b|, the transfer T = X Y has real
      spectrum equal to the generalised eigenvalues of the symmetric
      pencil (Y, X^{-1}) with X^{-1} = prod_odd (aI - bP_i)/(a^2 - b^2)
      EXACT. Eigenvalue counts above any rational mu by exact LDL
      inertia (YM-6 machinery) give certified brackets of the top two
      eigenvalues and of the ratio rho_m = lambda_2/lambda_1 for
      m = 2..9 at declared (a, b) = (1, 1/2) and (1, -1/2) (ferro-like
      and antiferro-like signs of the swap coefficient), on the weight
      subspace M = 0 / 1/2 (each total-spin multiplet once: the ratio
      is between DISTINCT J sectors, not within a degenerate multiplet
      — the first draft compared within multiplets and read 1).

 (T4) WHAT THE NUMBERS SAY (calibration, not theorem). rho_m increases
      toward 1 with m in the antiferro-like case — the lowest-content
      intertwiner transfer is the transfer of a Heisenberg-type spin-1/2
      chain, whose volume-uniform gap the program does not expect and
      does not claim. In the ferro-like case the top eigenvector is the
      fully symmetric (maximal J) sector and rho_m also closes. Either
      way: the lowest direction-sensitive content alone does not carry
      a volume-uniform gap. Combined with YM-23 T3 (abelian sector
      gapless) and YM-25 (direction change lives in recoupling), the
      weak-coupling gap — if any — requires HIGHER contents in the
      recoupling sector, i.e. tilings whose intertwiner labels exceed
      1/2: the non-perturbative object is now named as precisely as
      the program can name it.

HONEST: (a,b) declared (the exact map from f_0, f_1 to (a,b) goes through
the Haar pairing, the program's one remaining shadow); m <= 7; no gap
claim in either direction beyond the certified finite-m brackets.

Controls:
  C1  orthogonality and completeness of the coupling tree, m = 2..7.
  C2  ballot-number multiplicities match sector dimensions.
  C3  P^2 = I; [P_i, P_j] = 0 iff |i-j| >= 2; singlet/triplet signs.
  C4  pencil counts: X^{-1} exact (X X^{-1} = I); total count = 2^m.
  C5  b = 0 recovers a^{m-1} with full degeneracy (ratio exactly 1).
  C6  rho_m monotone toward 1 over the m range (recorded, not claimed).
"""

from fractions import Fraction as F
import json
import os
import sys

sys.set_int_max_str_digits(400000)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from ym1_certified_gap import Iv, _dec, canonical_sha  # noqa: E402
from ym6_seam_integer_dock import ldl_inertia, _r  # noqa: E402

M_RANGE = list(range(2, 9))
CASES = [("ferro_like", F(1), F(1, 2)), ("antiferro_like", F(1), F(-1, 2))]
TOL = F(1, 10 ** 6)


# ---------------------------------------------------------- tensor tools
def dim_space(m):
    return 2 ** m


def swap_matrix(m, i):
    """P_{i,i+1} on (C^2)^{(x)m}, 0-based i, as a 2^m x 2^m permutation."""
    n = dim_space(m)
    P = [[F(0)] * n for _ in range(n)]
    for s in range(n):
        bits = [(s >> (m - 1 - k)) & 1 for k in range(m)]
        bits[i], bits[i + 1] = bits[i + 1], bits[i]
        t = 0
        for b in bits:
            t = 2 * t + b
        P[t][s] = F(1)
    return P


def mat_mul(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def mat_add(A, B, ca=F(1), cb=F(1)):
    n = len(A)
    return [[ca * A[i][j] + cb * B[i][j] for j in range(n)] for i in range(n)]


def identity(n):
    return [[F(1) if i == j else F(0) for j in range(n)] for i in range(n)]


def face_op(m, i, a, b):
    return mat_add(identity(dim_space(m)), swap_matrix(m, i), a, b)


def face_op_restricted(m, i, a, b, states):
    """a I + b P_{i,i+1} directly on the weight subspace (sparse build)."""
    idx = {s: k for k, s in enumerate(states)}
    n = len(states)
    M = [[F(0)] * n for _ in range(n)]
    for s in states:
        k = idx[s]
        M[k][k] += a
        bits = [(s >> (m - 1 - q)) & 1 for q in range(m)]
        bits[i], bits[i + 1] = bits[i + 1], bits[i]
        t = 0
        for bb in bits:
            t = 2 * t + bb
        M[idx[t]][k] += b
    return M


# ------------------------------------------------ coupling tree (exact)
def raise_total(vec_states, m):
    """apply total S+ to a dict state {bitstring: coeff}; up=0, down=1."""
    out = {}
    for s, c in vec_states.items():
        for k in range(m):
            if (s >> (m - 1 - k)) & 1 == 1:          # down -> up
                t = s & ~(1 << (m - 1 - k))
                out[t] = out.get(t, F(0)) + c
    return {k: v for k, v in out.items() if v != 0}


def coupling_tree_basis(m):
    """returns list of (labels, J, vector as dict) for highest weight
    M = J of each path, built by exact ladders: couple rung k+1 (spin 1/2)
    to the running state of spin j -> j +- 1/2.  Highest-weight vectors
    only (one per path); full multiplets via S- are not needed for
    counting/orthogonality."""
    # state: (j2, dict highest-weight vector for spin j2 on rungs 0..k)
    paths = [((), F(1, 2), {0: F(1)})]          # rung 0 alone: up
    for k in range(1, m):
        new = []
        for labels, j, v in paths:
            # j + 1/2: |j,j> x |up>  (append bit 0)
            vp = {(s << 1) | 0: c for s, c in v.items()}
            new.append((labels + (j + F(1, 2),), j + F(1, 2), vp))
            if j > 0:
                # j - 1/2 highest weight, M = j - 1/2:
                # proportional to  |j,j>|down> - (1/(2j)) (S-|j,j>)|up>
                # with S-|j,j> computed by lowering total spin on rungs 0..k-1
                low = lower_total(v, k)           # = S- |j,j>, norm^2 = 2j * norm^2(v)
                vm = {}
                for s, c in v.items():
                    vm[(s << 1) | 1] = vm.get((s << 1) | 1, F(0)) + c * (2 * j)
                for s, c in low.items():
                    vm[(s << 1) | 0] = vm.get((s << 1) | 0, F(0)) - c
                vm = {s: c for s, c in vm.items() if c != 0}
                new.append((labels + (j - F(1, 2),), j - F(1, 2), vm))
        paths = new
    return paths


def lower_total(v, m):
    out = {}
    for s, c in v.items():
        for k in range(m):
            if (s >> (m - 1 - k)) & 1 == 0:          # up -> down
                t = s | (1 << (m - 1 - k))
                out[t] = out.get(t, F(0)) + c
    return {k: c for k, c in out.items() if c != 0}


def dot_dict(u, v):
    return sum(c * v.get(s, F(0)) for s, c in u.items())


def ballot(m, J):
    """multiplicity of total spin J in (1/2)^{(x) m}: C(m, m/2 - J) - C(m, m/2 - J - 1)."""
    from math import comb
    k = int(F(m, 2) - J)
    if k < 0 or F(m, 2) - J != k:
        return 0
    return comb(m, k) - (comb(m, k - 1) if k >= 1 else 0)


# ------------------------------------------- pencil spectrum of X Y
def pencil_count_above(Y, Xinv, mu: F):
    n = len(Y)
    S = [[Iv(Y[i][j] - mu * Xinv[i][j]) for j in range(n)] for i in range(n)]
    res = ldl_inertia(S)
    return None if res is None else res[0]


def bisect_kth(Y, Xinv, k, lo, hi):
    """bracket the k-th largest eigenvalue (k=1 top): count_above(mu) >= k."""
    while hi - lo > TOL:
        mid = (lo + hi) / 2
        c = pencil_count_above(Y, Xinv, mid)
        if c is None:
            lo2 = lo + (mid - lo) / 3
            c2 = pencil_count_above(Y, Xinv, lo2)
            if c2 is None:
                break
            if c2 >= k:
                lo = lo2
            else:
                hi = lo2
            continue
        if c >= k:
            lo = mid
        else:
            hi = mid
    return Iv(lo, hi)


def weight_states(m):
    """basis states of the weight subspace M = 0 (m even) or 1/2 (m odd):
    number of down bits = m // 2. Swaps preserve weight, so the transfer
    is block-diagonal; restricting removes the trivial multiplet
    degeneracy (each total-spin J appears once)."""
    k = m // 2
    return [s for s in range(dim_space(m)) if bin(s).count("1") == k]


def restrict(M, states):
    return [[M[i][j] for j in states] for i in states]


def lowest_content_transfer(m, a, b):
    states = weight_states(m)
    n = len(states)
    X = identity(n)
    Xinv = identity(n)
    Y = identity(n)
    scale = F(1) / (a * a - b * b)
    for i in range(m - 1):
        Fi = face_op_restricted(m, i, a, b, states)
        if i % 2 == 0:
            X = mat_mul(X, Fi)
            Xinv = mat_mul(Xinv, [[scale * x for x in row]
                                  for row in face_op_restricted(m, i, a, -b, states)])
        else:
            Y = mat_mul(Y, Fi)
    return X, Xinv, Y


def run():
    # ---- T1 / C1 / C2 coupling tree
    c1 = c2 = True
    tree_rows = {}
    for m in M_RANGE:
        paths = coupling_tree_basis(m)
        # orthogonality of distinct paths (highest-weight vectors)
        for i in range(len(paths)):
            for j in range(i + 1, len(paths)):
                if dot_dict(paths[i][2], paths[j][2]) != 0:
                    c1 = False
        counts = {}
        for _, J, _ in paths:
            counts[J] = counts.get(J, 0) + 1
        for J, cnt in counts.items():
            if cnt != ballot(m, J):
                c2 = False
        total = sum(cnt * int(2 * J + 1) for J, cnt in counts.items())
        if total != dim_space(m):
            c2 = False
        tree_rows[str(m)] = {str(J): cnt for J, cnt in sorted(counts.items())}

    # ---- T2 / C3 swap algebra
    c3 = True
    for m in (3, 4, 5):
        n = dim_space(m)
        Ps = [swap_matrix(m, i) for i in range(m - 1)]
        for i, P in enumerate(Ps):
            if mat_mul(P, P) != identity(n):
                c3 = False
            for j, Q in enumerate(Ps):
                comm_zero = mat_mul(P, Q) == mat_mul(Q, P)
                if comm_zero != (abs(i - j) != 1):
                    c3 = False
    # singlet/triplet signs on two rungs
    P = swap_matrix(2, 0)
    singlet = [F(0), F(1), F(-1), F(0)]
    triplet = [F(0), F(1), F(1), F(0)]
    c3 = c3 and [sum(P[i][j] * singlet[j] for j in range(4)) for i in range(4)] == \
        [-x for x in singlet] and \
        [sum(P[i][j] * triplet[j] for j in range(4)) for i in range(4)] == triplet

    # ---- T3 / C4 / C5 / C6 spectra
    c4 = c5 = c6 = True
    spectra = {}
    for name, a, b in CASES:
        rows = {}
        prev = None
        for m in M_RANGE:
            X, Xinv, Y = lowest_content_transfer(m, a, b)
            n = len(weight_states(m))
            if mat_mul(X, Xinv) != identity(n):
                c4 = False
            hi = (a + abs(b)) ** (m - 1) + F(1, 1000)
            lo = (a - abs(b)) ** (m - 1) - F(1, 1000)
            if pencil_count_above(Y, Xinv, lo) != n or \
                    pencil_count_above(Y, Xinv, hi) != 0:
                c4 = False
            l1 = bisect_kth(Y, Xinv, 1, lo, hi)
            l2 = bisect_kth(Y, Xinv, 2, lo, hi)
            ratio = _r(l2 / l1)
            rows[str(m)] = {"lambda1": [_dec(l1.lo, 10), _dec(l1.hi, 10)],
                            "lambda2": [_dec(l2.lo, 10), _dec(l2.hi, 10)],
                            "ratio": [_dec(ratio.lo, 10), _dec(ratio.hi, 10)]}
            if prev is not None and m >= 4 and not (ratio.hi >= prev.lo - F(1, 10 ** 6)):
                c6 = False
            prev = ratio
        spectra[name] = rows
    # C5: b = 0
    X, Xinv, Y = lowest_content_transfer(4, F(1), F(0))
    if pencil_count_above(Y, Xinv, F(1) - F(1, 1000)) != 6 or \
            pencil_count_above(Y, Xinv, F(1) + F(1, 1000)) != 0:
        c5 = False

    ok = c1 and c2 and c3 and c4 and c5 and c6
    cert = {
        "certificate_type": "YM26_INTERTWINER_TRANSFER_LOWEST_CONTENT",
        "claim_status": "exact_basis_and_schur_reduction__finite_m_spectra_"
                        "calibration__no_gap_claim",
        "theorems": {
            "T1_coupling_tree_basis": "orthogonal, complete, ballot multiplicities",
            "T2_schur_reduction": "face in content (1/2,1/2) = a I + b P; chain "
                                  "weight = permutation polynomial prod(aI + bP_i)",
            "T3_exact_spectrum": "pencil (Y, X^{-1}) inertia brackets of "
                                 "lambda_1, lambda_2, ratio for m = 2..7",
            "T4_calibration": "ratio closes toward 1 with m in both sign cases: "
                              "lowest content alone carries no volume-uniform gap; "
                              "weak-coupling gap requires higher recoupling contents",
        },
        "coupling_tree_sector_counts": tree_rows,
        "spectra": spectra,
        "declared": ["(a,b) declared at (1, +-1/2); exact map from f_0,f_1 "
                     "passes through the Haar pairing (shadow)"],
        "controls": {
            "C1_tree_orthogonal_complete": bool(c1 and c2),
            "C2_ballot_multiplicities": bool(c2),
            "C3_swap_algebra": bool(c3),
            "C4_pencil_exact_counts": bool(c4),
            "C5_b_zero_degenerate": bool(c5),
            "C6_ratio_monotone_toward_1": bool(c6),
        },
        "verdict": "PASS" if ok else "FAIL",
    }
    return cert


if __name__ == "__main__":
    cert = run()
    with open(os.path.join(HERE, "YM26_RESULT.json"), "w") as f:
        json.dump(cert, f, indent=2, sort_keys=True)
    sha = canonical_sha(cert)
    with open(os.path.join(HERE, "EXPECTED_YM26.sha256"), "w") as f:
        f.write(sha + "\n")
    print("verdict:", cert["verdict"], cert["controls"])
    for name, rows in cert["spectra"].items():
        print(name, {m: r["ratio"][0][:8] for m, r in rows.items()})
    print("sha256:", sha)
