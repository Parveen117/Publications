"""RST-2: recognition-seam topology, the arithmetic and fusion spine.

Source: "Recognition-Seam Topology and Daggered Index Geometry" (M. Dabas,
July 2026, arXiv source v1), sections 8-14 — the FINITE, EXACT content.
Companion to RST-1 (sections 2-7). Everything here is exact rational or
integer arithmetic; every logarithm is a FORMAL symbol (a vector of prime
exponents), every square root is avoided by working with squares and signs,
and every infinite/analytic statement of the source is declared, cited and
NOT consumed.

BLOCKS

  T1  SAMPLED REFLECTION AND SPECTRAL BLINDNESS (8.2-8.3). The pair-swap
      J_N on exactly reflected rational sample points s_k^+- = 1/2 +- eps
      + i t_k, its symmetric/antisymmetric eigenvectors, the leakage
      criterion, and the finite blindness witness: the seam projection of
      a vector is independent of its anti-seam component.

  T2  WEIGHTS, LEAKAGE NORM, CROSSING VS AVOIDED (8.4, 8.7 + the finite
      content of 8.5). Reflection-compatible weight criterion both ways
      (equal pair weights <=> the swap is isometric, certified with
      squared norms). The pure leakage operator on a finite section has
      squared norm EXACTLY max a_k^2 (attained), and an unbounded profile
      has finite-section norms growing without bound — an exact growth
      witness; the infinite-carrier statement is declared. Exact vs
      avoided crossings: det D_{c,mu}(eta) = -((eta-c)^2 + mu^2) <= -mu^2
      < 0, so zero is never in the spectrum of the avoided family, while
      oriented scalar branches cross with their declared signs.

  T3  ARITHMETIC RECOGNITION (10). THE FORMAL-LOG MOVE: log n is the
      vector of prime exponents of n in the free abelian group on the
      primes, so Moebius inversion (L1)*mu = Lambda is certified as an
      exact identity of integer vectors for every n <= 60 — no logarithm
      is ever evaluated. Radical recognition: the Eye transfer of Lambda
      lands ENTIRELY in prime-radical classes (anti-seam current zero),
      the raw log clock leaks into class 6, the controlled perturbation
      Lambda_supp + eps delta_6 produces current EXACTLY eps e_6, the
      crossing pairing reverses sign with orientation, and the
      topological lift preserves both, giving the exact chain
      Lambda -> 0 -> 0 and contaminated -> eps e_6 -> sgn(eps).

  T4  PRIME-SHELL TRANSFER AND DEFECT CURRENTS (11). Exact finite
      covariance U_lam A = A P_ar, J A = A, P_- A = 0, J U* J = U for
      rational shell weights (the Gaussian damping FORM is declared; the
      theorems hold for any symmetric positive weights). Asymmetric
      transport gives the exact leakage formula P_-(Delta) e_p =
      ((w_R - w_L)/2)(e_R - e_L), vanishing iff the reflected weights
      agree shell by shell. Exact transfer => zero defect current, and
      the CONVERSE FAILS: a common seam-sector error has zero anti-seam
      current — certified witness. Aligned/reversed/avoided crossing
      pairings are +alpha, -alpha, 0 exactly; the damping-stable
      obstruction keeps its sign with nonincreasing magnitude; and the
      local reflected displacement prototype has the exact trichotomy
      eps = 0 / > 0 / < 0  =>  O = 0 / > 0 / < 0.

  T5  ZERO-IDENTIFICATION BOUNDARY (13). The finite paired determinant
      D_N(E) = prod (1 - E^2/lambda_j^2): D_N(0) = 1, evenness as a
      statement about coefficients, zeros exactly +-lambda_j WITH
      multiplicity (counted by exact polynomial division), and the
      logarithmic-derivative identity certified in cross-multiplied
      polynomial form (no division, no evaluation gaps). The
      nonvanishing-factor lemma on polynomial models: multiplying by a
      unit preserves the zero divisor, with multiplicity. A working
      obstruction-class validator: it ACCEPTS a clean synthetic ledger
      and DETECTS extra zeros (Z1), missing zeros (Z2) and multiplicity
      mismatch (Z3) — and the capsule states that passing the validator
      provides NO evidence that any synthetic set equals Z_xi.

  T6  SURYA/WEIL FUSION ON RATIONAL MODELS (14). The Feshbach congruence
      M = S^T diag(A, F) S with S = [[I, A^{-1}C^T],[0, I]] as an exact
      matrix identity, giving W_- >= 0 <=> F >= 0 when A > 0 — both
      directions certified with principal minors and explicit negative
      witnesses transported through the congruence. The rectified Surya
      correction handled WITHOUT square roots: the squared amplitude
      ssq(x) = x^2/(1+x^2) for x < 0 and 0 otherwise satisfies
      ssq(x) = 0 <=> x >= 0, so on rational-eigenvalue models
      A_SS = 0 <=> A_SS^2 = 0 <=> W_0 >= 0 — the two terminal signs are
      one obstruction. Seam charge: k_Sigma = #{relative eigenvalues
      > 1} is additive across the even/odd direct sum, equals MINUS the
      spectral flow of I - t K (every crossing t = 1/mu is downward,
      slope -mu < 0: no cancellation is possible), and W_eta =
      A^{1/2}(I - K)A^{1/2} with rational square roots certifies
      k_Sigma = 0 <=> W_eta >= 0. The Surya angle indicator equality
      1_{(0,infty)}(c_SS(K)) = 1_{(1,infty)}(K) on rational spectra:
      four presentations, one obstruction. Principal holonomy is
      insufficient: m mod 1 = 0 for EVERY integer m while the lift m
      itself need not vanish — the same projection blindness certified
      in UGD-1 T4.

DECLARED, NOT CERTIFIED: Thm 8.8 (entire trace-class damped family),
Thm 9.1 (undamped strip classification), Thm 8.5 (Schatten thresholds,
infinite part), Thm 9.4 (archimedean defect), the Gaussian damping form
e^{-lambda (log p)^2}, the 1/sqrt(2) transfer normalization (a positive
constant affecting no vanishing statement), the completed Weil operator
W_3 and the Weil-to-RH equivalence (imported in the source from the MP
X3 stack), the odd compact reference (open gate F1), and every
completed-carrier statement. The open fusion gates F1-F5 remain open.
Riemann hypothesis: ABSTAIN.
"""

import hashlib
import json
import os
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))

# ----------------------------------------------------------------------
# exact linear algebra helpers
# ----------------------------------------------------------------------


def meye(n):
    return [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def mzero(n, m=None):
    m = n if m is None else m
    return [[Fr(0)] * m for _ in range(n)]


def madd(A, B):
    return [[A[i][j] + B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def msub(A, B):
    return [[A[i][j] - B[i][j] for j in range(len(A[0]))]
            for i in range(len(A))]


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
            for i in range(n)]


def mT(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def mscale(c, A):
    return [[c * x for x in row] for row in A]


def mv(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x))) for i in range(len(A))]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def det(A):
    n = len(A)
    if n == 1:
        return A[0][0]
    out = Fr(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        out += (-1) ** j * A[0][j] * det(minor)
    return out


def inv(A):
    n = len(A)
    d = det(A)
    assert d != 0
    out = mzero(n)
    for i in range(n):
        for j in range(n):
            minor = [row[:i] + row[i + 1:]
                     for k, row in enumerate(A) if k != j]
            out[i][j] = (-1) ** (i + j) * (det(minor) if n > 1 else Fr(1)) / d
    return out


def principal_minors(A):
    """All principal minors of a square matrix, exactly."""
    from itertools import combinations
    n = len(A)
    out = []
    for r in range(1, n + 1):
        for idx in combinations(range(n), r):
            sub = [[A[i][j] for j in idx] for i in idx]
            out.append(det(sub))
    return out


def is_psd(A):
    """Symmetric rational matrix PSD <=> every principal minor >= 0."""
    assert A == mT(A)
    return all(m >= 0 for m in principal_minors(A))


def block2(A, B, C, D):
    top = [ra + rb for ra, rb in zip(A, B)]
    bot = [rc + rd for rc, rd in zip(C, D)]
    return top + bot


# ----------------------------------------------------------------------
# exact polynomials over Q
# ----------------------------------------------------------------------


def ptrim(p):
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def padd(p, q):
    m = max(len(p), len(q))
    return ptrim([(p[i] if i < len(p) else Fr(0))
                  + (q[i] if i < len(q) else Fr(0)) for i in range(m)])


def pmul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return ptrim(out)


def pscale(c, p):
    return ptrim([c * a for a in p])


def pderiv(p):
    return ptrim([Fr(i) * p[i] for i in range(1, len(p))]) or [Fr(0)]


def peval(p, x):
    out = Fr(0)
    for c in reversed(p):
        out = out * x + c
    return out


def pdivmod(p, q):
    p, q = ptrim(p), ptrim(q)
    assert q != [Fr(0)]
    quot = [Fr(0)] * max(1, len(p) - len(q) + 1)
    rem = list(p)
    while len(rem) >= len(q) and rem != [Fr(0)]:
        c = rem[-1] / q[-1]
        d = len(rem) - len(q)
        quot[d] = c
        rem = [rem[i] - (c * q[i - d] if 0 <= i - d < len(q) else Fr(0))
               for i in range(len(rem))]
        assert rem[-1] == 0
        rem.pop()
        rem = rem or [Fr(0)]
    return ptrim(quot), ptrim(rem)


def root_multiplicity(p, r):
    m = 0
    while peval(p, r) == 0:
        p, rem = pdivmod(p, [-r, Fr(1)])
        assert rem == [Fr(0)]
        m += 1
    return m


# ----------------------------------------------------------------------
# T1  sampled reflection and spectral blindness (8.2-8.3)
# ----------------------------------------------------------------------

EPS = Fr(1, 5)
HEIGHTS = (Fr(0), Fr(2), Fr(-7, 3))


def certify_T1():
    # exactly reflected sample points: rho(s) = 1 - s swaps them
    pairs = []
    for t in HEIGHTS:
        sp = (Fr(1, 2) + EPS, t)
        sm = (Fr(1, 2) - EPS, t)
        assert (Fr(1) - sp[0], sp[1]) == sm
        assert (Fr(1) - sm[0], sm[1]) == sp
        pairs.append((sp, sm))

    # pair-swap on one reflected pair: basis (e+, e-)
    J = [[Fr(0), Fr(1)], [Fr(1), Fr(0)]]
    I2 = meye(2)
    assert mmul(J, J) == I2 and mT(J) == J
    u = [Fr(1), Fr(1)]                       # symmetric (unnormalized)
    v = [Fr(1), Fr(-1)]                      # antisymmetric
    assert mv(J, u) == u
    assert mv(J, v) == [-x for x in v]
    assert dot(u, v) == 0

    Pp = mscale(Fr(1, 2), madd(I2, J))
    Pm = mscale(Fr(1, 2), msub(I2, J))

    # leakage criterion on the 2x2 grid
    vals = (Fr(-1), Fr(0), Fr(1), Fr(2))
    checked = 0
    for a in vals:
        for b in vals:
            for c in vals:
                for d in vals:
                    A = [[a, b], [c, d]]
                    commutes = mmul(A, J) == mmul(J, A)
                    leaks0 = (mmul(Pm, mmul(A, Pp)) == mzero(2)
                              and mmul(Pp, mmul(A, Pm)) == mzero(2))
                    assert commutes == leaks0
                    checked += 1

    # spectral blindness: seam projection independent of anti-seam part
    for c1 in (Fr(3), Fr(-2)):
        base = mv(Pp, [c1, Fr(0)])
        for c2 in (Fr(0), Fr(1), Fr(-5), Fr(7, 2)):
            vec = [c1 * u[0] + c2 * v[0], c1 * u[1] + c2 * v[1]]
            proj = mv(Pp, vec)
            assert proj == [c1 * u[0], c1 * u[1]]
        assert base == mv(Pp, [c1, Fr(0)])

    return {
        "statement": (
            "Reflected sample points s = 1/2 +- eps + i t are exchanged "
            "EXACTLY by rho(s) = 1 - s; the pair swap J_N is a "
            "self-adjoint involution with symmetric/antisymmetric "
            "eigenvectors; leakage vanishes iff [A, J] = 0 on the full "
            "256-matrix grid; and the seam projection of a vector is "
            "independent of its anti-seam component — visible seam data "
            "alone never certify absence of an anti-seam ledger. Model "
            "points only: no assertion that any sampled point is a zero "
            "of xi"),
        "eps": str(EPS),
        "heights": [str(t) for t in HEIGHTS],
        "grid_checked": checked,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  weights, leakage norm, crossing vs avoided (8.4, 8.7)
# ----------------------------------------------------------------------


def wnorm2(x, w):
    return sum(wk * xk * xk for wk, xk in zip(w, x))


def certify_T2():
    # weight criterion on two reflected pairs: coords (x1+,x1-,x2+,x2-)
    def swap(x):
        return [x[1], x[0], x[3], x[2]]

    grid = [[Fr(a), Fr(b), Fr(c), Fr(d)]
            for a in (-1, 0, 2) for b in (-1, 0, 2)
            for c in (0, 1) for d in (0, 3)]

    w_eq = [Fr(2), Fr(2), Fr(1, 3), Fr(1, 3)]
    assert all(wnorm2(swap(x), w_eq) == wnorm2(x, w_eq) for x in grid)

    w_neq = [Fr(2), Fr(1), Fr(1, 3), Fr(1, 3)]
    witness = [Fr(1), Fr(0), Fr(0), Fr(0)]
    assert wnorm2(swap(witness), w_neq) != wnorm2(witness, w_neq)

    # pure leakage operator on a finite section: A_a u_k = a_k v_k.
    # squared operator norm = max a_k^2, attained
    a = (Fr(3), Fr(-5), Fr(1, 2))
    w = (Fr(1), Fr(2), Fr(4))
    amax2 = max(ak * ak for ak in a)
    for x in ([Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(1), Fr(0)],
              [Fr(1), Fr(1), Fr(1)], [Fr(-2), Fr(3), Fr(5)]):
        img2 = sum(wk * (ak * xk) ** 2 for wk, ak, xk in zip(w, a, x))
        assert img2 <= amax2 * wnorm2(x, w)
    k_star = max(range(3), key=lambda k: a[k] * a[k])
    e = [Fr(0)] * 3
    e[k_star] = Fr(1)
    img2 = sum(wk * (ak * xk) ** 2 for wk, ak, xk in zip(w, a, e))
    assert img2 == amax2 * wnorm2(e, w)          # attained exactly

    # unbounded profile: finite-section norms grow without bound
    growth = []
    for N in (1, 4, 16, 64):
        profile_max2 = max(Fr(k) ** 2 for k in range(1, N + 1))
        growth.append((N, profile_max2))
    assert all(g2 == Fr(N) ** 2 for N, g2 in growth)
    assert [g2 for _, g2 in growth] == sorted(
        {g2 for _, g2 in growth})                # strictly increasing

    # exact vs avoided crossings (8.7), no square roots:
    # avoided family never has 0 in its spectrum
    for c in (Fr(0), Fr(1, 2)):
        for mu in (Fr(1, 3), Fr(2)):
            for eta in (Fr(-1), Fr(0), c, c + Fr(1, 7), Fr(3)):
                D = [[eta - c, mu], [mu, -(eta - c)]]
                d = det(D)
                assert d == -((eta - c) ** 2 + mu ** 2)
                assert d <= -mu ** 2 < 0         # gap^2 >= mu^2
    # oriented scalar branches cross with declared signs
    for eps_j in (1, -1):
        c = Fr(1, 2)
        before = eps_j * (Fr(0) - c)
        after = eps_j * (Fr(1) - c)
        assert (before < 0 < after) == (eps_j == 1)
        assert (after < 0 < before) == (eps_j == -1)

    return {
        "statement": (
            "Equal reflected pair weights are EXACTLY the condition for "
            "the swap to be isometric (both directions, squared norms); "
            "the finite pure-leakage operator has squared norm exactly "
            "max a_k^2, attained; an unbounded profile has strictly "
            "growing finite-section norms (exact growth witness — the "
            "infinite-carrier boundedness criterion is declared); the "
            "avoided two-level family has det = -((eta-c)^2 + mu^2) <= "
            "-mu^2 < 0 so zero is NEVER in its spectrum, while oriented "
            "scalar branches cross with their declared signs: index "
            "changes need true crossings, not close approaches"),
        "leakage_norm_squared": str(max(ak * ak for ak in a)),
        "growth_witness": [[N, str(g2)] for N, g2 in growth],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T3  arithmetic recognition with formal logarithms (10)
# ----------------------------------------------------------------------

N_ARITH = 60


def factorize(n):
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def logvec(n):
    """log n as a FORMAL vector of prime exponents (never evaluated)."""
    return factorize(n)


def vadd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
        if out[k] == 0:
            del out[k]
    return out


def vscale(c, a):
    return {k: c * v for k, v in a.items() if c * v != 0}


def mobius(n):
    f = factorize(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def Lambda_vec(n):
    f = factorize(n)
    if len(f) == 1:
        p = next(iter(f))
        return {p: 1}                          # log p, formally
    return {}


def rad(n):
    out = 1
    for p in factorize(n):
        out *= p
    return out if n > 1 else 1


def is_prime(n):
    return n > 1 and len(factorize(n)) == 1 and factorize(n)[
        next(iter(factorize(n)))] == 1 and rad(n) == n and \
        list(factorize(n).values()) == [1]


def eye_transfer(stream):
    """T_E f (q) = sum over n with rad(n) = q of f(n); f vector-valued."""
    out = {}
    for n, fv in stream.items():
        q = rad(n)
        out[q] = vadd(out.get(q, {}), fv)
    return {q: v for q, v in out.items() if v}


def eye_transfer_scalar(stream):
    out = {}
    for n, fv in stream.items():
        q = rad(n)
        out[q] = out.get(q, Fr(0)) + fv
    return {q: v for q, v in out.items() if v != 0}


def anti_seam_classes(N):
    return sorted(q for q in {rad(n) for n in range(2, N + 1)}
                  if not is_prime(q) and q > 1)


def certify_T3():
    # Moebius pressure extractor: (L1)*mu = Lambda, formally, all n
    for n in range(1, N_ARITH + 1):
        conv = {}
        for d in range(1, n + 1):
            if n % d == 0:
                conv = vadd(conv, vscale(mobius(n // d), logvec(d)))
        assert conv == Lambda_vec(n), n
    # and the classical source identity log n = sum_{d|n} Lambda(d)
    for n in range(1, N_ARITH + 1):
        s = {}
        for d in range(1, n + 1):
            if n % d == 0:
                s = vadd(s, Lambda_vec(d))
        assert s == logvec(n)

    # radical recognition: Lambda transfers entirely to prime classes
    lam_stream = {n: Lambda_vec(n) for n in range(2, N_ARITH + 1)}
    t_lam = eye_transfer(lam_stream)
    assert all(is_prime(q) for q in t_lam)
    anti = anti_seam_classes(N_ARITH)
    assert 6 in anti
    assert all(q not in t_lam for q in anti)     # anti-seam current 0

    # the raw log clock leaks: class 6 receives log 6 + log 12 + ...
    log_stream = {n: logvec(n) for n in range(2, N_ARITH + 1)}
    t_log = eye_transfer(log_stream)
    assert 6 in t_log and t_log[6] != {}
    contributors_6 = [n for n in range(2, N_ARITH + 1) if rad(n) == 6]
    assert contributors_6[0] == 6

    # controlled contamination: Lambda_supp + eps delta_6 -> eps e_6
    eps = Fr(3, 7)
    supp_stream = {n: (Fr(1) if Lambda_vec(n) else Fr(0))
                   for n in range(2, N_ARITH + 1)}
    clean = eye_transfer_scalar(supp_stream)
    assert all(is_prime(q) for q in clean)
    contaminated = dict(supp_stream)
    contaminated[6] = contaminated.get(6, Fr(0)) + eps
    t_cont = eye_transfer_scalar(contaminated)
    anti_current = {q: v for q, v in t_cont.items() if q in anti}
    assert anti_current == {6: eps}              # EXACTLY eps e_6

    # crossing pairing and orientation reversal; topological lift
    pair_plus = anti_current.get(6, Fr(0)) * Fr(1)
    pair_minus = anti_current.get(6, Fr(0)) * Fr(-1)
    assert pair_plus == eps and pair_minus == -eps
    lift = (Fr(0), anti_current.get(6, Fr(0)))   # (seam, anti-seam)
    assert lift[1] == eps and lift[0] == 0       # sign preserved

    return {
        "statement": (
            "With log n carried as the FORMAL vector of prime exponents "
            "(never evaluated), Moebius inversion (L1)*mu = Lambda and "
            "log n = sum_(d|n) Lambda(d) are exact integer-vector "
            "identities for all n <= 60. Radical recognition sends the "
            "transferred Lambda stream ENTIRELY into prime-radical seam "
            "classes (anti-seam current identically zero) while the raw "
            "log clock leaks into class 6; the controlled perturbation "
            "produces anti-seam current EXACTLY eps e_6, the crossing "
            "pairing reverses sign with orientation, and the isometric "
            "topological lift preserves both — the exact chain "
            "Lambda -> 0 -> 0 versus contaminated -> eps e_6 -> sgn(eps). "
            "A pressure stream is clean only when its REJECTED component "
            "is measured. The finite sign is a witness, NOT a Fredholm "
            "index (paper Remark 11.8)"),
        "N": N_ARITH,
        "eps": str(eps),
        "anti_seam_classes_below_N": anti[:8],
        "class_6_contributors": contributors_6,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  prime-shell transfer and defect currents (11)
# ----------------------------------------------------------------------

SHELLS = (2, 3, 5, 7, 11)


def shell_ops(weights_R, weights_L, weights_ar):
    """Matrices on basis (e_{p,R}, e_{p,L})_p; transfer from (e_p)_p.

    A is the UNNORMALIZED transfer e_p -> e_{p,R} + e_{p,L}; the 1/sqrt2
    normalization is a declared positive constant affecting no vanishing
    statement.
    """
    n = len(SHELLS)
    A = mzero(2 * n, n)
    U = mzero(2 * n)
    Par = mzero(n)
    J = mzero(2 * n)
    for i, p in enumerate(SHELLS):
        A[2 * i][i] = Fr(1)
        A[2 * i + 1][i] = Fr(1)
        U[2 * i][2 * i] = weights_R[p]
        U[2 * i + 1][2 * i + 1] = weights_L[p]
        Par[i][i] = weights_ar[p]
        J[2 * i][2 * i + 1] = Fr(1)
        J[2 * i + 1][2 * i] = Fr(1)
    return A, U, Par, J


def certify_T4():
    n = len(SHELLS)
    w = {2: Fr(1, 2), 3: Fr(1, 3), 5: Fr(1, 7), 7: Fr(2, 9), 11: Fr(1, 13)}
    A, U, Par, J = shell_ops(w, w, w)
    I2n = meye(2 * n)
    Pm = mscale(Fr(1, 2), msub(I2n, J))

    # Theorem 11.2, exactly
    assert mmul(U, A) == mmul(A, Par)            # covariance
    assert mmul(J, A) == A                       # image in seam sector
    assert mmul(Pm, A) == mzero(2 * n, n)        # zero anti-seam part
    assert mmul(J, mmul(mT(U), J)) == U          # dagger covariance

    # asymmetric transport: exact leakage formula (Prop 11.3)
    wR = dict(w)
    wL = dict(w)
    wR[5], wL[5] = Fr(1, 4), Fr(1, 6)            # mismatch on shell 5
    _, Uasym, _, _ = shell_ops(wR, wL, w)
    Delta = msub(mmul(Uasym, A), mmul(A, Par))
    PmD = mmul(Pm, Delta)
    for i, p in enumerate(SHELLS):
        col = [PmD[r][i] for r in range(2 * n)]
        expect = [Fr(0)] * (2 * n)
        coef = (wR[p] - wL[p]) / 2
        expect[2 * i] = coef
        expect[2 * i + 1] = -coef
        assert col == expect
        assert (coef == 0) == (wR[p] == wL[p])

    # exact transfer => zero current; CONVERSE FAILS witness:
    # a common seam-sector error has zero anti-seam current
    w_common = dict(w)
    w_common[3] = Fr(1, 2)                       # both sectors wrong, equally
    _, Ucom, _, _ = shell_ops(w_common, w_common, w)
    Delta_c = msub(mmul(Ucom, A), mmul(A, Par))
    assert Delta_c != mzero(2 * n, n)            # routes differ
    assert mmul(Pm, Delta_c) == mzero(2 * n, n)  # but J_def = 0

    # crossing pairings (11.6-11.7): aligned/reversed/avoided
    alpha = Fr(-5, 3)
    c_minus = [Fr(1), Fr(-1)]                    # anti-seam direction
    J_def = [alpha * x for x in c_minus]
    norm2 = dot(c_minus, c_minus)
    assert dot(J_def, c_minus) == alpha * norm2
    assert dot(J_def, [-x for x in c_minus]) == -alpha * norm2
    c_avoided = [Fr(4), Fr(4)]                   # seam sector
    assert dot(J_def, c_avoided) == 0

    # damping-stable obstruction: sign constant, magnitude nonincreasing
    g = (Fr(1), Fr(2, 3), Fr(1, 2), Fr(1, 5))    # nonincreasing profile
    eps = Fr(2, 7)
    Os = [eps * gk * norm2 for gk in g]
    assert all(o > 0 for o in Os)
    assert all(a >= b for a, b in zip(Os, Os[1:]))
    # trichotomy of the displacement prototype (Thm 11.10)
    for e, expect in ((Fr(0), 0), (Fr(1, 9), 1), (Fr(-4), -1)):
        O = e * g[0] * norm2
        sgn = 0 if O == 0 else (1 if O > 0 else -1)
        assert sgn == expect

    return {
        "statement": (
            "Exact finite prime-shell transfer: U A = A P_ar, J A = A, "
            "P_- A = 0 and J U* J = U for rational shell weights (the "
            "Gaussian damping FORM is declared; the theorems hold for "
            "any symmetric positive weights). Asymmetric transport "
            "leaks EXACTLY ((w_R - w_L)/2)(e_R - e_L) per shell, "
            "vanishing iff the reflected weights agree. Exact transfer "
            "gives zero defect current and the CONVERSE FAILS: a common "
            "seam-sector error is invisible to the anti-seam projection "
            "— the anti-seam current detects only the antisymmetric "
            "discrepancy, so the seam component needs its own audit. "
            "Aligned/reversed/avoided pairings are +alpha, -alpha, 0; "
            "the damped obstruction keeps its sign with nonincreasing "
            "magnitude; the displacement prototype has the exact sign "
            "trichotomy. No prime-to-cycle equivalence is inferred"),
        "shells": list(SHELLS),
        "mismatch_shell": 5,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  zero-identification boundary (13)
# ----------------------------------------------------------------------

LEVELS = (Fr(1), Fr(3, 2), Fr(3, 2), Fr(2))      # with a double level


def paired_determinant(levels):
    p = [Fr(1)]
    for lam in levels:
        p = pmul(p, [Fr(1), Fr(0), Fr(-1) / (lam * lam)])
    return p


def certify_T5():
    D = paired_determinant(LEVELS)
    assert peval(D, Fr(0)) == 1                  # D(0) = 1
    assert all(c == 0 for i, c in enumerate(D) if i % 2 == 1)  # even

    # zeros exactly +-lambda_j WITH multiplicity
    for lam in set(LEVELS):
        mult = sum(1 for x in LEVELS if x == lam)
        assert root_multiplicity(D, lam) == mult
        assert root_multiplicity(D, -lam) == mult
    for probe in (Fr(1, 2), Fr(7, 4), Fr(5)):
        assert peval(D, probe) != 0

    # logarithmic-derivative identity, cross-multiplied (no division):
    # D' * prod_j (E^2 - lam_j^2) == D * sum_j 2E prod_{k!=j}(...)
    factors = [[-(lam * lam), Fr(0), Fr(1)] for lam in LEVELS]
    prod_all = [Fr(1)]
    for f in factors:
        prod_all = pmul(prod_all, f)
    rhs_sum = [Fr(0)]
    for j in range(len(LEVELS)):
        term = [Fr(0), Fr(2)]                    # 2E
        for k, f in enumerate(factors):
            if k != j:
                term = pmul(term, f)
        rhs_sum = padd(rhs_sum, term)
    lhs = pmul(pderiv(D), prod_all)
    rhs = pmul(D, rhs_sum)
    assert lhs == rhs                            # exact polynomial identity

    # nonvanishing-factor lemma on polynomial models:
    g = [Fr(1), Fr(0), Fr(1)]                    # E^2 + 1, no rational zeros
    F = pmul(g, D)
    for lam in set(LEVELS):
        assert root_multiplicity(F, lam) == root_multiplicity(D, lam)
        assert root_multiplicity(F, -lam) == root_multiplicity(D, -lam)
    for probe in (Fr(0), Fr(1, 2), Fr(3)):
        assert (peval(F, probe) == 0) == (peval(D, probe) == 0)

    # obstruction-class validator on synthetic multiset ledgers
    def audit(z_det, z_target):
        extra = {r: m for r, m in z_det.items()
                 if m > z_target.get(r, 0)}
        missing = {r: m for r, m in z_target.items()
                   if m > z_det.get(r, 0)}
        report = []
        if extra:
            report.append("Z1_extra")
        if missing:
            report.append("Z2_missing")
        if not extra and not missing and z_det != z_target:
            report.append("Z3_multiplicity")
        return report or ["CLEAN"]

    target = {str(lam): sum(1 for x in LEVELS if x == lam)
              for lam in set(LEVELS)}
    assert audit(dict(target), dict(target)) == ["CLEAN"]
    extra = dict(target)
    extra["5"] = 1
    assert audit(extra, target) == ["Z1_extra"]
    missing = dict(target)
    del missing["1"]
    assert audit(missing, target) == ["Z2_missing"]
    mism = dict(target)
    mism["3/2"] = 1
    assert "Z2_missing" in audit(mism, target)   # lowered multiplicity

    return {
        "statement": (
            "The finite paired determinant D_N(E) = prod(1 - "
            "E^2/lambda_j^2) satisfies D_N(0) = 1, is even at the "
            "coefficient level, and has zeros EXACTLY +-lambda_j with "
            "multiplicity counted by exact polynomial division (a double "
            "level gives a double zero); the logarithmic-derivative "
            "identity is certified in cross-multiplied polynomial form. "
            "Multiplying by a nonvanishing polynomial factor preserves "
            "the zero divisor with multiplicity (finite nonvanishing-"
            "factor lemma). The obstruction validator accepts a clean "
            "synthetic ledger and detects extra zeros, missing zeros and "
            "multiplicity mismatch — and passing it provides NO evidence "
            "that any synthetic set equals Z_xi. Zmodel -> Zxi and "
            "Zspectral -> Zxi remain FORBIDDEN promotions (13.7); "
            "finite symmetry is not completion (Remark 13.6)"),
        "levels": [str(x) for x in LEVELS],
        "determinant_coeffs": [str(c) for c in D],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T6  Surya/Weil fusion on rational models (14)
# ----------------------------------------------------------------------


def feshbach_congruence(A, C, Dm):
    """M = [[A, C^T],[C, D]];  S = [[I, A^{-1}C^T],[0, I]];
    returns (M, S, diag(A, F), F)."""
    p, q = len(A), len(Dm)
    Ct = mT(C)
    M = block2(A, Ct, C, Dm)
    Ainv = inv(A)
    F = msub(Dm, mmul(C, mmul(Ainv, Ct)))
    S = block2(meye(p), mmul(Ainv, Ct), mzero(q, p), meye(q))
    Dg = block2(A, mzero(p, q), mzero(q, p), F)
    return M, S, Dg, F


def ssq(x):
    """Squared rectified Surya amplitude: sin(arctan(-x))_+^2, exactly."""
    return x * x / (1 + x * x) if x < 0 else Fr(0)


def csq(mu):
    """Squared Surya angle amplitude at relative eigenvalue mu, with the
    sign carried separately: returns (sign, squared amplitude)."""
    d = mu - 1
    s = 0 if d == 0 else (1 if d > 0 else -1)
    return s, d * d / (1 + d * d)


def certify_T6():
    # ---- Feshbach congruence as an exact identity, PSD both ways ----
    A = [[Fr(2), Fr(1)], [Fr(1), Fr(2)]]         # A > 0
    assert all(m > 0 for m in principal_minors(A))
    C = [[Fr(1), Fr(-1)], [Fr(0), Fr(2)]]

    D_pos = [[Fr(4), Fr(0)], [Fr(0), Fr(6)]]
    M, S, Dg, F = feshbach_congruence(A, C, D_pos)
    assert M == mmul(mT(S), mmul(Dg, S))         # THE congruence, exactly
    assert is_psd(F) and is_psd(M)               # F >= 0  =>  M >= 0

    D_neg = [[Fr(1), Fr(0)], [Fr(0), Fr(1)]]
    M2, S2, Dg2, F2 = feshbach_congruence(A, C, D_neg)
    assert M2 == mmul(mT(S2), mmul(Dg2, S2))
    assert not is_psd(F2)
    z = None
    for cand in ([Fr(1), Fr(0)], [Fr(0), Fr(1)], [Fr(1), Fr(1)],
                 [Fr(1), Fr(-1)], [Fr(2), Fr(1)]):
        if dot(cand, mv(F2, cand)) < 0:
            z = cand
            break
    assert z is not None
    x = mv(inv(S2), [Fr(0), Fr(0)] + z)          # transported witness
    val = dot(x, mv(M2, x))
    assert val == dot(z, mv(F2, z)) < 0          # M2 not PSD, exactly
    assert not is_psd(M2)

    # ---- rectified Surya without square roots ----
    for xr in (Fr(-3), Fr(-1, 2), Fr(0), Fr(1, 4), Fr(5)):
        assert (ssq(xr) == 0) == (xr >= 0)
    # monotone amplitude on the forbidden side
    assert ssq(Fr(-3)) > ssq(Fr(-1)) > ssq(Fr(-1, 2)) > 0

    # operator version on a rational-eigenvalue model: W0 = Q L Q^T,
    # Q the 3-4-5 rotation (orthogonal, rational)
    Q = [[Fr(3, 5), Fr(-4, 5)], [Fr(4, 5), Fr(3, 5)]]
    assert mmul(mT(Q), Q) == meye(2)
    for spec, nonneg in (((Fr(2), Fr(1, 3)), True),
                         ((Fr(2), Fr(-1, 4)), False),
                         ((Fr(0), Fr(5)), True)):
        L = [[spec[0], Fr(0)], [Fr(0), spec[1]]]
        W0 = mmul(Q, mmul(L, mT(Q)))
        Ass_sq = mmul(Q, mmul(
            [[ssq(spec[0]), Fr(0)], [Fr(0), ssq(spec[1])]], mT(Q)))
        vanish = Ass_sq == mzero(2)
        assert vanish == nonneg == is_psd(W0)
        # A_SS itself is PSD with A_SS = 0 <=> A_SS^2 = 0, so the
        # squared amplitude decides the terminal sign exactly

    # ---- seam charge, additivity, spectral flow (14.4-14.5) ----
    spec_even = (Fr(2), Fr(1, 2))                # k+ = 1
    spec_odd = (Fr(3), Fr(-1), Fr(1))            # k3 = 1 (mu = 1 NOT > 1)
    k_plus = sum(1 for m in spec_even if m > 1)
    k_odd = sum(1 for m in spec_odd if m > 1)
    spec_sum = spec_even + spec_odd
    k_sigma = sum(1 for m in spec_sum if m > 1)
    assert k_sigma == k_plus + k_odd == 2        # additive, no cancellation

    # spectral flow of 1 - t mu: crossing iff mu > 1, always downward
    sf = 0
    for mu in spec_sum:
        if mu != 0:
            t_cross = Fr(1) / mu
            crosses = Fr(0) < t_cross < Fr(1)
            assert crosses == (mu > 1)
            if crosses:
                slope = -mu
                assert slope < 0                 # one-sided orientation
                sf -= 1
    assert sf == -k_sigma

    # W_eta = A^{1/2}(I - K)A^{1/2} with rational square roots:
    # inertia transported by congruence, k_sigma = 0 <=> W >= 0
    Ahalf = [[Fr(2), Fr(0)], [Fr(0), Fr(3)]]     # A = diag(4, 9) > 0
    for spec, nonneg in (((Fr(1, 2), Fr(-1)), True),
                         ((Fr(2), Fr(1, 2)), False)):
        K = [[spec[0], Fr(0)], [Fr(0), spec[1]]]
        ImK = msub(meye(2), K)
        W = mmul(Ahalf, mmul(ImK, Ahalf))
        k_here = sum(1 for m in spec if m > 1)
        assert is_psd(W) == is_psd(ImK) == (k_here == 0) == nonneg

    # ---- Surya angle indicator equality (14.6) ----
    for mu in (Fr(3), Fr(1), Fr(1, 2), Fr(-2), Fr(7, 5)):
        s, amp2 = csq(mu)
        positive_part_nonzero = (s > 0 and amp2 > 0)
        assert positive_part_nonzero == (mu > 1)
        assert (s == 0) == (mu == 1)

    # ---- principal holonomy is insufficient (14.8) ----
    for m in (-3, -1, 0, 1, 2, 5):
        assert m % 1 == 0                        # principal phase trivial
    assert any(m != 0 for m in (-3, -1, 0, 1, 2, 5))
    # the lift carries strictly more than its exponential image — the
    # same projection blindness certified in UGD-1 T4

    return {
        "statement": (
            "The Feshbach congruence M = S^T diag(A, F) S is an EXACT "
            "matrix identity, so with A > 0: W_- >= 0 <=> F >= 0, both "
            "directions certified (principal minors + a negative "
            "witness transported through the congruence with exact "
            "equality of quadratic forms). The rectified Surya "
            "correction is handled without square roots via the squared "
            "amplitude ssq(x) = x^2/(1+x^2) on x < 0: ssq = 0 <=> "
            "x >= 0, and on rational-eigenvalue models A_SS = 0 <=> "
            "A_SS^2 = 0 <=> W_0 >= 0 — the two terminal Weil signs are "
            "ONE operator obstruction. The seam charge k_Sigma = "
            "#{mu > 1} is additive across the even/odd sum (mu = 1 "
            "counts strictly, NOT in the charge), equals minus the "
            "spectral flow of I - tK with every crossing t = 1/mu "
            "downward (slope -mu < 0: cancellation is impossible), and "
            "W_eta = A^(1/2)(I-K)A^(1/2) with rational square roots "
            "gives k_Sigma = 0 <=> W_eta >= 0. Surya angle indicator: "
            "positive part of c_SS nonzero <=> mu > 1. Principal "
            "holonomy m mod 1 = 0 for every integer m while the lift m "
            "need not vanish: the no-leakage theorem must kill the "
            "INTEGER, not its exponential image. The completed W_3, the "
            "Weil-to-RH equivalence and the odd compact reference are "
            "IMPORTED/OPEN in the source and are NOT certified here"),
        "k_plus": k_plus,
        "k_odd": k_odd,
        "k_sigma": k_sigma,
        "spectral_flow": sf,
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "RST-2: recognition-seam topology, arithmetic and "
                   "fusion spine",
        "source": {
            "paper": ("Recognition-Seam Topology and Daggered Index "
                      "Geometry (M. Dabas, July 2026, arXiv source v1)"),
            "sections_certified": "8-14 (finite exact content)",
            "companion": "RST-1 (sections 2-7), same folder",
        },
        "T1_sampled_reflection_and_blindness": certify_T1(),
        "T2_weights_leakage_norm_crossings": certify_T2(),
        "T3_arithmetic_recognition_formal_logs": certify_T3(),
        "T4_prime_shell_transfer_and_currents": certify_T4(),
        "T5_zero_identification_boundary": certify_T5(),
        "T6_surya_weil_fusion": certify_T6(),
        "claim_boundary": {
            "analytic_layer": (
                "DECLARED, NOT CERTIFIED: Thm 8.8 (entire trace-class "
                "damped family), Thm 9.1 (undamped strip "
                "classification), the infinite part of Thm 8.5 "
                "(Schatten thresholds), Thm 9.4 (archimedean defect), "
                "and every completed-carrier statement"),
            "declared_forms": (
                "the Gaussian damping e^{-lambda(log p)^2} and the "
                "1/sqrt(2) transfer normalization are DECLARED forms; "
                "every certified identity holds for arbitrary symmetric "
                "positive rational weights and is normalization-"
                "independent for vanishing statements"),
            "imported_not_certified": (
                "the completed Weil operator W_3, the Weil-to-RH "
                "equivalence, the parity/Feshbach reduction ON THE "
                "COMPLETED CARRIER, and the even compact shifted "
                "relative operator are imported by the source from the "
                "MP X3 stack; RST-2 certifies their finite rational "
                "SHAPE only"),
            "open_gates": (
                "F1 odd compact reference, F2 even charge, F3 odd "
                "charge, F4 completion stability as eta -> 0, F5 the "
                "no-leakage theorem: ALL OPEN, exactly as the source "
                "states"),
            "forbidden_promotions": (
                "Zmodel -> Zxi and Zspectral -> Zxi are FORBIDDEN; the "
                "finite sign is a witness, not a Fredholm index; the "
                "clean validator run is no evidence about actual zeros"),
            "riemann_hypothesis": "ABSTAIN",
        },
        "provenance": {
            "prior_executable_version": (
                "NONE found for the paper's arithmetic/fusion spine; "
                "first certification"),
            "companions": ("RST-1; EMK-1, EMK-2, UGD-1, EMK-T1 "
                           "(papers/emk-ugd-algebra); EMK-G1 "
                           "(papers/emk-recognition-geometry)"),
            "projection_blindness_thread": (
                "UGD-1 T4 (classical projection blind to seam), EMK-T1 "
                "T5 (clock equality is not time closure), RST-2 T1 "
                "(seam projection blind to anti-seam part) and RST-2 "
                "principal-holonomy insufficiency are one shape: the "
                "visible coordinate never carries the full state"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "RST2_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    print("RST-2 certificate written:", out)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
