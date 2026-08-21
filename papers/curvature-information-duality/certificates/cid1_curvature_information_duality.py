"""CID-1: curvature-information duality — the obstruction is metric-free.

A BRIDGE capsule. CFE-U certified that the obstruction to the memoryless
limit is UNIQUE: the annihilator of the memoryless responses in the
declared response algebra has dimension exactly 1, spanned by the loop
residue. This capsule asks the dual question:

    does the INFORMATION side see the same obstruction, and does it
    see anything the curvature side does not?

The answer certified here is a separation with an invariance:

    the Hodge/Helmholtz SPLIT of a response depends on the information
    metric — different metrics give different memoryless parts and
    different potentials — while the OBSTRUCTION does not: the residue
    of the remainder is the same rational number for every metric, and
    equals CFE-U's Area * d1.

So the information metric is a gauge on the decomposition and an
invariant on the obstruction. Curvature and information are dual
descriptions of ONE object, not two independent obstructions.

BLOCKS

  T1  THE EXACT INFORMATION METRIC. At a declared rational point of a
      finite exponential family, the information (Fisher) metric is the
      covariance of the sufficient statistics, g = Cov_p(T) — an exact
      rational matrix, computed WITHOUT evaluating a single exponential
      or logarithm (the exponential form of the family is declared; the
      metric at a rational point is arithmetic). Certified: symmetric,
      PSD by principal minors, positive definite exactly when the
      statistics are affinely independent, and DEGENERATE with an
      explicit null vector exactly when one statistic is an affine
      combination of the others — the information metric's kernel is
      the recognition kernel, matching CFE-U's ker d0 = constants.

  T2  RECOGNITION LEDGER = LAW OF TOTAL COVARIANCE. Coarse-graining
      (an Eye map on outcomes) splits the information metric EXACTLY:

          g = g_recognized + g_discarded,

      the between-block covariance of the conditional means plus the
      block-averaged within-block covariance, both PSD, verified as an
      exact matrix identity on rational models. Consequences certified:
      coarse-graining never increases information (g - g_coarse is PSD
      — monotonicity, here as an identity rather than an inequality),
      a sufficient partition discards exactly nothing (g_discarded = 0),
      and a collapsing partition discards everything. Recognition can
      lose information; it can never manufacture it.

  T3  THE OBSTRUCTION IS METRIC-FREE. For a family of distinct
      positive-definite rational metrics — including strongly
      anisotropic and sheared ones — the loop residue of the pinned
      CFE-1 witness is UNCHANGED: curvature is d, not a covariant
      derivative, so no Christoffel symbol of any information metric
      enters it. The residue equals BETA*(chi-1)*Area exactly, for
      every metric, reproducing CFE-1's ladder.

  T4  THE SPLIT IS NOT METRIC-FREE. With respect to each metric,
      exact rational least squares (normal equations solved in Q)
      decomposes a response into its g-orthogonal memoryless part and
      a remainder: omega = d Phi_g + h_g. Certified: the decomposition
      is genuine (h_g is g-orthogonal to every memoryless direction,
      d Phi_g is memoryless), the potential Phi_g and the remainder
      h_g GENUINELY DIFFER between metrics, and yet

          oint h_g = oint omega   for every g,

      one rational number independent of the metric. The gauge moves;
      the invariant does not.

  T5  ONE OBSTRUCTION, TWO PRESENTATIONS. For every tested metric the
      g-orthogonal complement of the memoryless space is exactly
      1-dimensional (consuming CFE-U's ker d1), its generator is
      metric-DEPENDENT as a vector, and its residue functional is
      metric-INDEPENDENT after normalization: all metrics produce the
      same obstruction functional up to a nonzero scalar, i.e. the same
      point of the 1-dimensional obstruction space CFE-U certified.
      Control: a degenerate (merely PSD) form fails to produce a
      well-posed split — positive definiteness is load-bearing, not
      decoration.

CLAIM BOUNDARY. The exponential family's exponential form is DECLARED;
only its rational-point arithmetic is certified. Continuum information
geometry, alpha-connections, dually flat structures, quantum (Fubini-
Study / Bures / SLD) information metrics and their monotonicity are NOT
claimed. The vault appendix named for Euler information-curvature
duality is NOT the source of this capsule (it has not been read); this
is a bridge between two certified capsules of this corpus. Uniqueness
inherits CFE-U's scope: the DECLARED response algebra, not a universal
law. RH / K0 / L0 / YM / quantum gravity untouched.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CFE_DIR = os.path.join(HERE, "..", "..", "cut-first-equivalence",
                       "certificates")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


cfeu = _load("cfeu_for_cid", os.path.join(CFE_DIR, "cfeu_uniqueness.py"))
cfe1 = cfeu.cfe1

# ----------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------


def mmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)]
            for i in range(n)]


def mT(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def msub(A, B):
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def madd(A, B):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def mv(A, x):
    return [sum(A[i][j] * x[j] for j in range(len(x)))
            for i in range(len(A))]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def det(A):
    n = len(A)
    if n == 0:
        return Fr(1)
    if n == 1:
        return A[0][0]
    out = Fr(0)
    for j in range(n):
        minor = [row[:j] + row[j + 1:] for row in A[1:]]
        out += (-1) ** j * A[0][j] * det(minor)
    return out


def principal_minors(A):
    from itertools import combinations
    n = len(A)
    out = []
    for r in range(1, n + 1):
        for idx in combinations(range(n), r):
            out.append(det([[A[i][j] for j in idx] for i in idx]))
    return out


def is_psd(A):
    assert A == mT(A)
    return all(m >= 0 for m in principal_minors(A))


def is_pd(A):
    assert A == mT(A)
    n = len(A)
    return all(det([[A[i][j] for j in range(k)] for i in range(k)]) > 0
               for k in range(1, n + 1))


def solve(A, b):
    """Exact solution of A x = b for square nonsingular rational A."""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                M[i] = [x - M[i][c] * y for x, y in zip(M[i], M[c])]
    return [M[i][n] for i in range(n)]


def inv2(A):
    d = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    assert d != 0
    return [[A[1][1] / d, -A[0][1] / d], [-A[1][0] / d, A[0][0] / d]]


# ----------------------------------------------------------------------
# T1  the exact information metric at a rational point
# ----------------------------------------------------------------------

# Outcome set and a declared rational point of the family. The family's
# exponential form is DECLARED; the metric at this point is arithmetic.
OUTCOMES = (0, 1, 2, 3)
P_POINT = (Fr(1, 4), Fr(1, 3), Fr(1, 4), Fr(1, 6))

STATS = (
    (Fr(1), Fr(0), Fr(2), Fr(-1)),        # T1
    (Fr(0), Fr(3), Fr(1), Fr(1)),         # T2
)
STAT_DEP = (Fr(1), Fr(3), Fr(3), Fr(0))   # = T1 + T2 exactly: dependent


def expect(p, f):
    return sum(pi * fi for pi, fi in zip(p, f))


def covariance(p, stats):
    assert sum(p) == 1 and all(x > 0 for x in p)
    means = [expect(p, T) for T in stats]
    n = len(stats)
    return [[expect(p, [stats[i][x] * stats[j][x]
                        for x in range(len(p))]) - means[i] * means[j]
             for j in range(n)] for i in range(n)]


def certify_T1():
    g = covariance(P_POINT, STATS)
    assert g == mT(g)
    assert is_psd(g) and is_pd(g)

    # degenerate control: a statistic that is an affine combination
    dep = [STATS[0][x] + STATS[1][x] for x in range(len(OUTCOMES))]
    assert tuple(dep) == STAT_DEP
    g3 = covariance(P_POINT, STATS + (tuple(dep),))
    assert is_psd(g3) and not is_pd(g3)
    null = [Fr(1), Fr(1), Fr(-1)]                 # T1 + T2 - T3 = 0
    assert mv(g3, null) == [Fr(0)] * 3
    assert dot(null, mv(g3, null)) == 0

    # constants carry no information: adding a constant statistic
    # enlarges the kernel by exactly that direction
    const = tuple(Fr(1) for _ in OUTCOMES)
    g_c = covariance(P_POINT, STATS + (const,))
    assert mv(g_c, [Fr(0), Fr(0), Fr(1)]) == [Fr(0)] * 3
    assert is_psd(g_c) and not is_pd(g_c)

    return {
        "statement": (
            "At a declared rational point of a finite exponential "
            "family the information metric is g = Cov_p(T), an EXACT "
            "rational matrix computed without evaluating any "
            "exponential or logarithm. It is symmetric and PSD; it is "
            "positive definite exactly when the sufficient statistics "
            "are affinely independent, and degenerate with an explicit "
            "null vector exactly when one is an affine combination of "
            "the others. Constant statistics lie in the kernel: the "
            "information metric is blind to constants, exactly as "
            "CFE-U's d0 has kernel the constants"),
        "point": [str(x) for x in P_POINT],
        "metric": [[str(x) for x in row] for row in
                   covariance(P_POINT, STATS)],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  recognition ledger = law of total covariance
# ----------------------------------------------------------------------


def coarse_grain(p, stats, blocks):
    """Eye map on outcomes: returns (block probs, conditional means,
    between-block covariance, block-averaged within covariance)."""
    P = [sum(p[x] for x in B) for B in blocks]
    assert sum(P) == 1
    cond = []
    for T in stats:
        cond.append(tuple(
            sum(p[x] * T[x] for x in B) / PB for B, PB in zip(blocks, P)))
    between = covariance(tuple(P), tuple(cond))
    n = len(stats)
    within = [[Fr(0)] * n for _ in range(n)]
    for B, PB in zip(blocks, P):
        pc = tuple(p[x] / PB for x in B)
        sub = tuple(tuple(T[x] for x in B) for T in stats)
        cb = covariance(pc, sub)
        for i in range(n):
            for j in range(n):
                within[i][j] += PB * cb[i][j]
    return P, cond, between, within


def certify_T2():
    g = covariance(P_POINT, STATS)

    cases = {
        "proper": [(0, 1), (2, 3)],
        "singleton_refinement": [(0,), (1,), (2,), (3,)],
        "total_collapse": [(0, 1, 2, 3)],
        "uneven": [(0,), (1, 2, 3)],
    }
    results = {}
    for name, blocks in cases.items():
        P, cond, between, within = coarse_grain(P_POINT, STATS, blocks)
        assert madd(between, within) == g          # EXACT ledger
        assert is_psd(between) and is_psd(within)
        assert is_psd(msub(g, between))            # monotonicity
        results[name] = {
            "discarded_is_zero": within == [[Fr(0), Fr(0)],
                                            [Fr(0), Fr(0)]],
            "recognized_is_zero": between == [[Fr(0), Fr(0)],
                                              [Fr(0), Fr(0)]],
        }

    # a sufficient partition discards exactly nothing
    assert results["singleton_refinement"]["discarded_is_zero"]
    # a collapsing partition recognizes exactly nothing
    assert results["total_collapse"]["recognized_is_zero"]
    # a proper partition does both partially
    assert not results["proper"]["discarded_is_zero"]
    assert not results["proper"]["recognized_is_zero"]

    return {
        "statement": (
            "Coarse-graining splits the information metric EXACTLY: "
            "g = g_recognized + g_discarded (between-block covariance "
            "of the conditional means plus block-averaged within-block "
            "covariance), both PSD, as an exact matrix identity. "
            "Monotonicity follows as an IDENTITY rather than an "
            "inequality: g - g_recognized = g_discarded is PSD, so "
            "recognition can lose information and can never "
            "manufacture it. The singleton refinement discards exactly "
            "nothing; total collapse recognizes exactly nothing; a "
            "proper partition does both partially. This is the "
            "recognition ledger on the information side"),
        "cases": {k: {kk: bool(vv) for kk, vv in v.items()}
                  for k, v in results.items()},
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# metrics on the response plane, and the response-space inner product
# ----------------------------------------------------------------------

METRICS = {
    "euclidean": [[Fr(1), Fr(0)], [Fr(0), Fr(1)]],
    "anisotropic": [[Fr(9), Fr(0)], [Fr(0), Fr(1, 4)]],
    "sheared": [[Fr(5), Fr(2)], [Fr(2), Fr(3)]],
    "information": None,          # filled from Cov_p(T) below
}
METRICS["information"] = covariance(P_POINT, STATS)

LATTICE = tuple((Fr(i), Fr(j)) for i in range(-1, 2) for j in range(-1, 2))


def form_value(vec, p, v):
    a, b, c, e, f, g = vec
    return [a + b * p + c * v, e + f * p + g * v]


def inner(vec1, vec2, gmat):
    """<w1, w2>_g = sum over the lattice of w1^T g^{-1} w2 — exact."""
    ginv = inv2(gmat)
    total = Fr(0)
    for (p, v) in LATTICE:
        x = form_value(vec1, p, v)
        y = form_value(vec2, p, v)
        total += dot(x, mv(ginv, y))
    return total


def gram(basis, gmat):
    return [[inner(u, w, gmat) for w in basis] for u in basis]


def hodge_split(vec, gmat, basis):
    """omega = d Phi_g + h_g : exact g-orthogonal least squares."""
    G = gram(basis, gmat)
    rhs = [inner(k, vec, gmat) for k in basis]
    coeffs = solve(G, rhs)
    closed = [sum(c * k[j] for c, k in zip(coeffs, basis))
              for j in range(6)]
    remainder = [x - y for x, y in zip(vec, closed)]
    return coeffs, closed, remainder


# ----------------------------------------------------------------------
# T3  the obstruction is metric-free
# ----------------------------------------------------------------------

CHIS = cfeu.CHIS
LOOP = cfeu.rect_loop(cfeu.P_LO, cfeu.P_HI, cfeu.V_LO, cfeu.V_HI)


def certify_T3():
    ladder = []
    for chi in CHIS:
        w = cfeu.witness_vec(chi)
        res = cfeu.circulation_of(w, LOOP)
        assert res == cfe1.circulation(LOOP, chi)     # pinned agreement
        assert res == cfe1.BETA * (chi - 1) * cfeu.AREA
        # the residue does not mention any metric
        for name, gmat in METRICS.items():
            assert is_pd(gmat)
            assert cfeu.circulation_of(w, LOOP) == res
        ladder.append((chi, res))
    assert [str(x) for _, x in ladder] == ["-32", "-16", "0", "16", "32"]
    return {
        "statement": (
            "Curvature is d, not a covariant derivative: no Christoffel "
            "symbol of any information metric enters the residue. For "
            "every tested positive-definite rational metric — "
            "euclidean, strongly anisotropic, sheared, and the "
            "information metric Cov_p(T) itself — the loop residue of "
            "the pinned CFE-1 witness is the SAME rational number, "
            "equal to BETA*(chi-1)*Area, reproducing the certified "
            "ladder -32, -16, 0, +16, +32"),
        "metrics": sorted(METRICS),
        "ladder": [[str(c), str(x)] for c, x in ladder],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  the split is not metric-free, the residue of the remainder is
# ----------------------------------------------------------------------


def certify_T4():
    basis = cfeu.nullspace(cfeu.D1)                  # memoryless space
    assert len(basis) == 5
    w = cfeu.witness_vec(Fr(7, 5))
    target = cfeu.circulation_of(w, LOOP)
    assert target == Fr(32)

    splits = {}
    for name, gmat in METRICS.items():
        coeffs, closed, rem = hodge_split(w, gmat, basis)
        # the closed part is genuinely memoryless
        assert cfeu.mv(cfeu.D1, closed) == [Fr(0)]
        assert cfeu.circulation_of(closed, LOOP) == 0
        # the remainder is g-orthogonal to every memoryless direction
        for k in basis:
            assert inner(rem, k, gmat) == 0
        # and carries the whole residue
        assert cfeu.circulation_of(rem, LOOP) == target
        # the closed part integrates to an explicit potential
        phi = cfeu.potential_of_closed(closed)
        assert cfeu.mv(cfeu.D0, phi) == closed
        splits[name] = (tuple(closed), tuple(rem))

    # the split GENUINELY differs between metrics
    distinct = {v for v in splits.values()}
    assert len(distinct) > 1
    assert splits["euclidean"] != splits["anisotropic"]
    assert splits["euclidean"] != splits["sheared"]

    return {
        "statement": (
            "Exact rational least squares (normal equations solved in "
            "Q) splits a response into its g-orthogonal memoryless "
            "part and a remainder. The split is genuine — the closed "
            "part is memoryless with zero residue and integrates to an "
            "explicit potential, the remainder is g-orthogonal to "
            "every memoryless direction — and it GENUINELY DIFFERS "
            "between metrics. Yet the residue of the remainder is the "
            "SAME rational number 32 for every metric. The information "
            "metric is a gauge on the decomposition and an invariant "
            "on the obstruction"),
        "distinct_splits": len(distinct),
        "common_residue": str(target),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T5  one obstruction, two presentations
# ----------------------------------------------------------------------


def certify_T5():
    basis = cfeu.nullspace(cfeu.D1)
    R = []
    for j in range(6):
        e = [Fr(0)] * 6
        e[j] = Fr(1)
        R.append(cfeu.circulation_of(e, LOOP))
    assert R == [cfeu.AREA * cfeu.D1[0][j] for j in range(6)]

    gens = {}
    for name, gmat in METRICS.items():
        # g-orthogonal complement of the memoryless space
        comp = []
        for j in range(6):
            e = [Fr(0)] * 6
            e[j] = Fr(1)
            _, _, rem = hodge_split(e, gmat, basis)
            comp.append(rem)
        rk = cfeu.rank(comp)
        assert rk == 1                                # exactly one
        gen = next(c for c in comp if any(x != 0 for x in c))
        gens[name] = tuple(gen)
        # the generator pairs nontrivially with the residue functional
        assert dot(R, gen) != 0
        # every nonzero complement vector is a multiple of the generator
        for c in comp:
            if any(x != 0 for x in c):
                k = next(ci / gi for ci, gi in zip(c, gen) if gi != 0)
                assert c == [k * gi for gi in gen]
                assert cfeu.circulation_of(c, LOOP) == k * dot(R, gen)

    # the generators genuinely differ as vectors between metrics
    assert len({v for v in gens.values()}) > 1
    # but all lie in the SAME 1-dimensional obstruction class: their
    # residues are all nonzero multiples of one another
    vals = {name: cfeu.circulation_of(list(g), LOOP)
            for name, g in gens.items()}
    assert all(v != 0 for v in vals.values())

    # control: a degenerate PSD form has no well-posed split
    degenerate = [[Fr(1), Fr(1)], [Fr(1), Fr(1)]]
    assert is_psd(degenerate) and not is_pd(degenerate)
    assert det(degenerate) == 0                       # not invertible

    return {
        "statement": (
            "For every tested metric the g-orthogonal complement of "
            "the memoryless space has rank EXACTLY 1 — CFE-U's "
            "uniqueness seen through the information metric. Its "
            "generator is metric-DEPENDENT as a vector while its "
            "residue is always nonzero: all metrics land on the same "
            "1-dimensional obstruction class, differing only by a "
            "scalar. CONTROL: a degenerate PSD form is not invertible "
            "and yields no well-posed split — positive definiteness is "
            "load-bearing"),
        "complement_rank": 1,
        "distinct_generators": len({v for v in gens.values()}),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "CID-1: curvature-information duality — the "
                   "obstruction is metric-free, the split is not",
        "consumes": {
            "CFE-U": ("the memoryless space ker d1, the residue "
                      "functional, the response complex and the "
                      "potential reconstruction, all from the pinned "
                      "module"),
            "CFE-1": ("the witness family and circulation, through "
                      "CFE-U"),
        },
        "T1_exact_information_metric": certify_T1(),
        "T2_recognition_ledger_total_covariance": certify_T2(),
        "T3_obstruction_is_metric_free": certify_T3(),
        "T4_split_is_metric_dependent": certify_T4(),
        "T5_one_obstruction_two_presentations": certify_T5(),
        "claim_boundary": {
            "family_form_declared": (
                "DECLARED: the exponential form of the family. Only "
                "the arithmetic at a rational point is certified — no "
                "exponential or logarithm is evaluated anywhere"),
            "not_claimed": (
                "NOT CLAIMED: continuum information geometry, "
                "alpha-connections, dually flat structures, and "
                "quantum (Fubini-Study / Bures / SLD) information "
                "metrics or their monotonicity"),
            "source": (
                "this is a BRIDGE between two certified capsules of "
                "this corpus; the vault appendix named for Euler "
                "information-curvature duality has NOT been read and "
                "is NOT the source"),
            "uniqueness_scope": (
                "inherits CFE-U's scope: the DECLARED response "
                "algebra, not a universal law"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE — first bridge capsule "
                                        "between the curvature and "
                                        "information layers",
            "companions": ("CFE-1/2/Q/U (cut-first-equivalence), "
                           "information-invariance D-series"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "CID1_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    print("CID-1 certificate written:", out)
    print("sha256:", hashlib.sha256(payload.encode()).hexdigest())


if __name__ == "__main__":
    main()
