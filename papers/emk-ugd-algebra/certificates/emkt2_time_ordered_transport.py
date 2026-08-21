"""EMK-T2: time-ordered transport — order is content.

Dynamic companion to EMK-T1. That capsule certified the master-tensor
appendix's STATIC spine: the channel curvature decomposition, the
necessity of the mixed term (T3), and the separation of time closure
from clock equality (T5). This capsule certifies the same statements
IN TRANSPORT: when the generators actually evolve a state, the order
of evolution is content, the mixed term obstructs order-independence,
and equal clocks with different histories are separated by an exact
holonomy that a lawful ledger closes to zero.

The machinery is the exact-exponential move of RST-1 T4: for nilpotent
generators every exponential is a FINITE polynomial, so every transport,
loop and discrepancy below is an identity of rational matrices — no
limit, no series truncation, no float.

BLOCKS

  T1  ORDER IS CONTENT. For K1 = E12, K2 = E23 and rational times a, b:

          e^{aK1} e^{bK2}  =  I + aK1 + bK2 + ab E13
          e^{bK2} e^{aK1}  =  I + aK1 + bK2
          e^{aK1 + bK2}    =  I + aK1 + bK2 + (ab/2) E13

      three EXACTLY DIFFERENT transports, and every pairwise
      discrepancy is an exact multiple of the commutator
      [K1, K2] = E13: the two orders differ by ab[K1,K2], and the
      unordered sum-exponential sits exactly halfway between them.
      Commuting control (K2' = E13): all three coincide identically.

  T2  REFINEMENT IS FREE, REVERSAL IS NOT. Within one constant
      generator, e^{K(a+b)} = e^{Ka} e^{Kb} exactly, so subdividing a
      time interval changes nothing: certified over several rational
      partitions of the same total time, all products equal. Reversing
      the order of two DISTINCT-generator segments changes the
      transport by exactly the T1 commutator term. Time-ordering, not
      time-slicing, is the load-bearing structure.

  T3  THE MIXED TERM OBSTRUCTS TRANSPORT (dynamic EMK-T1 T3). The
      pinned EMK-T1 channels — each individually flat — have
      channel-internal transport loops that close EXACTLY to I. The
      CROSS-channel loop

          L(a, b) = e^{aE12} e^{bE21} e^{-aE12} e^{-bE21}
                  = I + ab [E12, E21] + (exact remainder)

      is written in closed form as a rational matrix identity: the
      leading bilinear part is EXACTLY ab times the cross-commutator,
      L != I for ab != 0, and det L = 1. Channel-by-channel transport
      audit is provably insufficient IN DYNAMICS, not only in algebra.

  T4  SAME CLOCK, DIFFERENT HISTORY, EXACT LEDGER (dynamic EMK-T1
      T5). Two evolutions through the same generators for the same
      total time (identical clock tau = a + b) but in opposite order
      end at transports differing by exactly ab[K1, K2]: the clock
      agrees, the state does not. The temporal residue (consumed from
      the pinned EMK-T1 residue functional) is nonzero exactly then,
      and the lawful ledger — the exact inverse of the discrepancy —
      closes it to 0 identically. Fourth dynamic appearance of the
      projection-blindness shape (UGD-1 T4, EMK-T1 T5, RST-2 T1).

CLAIM BOUNDARY. Continuous-time ordered exponentials (P exp), Magnus
expansions, non-nilpotent generators and refinement limits are NOT
claimed — the source (topology paper 5.4) states the ordered
exponential as the lawful general object and this capsule certifies
only its exact nilpotent sector. The left-action convention (later
transport multiplies on the left) is declared. RH / K0 / L0 / YM /
quantum gravity untouched.
"""

import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(HERE, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


t1 = _load("emkt1_for_t2", "emkt1_master_tensor_and_time.py")

# ----------------------------------------------------------------------
# exact matrices and nilpotent exponentials
# ----------------------------------------------------------------------


def eye(n):
    return [[Fr(1) if i == j else Fr(0) for j in range(n)] for i in range(n)]


def zeros(n):
    return [[Fr(0)] * n for _ in range(n)]


def mm(A, B):
    n = len(A)
    return [[sum(A[i][k] * B[k][j] for k in range(n)) for j in range(n)]
            for i in range(n)]


def madd(A, B):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def msub(A, B):
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(c, A):
    return [[c * x for x in row] for row in A]


def comm(A, B):
    return msub(mm(A, B), mm(B, A))


def E(n, i, j):
    M = zeros(n)
    M[i][j] = Fr(1)
    return M


def nilpotent_exp(t, K):
    """e^{t K} for nilpotent K, as an EXACT finite polynomial sum."""
    n = len(K)
    out = eye(n)
    term = eye(n)
    fact = 1
    for k in range(1, n + 1):
        term = mm(term, K)
        if term == zeros(n):
            break
        fact *= k
        out = madd(out, mscale(Fr(t) ** k / fact, term))
    return out


def ordered_product(segments):
    """Transport of time-ordered segments [(t1,K1),(t2,K2),...]:
    later segments multiply on the LEFT (declared convention)."""
    n = len(segments[0][1])
    W = eye(n)
    for t, K in segments:
        W = mm(nilpotent_exp(t, K), W)
    return W


# ----------------------------------------------------------------------
# T1  order is content
# ----------------------------------------------------------------------

A_T, B_T = Fr(2, 3), Fr(5, 7)


def certify_T1():
    n = 3
    K1, K2 = E(n, 0, 1), E(n, 1, 2)
    C = comm(K1, K2)
    assert C == E(n, 0, 2)

    a, b = A_T, B_T
    W_12 = ordered_product([(a, K1), (b, K2)])       # K1 first, then K2
    W_21 = ordered_product([(b, K2), (a, K1)])       # K2 first, then K1
    W_sum = nilpotent_exp(1, madd(mscale(a, K1), mscale(b, K2)))

    base = madd(eye(n), madd(mscale(a, K1), mscale(b, K2)))
    # left-action: W_12 = e^{bK2} e^{aK1} = base (E23 E12 = 0), while
    # W_21 = e^{aK1} e^{bK2} = base + ab E13 (E12 E23 = E13)
    assert W_12 == base
    assert W_21 == madd(base, mscale(a * b, C))
    assert W_sum == madd(base, mscale(a * b / 2, C))

    # every pairwise discrepancy is an exact commutator multiple
    assert msub(W_21, W_12) == mscale(a * b, C)
    assert msub(W_21, W_sum) == mscale(a * b / 2, C)
    assert msub(W_sum, W_12) == mscale(a * b / 2, C)
    # and the sum-exponential is exactly halfway between the orders
    assert mscale(Fr(2), W_sum) == madd(W_12, W_21)

    # commuting control: all three coincide identically
    K2c = E(n, 0, 2)
    assert comm(K1, K2c) == zeros(n)
    U = ordered_product([(a, K1), (b, K2c)])
    V = ordered_product([(b, K2c), (a, K1)])
    S = nilpotent_exp(1, madd(mscale(a, K1), mscale(b, K2c)))
    assert U == V == S

    return {
        "statement": (
            "Three transports of the same generator content — the two "
            "time-orders and the unordered sum-exponential — are "
            "EXACTLY different rational matrices, every pairwise "
            "discrepancy an exact multiple of [K1,K2] (ab, ab/2, "
            "ab/2), and the sum-exponential sits exactly halfway "
            "between the two orders. Commuting generators collapse all "
            "three to one transport identically. Order is content, and "
            "the content is the commutator"),
        "times": [str(A_T), str(B_T)],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T2  refinement is free, reversal is not
# ----------------------------------------------------------------------


def certify_T2():
    n = 3
    K = madd(E(n, 0, 1), mscale(Fr(2), E(n, 1, 2)))   # nilpotent
    total = Fr(7, 4)
    partitions = (
        [total],
        [Fr(1), Fr(3, 4)],
        [Fr(1, 2)] * 3 + [Fr(1, 4)],
        [Fr(7, 20)] * 5,
    )
    W_ref = nilpotent_exp(total, K)
    for parts in partitions:
        assert sum(parts) == total
        W = ordered_product([(t, K) for t in parts])
        assert W == W_ref                             # refinement free

    # reversal of DISTINCT-generator segments is not free
    K1, K2 = E(n, 0, 1), E(n, 1, 2)
    a, b = Fr(1, 2), Fr(3)
    fwd = ordered_product([(a, K1), (b, K2)])
    rev = ordered_product([(b, K2), (a, K1)])
    assert msub(rev, fwd) == mscale(a * b, comm(K1, K2))
    assert fwd != rev

    return {
        "statement": (
            "Within one constant generator, every rational partition "
            "of the same total time yields the SAME transport — "
            "e^(K(a+b)) = e^(Ka) e^(Kb) exactly, so time-slicing is "
            "free. Reversing the order of two distinct-generator "
            "segments changes the transport by exactly ab[K1,K2]. "
            "Time-ORDERING, not time-slicing, is the load-bearing "
            "structure"),
        "partitions_checked": len(partitions),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T3  the mixed term obstructs transport (dynamic EMK-T1 T3)
# ----------------------------------------------------------------------

CH_UP = ([[Fr(0), Fr(1)], [Fr(0), Fr(0)]],
         [[Fr(0), Fr(2)], [Fr(0), Fr(0)]])            # EMK-T1's c1
CH_LO = ([[Fr(0), Fr(0)], [Fr(1), Fr(0)]],
         [[Fr(0), Fr(0)], [Fr(3), Fr(0)]])            # EMK-T1's c2


def group_loop(a, A, b, B):
    return mm(mm(nilpotent_exp(a, A), nilpotent_exp(b, B)),
              mm(nilpotent_exp(-a, A), nilpotent_exp(-b, B)))


def det2(M):
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


def certify_T3():
    # consume the pinned EMK-T1 channels and flatness facts
    assert t1.is_zero(t1.comm(CH_UP[0], CH_UP[1]))
    assert t1.is_zero(t1.comm(CH_LO[0], CH_LO[1]))
    mixed = t1.mixed_curvature([CH_UP, CH_LO])
    assert not t1.is_zero(mixed)

    a, b = Fr(1, 3), Fr(2, 5)
    # channel-internal loops close exactly
    assert group_loop(a, CH_UP[0], b, CH_UP[1]) == eye(2)
    assert group_loop(a, CH_LO[0], b, CH_LO[1]) == eye(2)

    # cross-channel loop: closed form, exact
    K_cross = comm(CH_UP[0], CH_LO[0])                # [E12, E21]
    assert K_cross == [[Fr(1), Fr(0)], [Fr(0), Fr(-1)]]
    L = group_loop(a, CH_UP[0], b, CH_LO[0])
    closed_form = [
        [1 + a * b + a * a * b * b, -a * a * b],
        [a * b * b, 1 - a * b],
    ]
    assert L == closed_form
    remainder = msub(msub(L, eye(2)), mscale(a * b, K_cross))
    assert remainder == [[a * a * b * b, -a * a * b],
                         [a * b * b, Fr(0)]]          # higher order only
    assert L != eye(2)
    assert det2(L) == 1

    # trying more rational times: never closes unless ab = 0
    for aa in (Fr(1), Fr(-1, 2), Fr(0)):
        for bb in (Fr(2), Fr(0), Fr(-3, 4)):
            Lx = group_loop(aa, CH_UP[0], bb, CH_LO[0])
            assert (Lx == eye(2)) == (aa * bb == 0)

    return {
        "statement": (
            "The pinned EMK-T1 channels are individually flat and "
            "their internal transport loops close EXACTLY to I; the "
            "cross-channel loop L(a,b) is given in closed form as a "
            "rational matrix identity with leading bilinear part "
            "EXACTLY ab[E12,E21], never closes unless ab = 0, and has "
            "det 1. The mixed term obstructs transport itself: a "
            "channel-by-channel audit is provably insufficient in "
            "dynamics, not only in the static curvature decomposition"),
        "cross_commutator": [[str(x) for x in r] for r in
                             comm(CH_UP[0], CH_LO[0])],
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# T4  same clock, different history, exact ledger (dynamic EMK-T1 T5)
# ----------------------------------------------------------------------


def certify_T4():
    n = 3
    K1, K2 = E(n, 0, 1), E(n, 1, 2)
    a, b = Fr(3, 2), Fr(1, 5)

    clock_A = a + b                                   # identical clocks
    clock_B = b + a
    assert clock_A == clock_B

    W_A = ordered_product([(b, K2), (a, K1)])         # K2 first
    W_B = ordered_product([(a, K1), (b, K2)])         # K1 first
    D = msub(W_A, W_B)
    assert D == mscale(a * b, comm(K1, K2))           # exact holonomy
    assert D != zeros(n)

    # temporal residue via the PINNED EMK-T1 functional (max-abs norm
    # of Theta_1 - transported Theta_0 - ledger), applied to the 2x2
    # corner blocks the discrepancy lives in
    blk = [[W_A[0][1], W_A[0][2]], [W_A[1][1], W_A[1][2]]]
    blk_t = [[W_B[0][1], W_B[0][2]], [W_B[1][1], W_B[1][2]]]
    no_ledger = [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
    res_open = t1.temporal_residue(blk, blk_t, no_ledger)
    assert res_open == a * b                          # nonzero, exact
    lawful = [[Fr(0), a * b], [Fr(0), Fr(0)]]         # THE discrepancy
    assert t1.temporal_residue(blk, blk_t, lawful) == 0

    # and at the transport level: the exact corrective loop closes it
    corrected = mm(nilpotent_exp(1, mscale(-a * b, comm(K1, K2))), W_A)
    assert corrected == W_B

    return {
        "statement": (
            "Two evolutions through the same generators for the same "
            "total time carry IDENTICAL clocks but transports "
            "differing by exactly ab[K1,K2]. The pinned EMK-T1 "
            "temporal-residue functional reads the open residue as "
            "exactly ab, the lawful ledger closes it to 0 identically, "
            "and the exact corrective transport maps one history onto "
            "the other. Clock equality is not time closure — now as a "
            "statement about transport, completing the static EMK-T1 "
            "T5 separation"),
        "clock": str(a + b),
        "open_residue": str(Fr(3, 2) * Fr(1, 5)),
        "verdict": "PASS",
    }


# ----------------------------------------------------------------------
# certificate assembly
# ----------------------------------------------------------------------


def build_certificate():
    return {
        "capsule": "EMK-T2: time-ordered transport — order is content",
        "consumes": {
            "EMK-T1": ("channel matrices, comm/is_zero/mixed_curvature "
                       "and the temporal_residue functional consumed "
                       "from the pinned module, not re-declared"),
            "RST-1": ("the exact nilpotent-exponential machinery of "
                      "its T4 flow-generator block (cited; this module "
                      "is self-contained to keep folder CI local)"),
        },
        "T1_order_is_content": certify_T1(),
        "T2_refinement_free_reversal_not": certify_T2(),
        "T3_mixed_term_obstructs_transport": certify_T3(),
        "T4_same_clock_different_history_ledger": certify_T4(),
        "claim_boundary": {
            "continuous_time": (
                "NOT CLAIMED: P exp ordered exponentials, Magnus "
                "expansions, non-nilpotent generators and refinement "
                "limits — the source (topology paper 5.4) states the "
                "general object; only the exact nilpotent sector is "
                "certified"),
            "convention": (
                "DECLARED: later transport multiplies on the left; "
                "reversing the convention swaps the roles of the two "
                "orders and negates the discrepancy"),
            "RH_K0_L0": "not touched",
            "yang_mills_quantum_gravity": "not touched",
        },
        "provenance": {
            "prior_executable_version": "NONE — first dynamic capsule "
                                        "of the tensor/time layer",
            "threads": (
                "mixed-term necessity now static (EMK-T1 T3) AND "
                "dynamic (T3 here); projection blindness now includes "
                "a transport form (T4 here) alongside UGD-1 T4, "
                "EMK-T1 T5 and RST-2 T1"),
        },
    }


def main():
    cert = build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    out = os.path.join(HERE, "EMKT2_RESULT.json")
    with open(out, "w") as f:
        f.write(payload)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    print("EMK-T2 certificate written:", out)
    print("sha256:", digest)


if __name__ == "__main__":
    main()
