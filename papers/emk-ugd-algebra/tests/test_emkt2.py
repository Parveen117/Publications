import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CERT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


t2 = _load("emkt2", "emkt2_time_ordered_transport.py")


# ---------------- T1: order is content ----------------

def test_nilpotent_exp_is_exact_group_homomorphism_in_time():
    K = t2.madd(t2.E(3, 0, 1), t2.E(3, 1, 2))
    for a in (Fr(1, 2), Fr(-2), Fr(3, 7)):
        for b in (Fr(1), Fr(-1, 3)):
            lhs = t2.nilpotent_exp(a + b, K)
            rhs = t2.mm(t2.nilpotent_exp(a, K), t2.nilpotent_exp(b, K))
            assert lhs == rhs


def test_two_orders_and_sum_exponential_all_differ_by_commutator():
    n = 3
    K1, K2 = t2.E(n, 0, 1), t2.E(n, 1, 2)
    C = t2.comm(K1, K2)
    for a in (Fr(1), Fr(-1, 2)):
        for b in (Fr(2), Fr(1, 3)):
            W12 = t2.ordered_product([(a, K1), (b, K2)])
            W21 = t2.ordered_product([(b, K2), (a, K1)])
            Ws = t2.nilpotent_exp(
                1, t2.madd(t2.mscale(a, K1), t2.mscale(b, K2)))
            assert t2.msub(W21, W12) == t2.mscale(a * b, C)
            assert t2.mscale(Fr(2), Ws) == t2.madd(W12, W21)


def test_commuting_generators_collapse_all_three():
    n = 3
    K1, K3 = t2.E(n, 0, 1), t2.E(n, 0, 2)
    assert t2.comm(K1, K3) == t2.zeros(n)
    a, b = Fr(5, 4), Fr(-2, 3)
    U = t2.ordered_product([(a, K1), (b, K3)])
    V = t2.ordered_product([(b, K3), (a, K1)])
    S = t2.nilpotent_exp(1, t2.madd(t2.mscale(a, K1), t2.mscale(b, K3)))
    assert U == V == S


# ---------------- T2: refinement vs reversal ----------------

def test_partition_invariance_within_constant_generator():
    K = t2.madd(t2.E(3, 0, 1), t2.mscale(Fr(2), t2.E(3, 1, 2)))
    total = Fr(3, 2)
    ref = t2.nilpotent_exp(total, K)
    for parts in ([total], [Fr(1), Fr(1, 2)], [Fr(1, 2)] * 3,
                  [Fr(3, 10)] * 5):
        assert sum(parts) == total
        assert t2.ordered_product([(p, K) for p in parts]) == ref


def test_reversal_changes_transport_by_exact_commutator():
    n = 3
    K1, K2 = t2.E(n, 0, 1), t2.E(n, 1, 2)
    a, b = Fr(1, 2), Fr(3)
    fwd = t2.ordered_product([(a, K1), (b, K2)])
    rev = t2.ordered_product([(b, K2), (a, K1)])
    assert t2.msub(rev, fwd) == t2.mscale(a * b, t2.comm(K1, K2))


# ---------------- T3: mixed term obstructs transport ----------------

def test_channels_are_the_pinned_emkt1_channels_and_flat():
    assert t2.t1.is_zero(t2.t1.comm(t2.CH_UP[0], t2.CH_UP[1]))
    assert t2.t1.is_zero(t2.t1.comm(t2.CH_LO[0], t2.CH_LO[1]))
    assert not t2.t1.is_zero(
        t2.t1.mixed_curvature([t2.CH_UP, t2.CH_LO]))


def test_channel_internal_loops_close_exactly():
    for a in (Fr(1), Fr(-2, 3)):
        for b in (Fr(1, 5), Fr(4)):
            assert t2.group_loop(a, t2.CH_UP[0], b, t2.CH_UP[1]) == \
                t2.eye(2)
            assert t2.group_loop(a, t2.CH_LO[0], b, t2.CH_LO[1]) == \
                t2.eye(2)


def test_cross_channel_loop_closed_form_and_leading_term():
    Kc = t2.comm(t2.CH_UP[0], t2.CH_LO[0])
    assert Kc == [[Fr(1), Fr(0)], [Fr(0), Fr(-1)]]
    for a in (Fr(1, 3), Fr(2)):
        for b in (Fr(2, 5), Fr(-1)):
            L = t2.group_loop(a, t2.CH_UP[0], b, t2.CH_LO[0])
            assert L == [[1 + a * b + a * a * b * b, -a * a * b],
                         [a * b * b, 1 - a * b]]
            assert t2.det2(L) == 1
            assert (L == t2.eye(2)) == (a * b == 0)


def test_cross_loop_closes_only_when_a_time_vanishes():
    assert t2.group_loop(Fr(0), t2.CH_UP[0], Fr(7), t2.CH_LO[0]) == \
        t2.eye(2)
    assert t2.group_loop(Fr(7), t2.CH_UP[0], Fr(0), t2.CH_LO[0]) == \
        t2.eye(2)


# ---------------- T4: same clock, different history ----------------

def test_identical_clock_different_transport_exact_holonomy():
    n = 3
    K1, K2 = t2.E(n, 0, 1), t2.E(n, 1, 2)
    a, b = Fr(3, 2), Fr(1, 5)
    W_A = t2.ordered_product([(b, K2), (a, K1)])
    W_B = t2.ordered_product([(a, K1), (b, K2)])
    assert a + b == b + a                            # clocks agree
    assert t2.msub(W_A, W_B) == t2.mscale(a * b, t2.comm(K1, K2))


def test_pinned_temporal_residue_reads_ab_then_ledger_closes():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKT2_RESULT.json")))
    t4 = cert["T4_same_clock_different_history_ledger"]
    assert t4["open_residue"] == str(Fr(3, 2) * Fr(1, 5))
    assert t4["verdict"] == "PASS"


def test_corrective_transport_maps_history_onto_history():
    n = 3
    K1, K2 = t2.E(n, 0, 1), t2.E(n, 1, 2)
    a, b = Fr(3, 2), Fr(1, 5)
    W_A = t2.ordered_product([(b, K2), (a, K1)])
    W_B = t2.ordered_product([(a, K1), (b, K2)])
    corr = t2.mm(t2.nilpotent_exp(
        1, t2.mscale(-a * b, t2.comm(K1, K2))), W_A)
    assert corr == W_B


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = t2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKT2.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_continuous_time_not_claimed():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKT2_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["continuous_time"].startswith("NOT CLAIMED")
    assert cb["convention"].startswith("DECLARED")
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "emkt2_time_ordered_transport.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKT2_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 4
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
