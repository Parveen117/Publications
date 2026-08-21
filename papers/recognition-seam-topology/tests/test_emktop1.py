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


tp = _load("emktop1", "emktop1_vault_topology_appendix.py")


# ---------------- T1: KIR holonomy ----------------

def test_alphabet_relations_exact():
    assert tp.cmmul(tp.R, tp.R) == tp.cmscale(tp.GNEG, tp.I2)
    assert tp.cmmul(tp.K, tp.K) == tp.I2
    assert tp.cmmul(tp.R, tp.K) == tp.cmscale(
        tp.GNEG, tp.cmmul(tp.K, tp.R))


def test_operator_system_is_closed_of_order_eight():
    seen, frontier = set(), [tp.I2]
    while frontier:
        W = frontier.pop()
        key = tuple(tuple(r) for r in W)
        if key in seen:
            continue
        seen.add(key)
        for G in (tp.R, tp.K):
            frontier.append(tp.cmmul(W, G))
    assert len(seen) == 8


def test_visibly_closed_word_with_nontrivial_holonomy():
    w = tp.word("RR")
    assert w == tp.cmscale(tp.GNEG, tp.I2)
    assert tp.shadow(w) == tp.shadow(tp.I2)
    assert w != tp.I2


def test_word_closing_in_both_sectors():
    assert tp.word("RRRR") == tp.I2
    assert tp.shadow(tp.word("RRRR")) == tp.shadow(tp.I2)


def test_mixed_word_squares_to_identity():
    assert tp.cmmul(tp.RK, tp.RK) == tp.I2
    assert tp.shadow(tp.RK) != tp.shadow(tp.I2)


def test_ledger_closes_the_open_holonomy():
    w = tp.word("RR")
    assert tp.cmmul(tp.cmscale(tp.GNEG, tp.I2), w) == tp.I2


# ---------------- T2: cochains and seam charge ----------------

def test_coboundary_composition_vanishes():
    for v in tp.VERTS:
        f = {u: (Fr(1) if u == v else Fr(0)) for u in tp.VERTS}
        assert tp.d1(tp.d0(f)) == [Fr(0)] * len(tp.FACES)


def test_cycles_are_genuine_cycles():
    for z in ([Fr(1), Fr(1), Fr(1), Fr(1), Fr(0)],
              [Fr(0), Fr(0), Fr(1), Fr(1), Fr(1)],
              [Fr(1), Fr(1), Fr(0), Fr(0), Fr(-1)]):
        assert all(v == 0 for v in tp.boundary(z).values())


def test_closed_cochain_charge_equal_on_homologous_cycles():
    alpha = [Fr(2), Fr(-1), Fr(3), Fr(1, 2), Fr(1)]
    assert tp.d1(alpha) == [Fr(0)]
    z1 = [Fr(1), Fr(1), Fr(1), Fr(1), Fr(0)]
    z2 = [Fr(0), Fr(0), Fr(1), Fr(1), Fr(1)]
    assert tp.pair(alpha, z1) == tp.pair(alpha, z2) != 0


def test_coboundary_has_zero_charge_on_every_cycle():
    beta = tp.d0({0: Fr(1), 1: Fr(-4), 2: Fr(2, 3), 3: Fr(9)})
    for z in ([Fr(1), Fr(1), Fr(1), Fr(1), Fr(0)],
              [Fr(0), Fr(0), Fr(1), Fr(1), Fr(1)]):
        assert tp.pair(beta, z) == 0


def test_nonclosed_cochain_is_not_a_class_function():
    gamma = [Fr(1), Fr(0), Fr(0), Fr(0), Fr(0)]
    assert tp.d1(gamma) != [Fr(0)]
    z1 = [Fr(1), Fr(1), Fr(1), Fr(1), Fr(0)]
    z2 = [Fr(0), Fr(0), Fr(1), Fr(1), Fr(1)]
    assert tp.pair(gamma, z1) != tp.pair(gamma, z2)


# ---------------- T3: determinant charge ----------------

def test_global_determinant_return_with_nonzero_subcycles():
    charges = (Fr(2), Fr(3), Fr(1, 2), Fr(1, 3))
    total = Fr(1)
    for c in charges:
        total *= c
    assert total == 1
    assert charges[0] * charges[1] == 6
    assert charges[2] * charges[3] == Fr(1, 6)


def test_refinement_preserves_the_loop_charge():
    charges = (Fr(2), Fr(3), Fr(1, 2), Fr(1, 3))
    refined = (Fr(4), Fr(1, 2)) + charges[1:]
    t1 = t2 = Fr(1)
    for c in charges:
        t1 *= c
    for c in refined:
        t2 *= c
    assert t1 == t2 == 1


def test_charge_is_realized_by_actual_matrices():
    M0 = [[Fr(1), Fr(0)], [Fr(0), Fr(1)]]
    M1 = [[Fr(2), Fr(0)], [Fr(0), Fr(1)]]
    assert tp.det2(M1) / tp.det2(M0) == Fr(2)


def test_ledger_neutralizes_the_determinant_charge():
    q = Fr(6)
    assert q * (Fr(1) / q) == 1


# ---------------- T4: the trichotomy ----------------

def test_three_witnesses_get_three_classes():
    exact = tp.d0({0: Fr(0), 1: Fr(3), 2: Fr(-1), 3: Fr(4)})
    closed_ne = [Fr(2), Fr(-1), Fr(3), Fr(1, 2), Fr(1)]
    open_r = [Fr(1), Fr(0), Fr(0), Fr(0), Fr(0)]
    assert tp.classify(exact) == "exact"
    assert tp.classify(closed_ne) == "closed_not_exact"
    assert tp.classify(open_r) == "open"


def test_closed_nonexact_survives_every_coboundary_correction():
    closed_ne = [Fr(2), Fr(-1), Fr(3), Fr(1, 2), Fr(1)]
    z = [Fr(1), Fr(1), Fr(1), Fr(1), Fr(0)]
    for shift in ({0: Fr(1), 1: Fr(0), 2: Fr(0), 3: Fr(0)},
                  {0: Fr(0), 1: Fr(-3), 2: Fr(8), 3: Fr(1, 5)}):
        corrected = [a - b for a, b in zip(closed_ne, tp.d0(shift))]
        assert tp.classify(corrected) == "closed_not_exact"
        assert tp.pair(corrected, z) == tp.pair(closed_ne, z)


def test_exact_residue_is_removed_by_its_own_beta():
    beta = {0: Fr(0), 1: Fr(3), 2: Fr(-1), 3: Fr(4)}
    exact = tp.d0(beta)
    assert [a - b for a, b in zip(exact, tp.d0(beta))] == [Fr(0)] * 5


# ---------------- T5: sheet rotation ----------------

def test_visible_phase_return_with_nonzero_sheet_class():
    st = (0, 0, 0)
    for (a, n) in ((3, 1), (3, 1)):
        st = tp.sheet_act(st, a, n)
    assert st[0] == 0
    assert (st[1], st[2]) == (2, 6)


def test_full_closure_control():
    st = tp.sheet_act(tp.sheet_act((0, 0, 0), 2, 0), 4, 0)
    assert st == (0, 0, 0)


def test_phase_stays_an_integer_exponent_mod_K():
    st = (0, 0, 0)
    for (a, n) in ((5, 2), (4, -3), (1, 1)):
        st = tp.sheet_act(st, a, n)
        assert isinstance(st[0], int) and 0 <= st[0] < tp.KMOD


# ---------------- T6: seam survival ----------------

def test_shadow_survival_only_when_other_sectors_open():
    r = {"visible": Fr(0), "KIR": Fr(1), "determinant": Fr(0),
         "projection": Fr(0), "memory": Fr(3, 4)}
    allact = {s: True for s in tp.SECTORS}
    assert not tp.survives(r, allact)
    assert tp.survives(r, {s: (s == "visible") for s in tp.SECTORS})


def test_exact_ledger_closes_and_wrong_ledger_does_not():
    r = {"visible": Fr(0), "KIR": Fr(1), "determinant": Fr(0),
         "projection": Fr(0), "memory": Fr(3, 4)}
    allact = {s: True for s in tp.SECTORS}
    assert tp.survives(r, allact, {"KIR": Fr(1), "memory": Fr(3, 4)})
    assert not tp.survives(r, allact, {"KIR": Fr(1)})


def test_deactivating_a_sector_only_helps():
    r = {"visible": Fr(0), "KIR": Fr(1), "determinant": Fr(0),
         "projection": Fr(0), "memory": Fr(3, 4)}
    allact = {s: True for s in tp.SECTORS}
    for s in tp.SECTORS:
        reduced = dict(allact)
        reduced[s] = False
        assert tp.survives(r, allact) <= tp.survives(r, reduced)


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = tp.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKTOP1.sha256")) as f:
        assert digest == f.read().strip()


def test_source_records_distinctness_from_rst1():
    cert = json.load(open(os.path.join(CERT_DIR,
                                       "EMKTOP1_RESULT.json")))
    assert "RST-1 certified the arXiv topology PAPER" in \
        cert["source"]["distinct_from_RST1"]


def test_claim_boundary_declares_protocol():
    cert = json.load(open(os.path.join(CERT_DIR,
                                       "EMKTOP1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["protocol_declared_not_certified"].startswith("DECLARED")
    assert cb["continuum_representative"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "emktop1_vault_topology_appendix.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR,
                                       "EMKTOP1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
