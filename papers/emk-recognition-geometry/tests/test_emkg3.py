import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CERT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


g3 = _load("emkg3", "emkg3_helical_sheet_memory.py")


# ---------------- T1: the mapping torus ----------------

def test_transition_is_a_bijection_with_exact_inverse():
    for f in g3.FIBRES:
        assert g3.rho_inv(g3.rho(f)) == f
        assert g3.rho(g3.rho_inv(f)) == f
    assert len({g3.rho(f) for f in g3.FIBRES}) == len(g3.FIBRES)


def test_monodromy_after_n_circuits_is_rho_to_the_n():
    for f in g3.FIBRES[:6]:
        acc = f
        for n in range(1, 6):
            acc = g3.rho(acc)
            assert acc == g3.rho_power(f, n)
        assert g3.rho_power(f, 0) == f
        assert g3.rho_power(g3.rho_power(f, 3), -3) == f


def test_deck_action_is_free_and_additive():
    pts = [(Fr(u), Fr(v), f) for u in (0, 1) for v in (0, 2)
           for f in g3.FIBRES[:4]]
    for p in pts:
        for n in range(-3, 4):
            if n != 0:
                assert g3.deck(p, n) != p
        for m in (-2, 1):
            for n in (-1, 3):
                assert g3.deck(g3.deck(p, m), n) == g3.deck(p, m + n)


def test_each_orbit_meets_the_fundamental_domain_once():
    pts = [(Fr(u), Fr(0), g3.FIBRES[0]) for u in (0, 1, 3)]
    for p in pts:
        hits = [n for n in range(-5, 6)
                if Fr(0) <= g3.deck(p, n)[0] < g3.L]
        assert len(hits) == 1


# ---------------- T2: pitch ----------------

def test_compensated_coordinates_are_exactly_invariant():
    for u in (Fr(0), Fr(2, 7), Fr(-3)):
        for phi in (Fr(0), Fr(4, 3)):
            for sigma in (Fr(0), Fr(-5, 2)):
                assert g3.compensated(u + g3.L, phi + g3.ALPHA,
                                      sigma + g3.BETA) == \
                    g3.compensated(u, phi, sigma)


def test_pitch_times_period_recovers_the_advances():
    assert g3.P_PHI * g3.L == g3.ALPHA
    assert g3.P_SIGMA * g3.L == g3.BETA


def test_wrong_pitch_fails_by_an_exact_amount():
    wrong = g3.P_PHI + Fr(1, 4)
    fail = g3.ALPHA - wrong * g3.L
    assert fail == -Fr(5, 4) != 0


# ---------------- T3: flat is not memoryless ----------------

def test_constant_pitch_connection_is_exactly_flat():
    a = g3.bp({(0, 0): g3.P_PHI})
    assert g3.abelian_curvature(a, g3.bp({})) == {}


def test_every_contractible_rectangle_has_zero_holonomy():
    a, b = g3.bp({(0, 0): g3.P_PHI}), g3.bp({})
    for r in ((Fr(0), Fr(2), Fr(-1), Fr(1)),
              (Fr(-3), Fr(1), Fr(0), Fr(5, 2)),
              (Fr(1, 3), Fr(4), Fr(-2), Fr(1, 2))):
        assert g3.loop_holonomy(a, b, *r) == 0


def test_noncontractible_circuit_holonomy_is_exactly_alpha():
    a = g3.bp({(0, 0): g3.P_PHI})
    for v in (Fr(0), Fr(2), Fr(-1, 3)):
        assert g3.bp_int_du(a, Fr(0), g3.L, v) == g3.ALPHA != 0
    for n in (2, 5, -3):
        assert g3.bp_int_du(a, Fr(0), n * g3.L, Fr(0)) == n * g3.ALPHA


def test_curved_connection_holonomy_equals_exact_flux():
    a = g3.bp({(0, 0): g3.P_PHI, (0, 1): Fr(1)})
    b = g3.bp({})
    assert g3.abelian_curvature(a, b) == {(0, 0): Fr(-1)}
    u0, u1, v0, v1 = Fr(0), Fr(2), Fr(-1), Fr(1)
    assert g3.loop_holonomy(a, b, u0, u1, v0, v1) == \
        Fr(-1) * (u1 - u0) * (v1 - v0)


# ---------------- T4: visible vs lifted return ----------------

def test_lifted_return_always_implies_visible_return():
    for al in (0, 2, 3, 6):
        for be in (Fr(0), Fr(2, 5)):
            for q in (0, 1, -2):
                if g3.lifted_returns(al, be, q):
                    assert g3.visible_returns(al)


def test_lifted_return_exactly_when_all_three_close():
    for al in (0, 2, 3, 6, 12):
        for be in (Fr(0), Fr(2, 5), Fr(-1)):
            for q in (0, 1, -2):
                assert g3.lifted_returns(al, be, q) == \
                    (al % g3.KMOD == 0 and be == 0 and q == 0)


def test_headline_witness_phase_closes_sheet_does_not():
    assert g3.visible_returns(g3.KMOD)
    assert not g3.lifted_returns(g3.KMOD, g3.BETA, g3.Q)
    assert g3.visible_returns(0) and not g3.lifted_returns(0, Fr(0), 2)


def test_sheet_increment_never_undone_by_more_circuits():
    f0 = (0, Fr(0), 0)
    for n in range(-6, 7):
        if n != 0:
            assert g3.rho_power(f0, n)[2] == n * g3.Q != 0


# ---------------- T5: the four classes ----------------

def test_each_class_has_its_witness():
    adm = {"sheet": True}
    cases = {
        "exact_return": ((0, Fr(0), 0, Fr(0)), (0, Fr(0), 0, Fr(0))),
        "lawfully_transported_return": (
            (2, Fr(3, 4), 0, Fr(0)), (2, Fr(3, 4), 0, Fr(0))),
        "memory_bearing_return": (
            (0, Fr(0), 3, Fr(0)), (0, Fr(0), 3, Fr(0))),
        "open_obstruction": (
            (0, Fr(0), 3, Fr(1, 2)), (0, Fr(0), 3, Fr(0))),
    }
    for expected, (res, led) in cases.items():
        assert g3.classify_return(res, led, adm) == expected
    assert len({g3.classify_return(r, l, adm)
                for (r, l) in cases.values()}) == 4


def test_memory_bearing_is_not_exact_return():
    adm = {"sheet": True}
    res, led = (0, Fr(0), 3, Fr(0)), (0, Fr(0), 3, Fr(0))
    assert res[2] != 0
    assert g3.classify_return(res, led, adm) != "exact_return"


def test_removing_the_ledger_opens_the_obstruction():
    adm = {"sheet": True}
    res = (0, Fr(0), 3, Fr(0))
    assert g3.classify_return(res, (0, Fr(0), 0, Fr(0)), adm) == \
        "open_obstruction"


def test_deadmitting_the_sheet_reclassifies():
    res, led = (0, Fr(0), 3, Fr(0)), (0, Fr(0), 3, Fr(0))
    assert g3.classify_return(res, led, {"sheet": False}) == \
        "lawfully_transported_return"


# ---------------- T6: no proper subset forces closure ----------------

def test_pinned_emkg2_seam_holonomy_is_zero():
    A = g3.g2.A_rational(Fr(3))
    assert g3.g2.holonomy_boundary(A, Fr(0), g3.L, Fr(0), Fr(0)) == 0


def test_flat_and_lc_trivial_yet_transition_open():
    res = {"base_point": Fr(0), "levi_civita": Fr(0),
           "helical": Fr(g3.Q), "rtc_curvature": Fr(0),
           "thermodynamic": Fr(0)}
    allact = {s: True for s in g3.SECTORS}
    assert g3.abelian_curvature(g3.bp({(0, 0): g3.P_PHI}),
                                g3.bp({})) == {}
    assert not g3.closed_on(res, allact)


def test_every_proper_subset_closes_without_forcing_closure():
    allact = {s: True for s in g3.SECTORS}
    count = 0
    for r in range(len(g3.SECTORS)):
        for S in combinations(g3.SECTORS, r):
            res = {s: (Fr(0) if s in S else Fr(1)) for s in g3.SECTORS}
            assert g3.closed_on(res, {s: (s in S) for s in g3.SECTORS})
            assert not g3.closed_on(res, allact)
            count += 1
    assert count == 2 ** len(g3.SECTORS) - 1


def test_nonfaithful_projection_closes_an_open_lift():
    res = {"base_point": Fr(0), "levi_civita": Fr(0),
           "helical": Fr(g3.Q), "rtc_curvature": Fr(0),
           "thermodynamic": Fr(0)}
    allact = {s: True for s in g3.SECTORS}
    dropped = dict(res)
    dropped["helical"] = Fr(0)
    assert not g3.closed_on(res, allact)
    assert g3.closed_on(dropped, allact)


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = g3.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKG3.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_declares_and_refuses():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG3_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["declared_structure"].startswith("DECLARED")
    assert cb["not_certified"].startswith("NOT CERTIFIED")
    assert cb["physical_identification"].startswith("NOT CLAIMED")
    assert "exponent-mod-K" in cb["phase_convention"]
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "emkg3_helical_sheet_memory.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG3_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
