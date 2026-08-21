import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CERT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ug = _load("ugdg1", "ugdg1_phase_scale_seam_geometry.py")


# ---------------- T1: heterogeneous closure ----------------

def test_three_sectors_have_three_rules():
    assert ug.phase_closes(ug.K) and ug.phase_closes(0)
    assert not ug.phase_closes(1)
    assert ug.scale_closes(Fr(0)) and not ug.scale_closes(Fr(ug.K))
    assert ug.sheet_closes(0) and not ug.sheet_closes(ug.K)


def test_mod_K_audit_wrongly_passes_an_open_sheet():
    open_state = (0, Fr(0), ug.K)
    assert not ug.ugd_closes(open_state)
    collapsed = (open_state[0] % ug.K == 0
                 and open_state[1] % ug.K == 0
                 and open_state[2] % ug.K == 0)
    assert collapsed


def test_exact_audit_wrongly_rejects_a_phase_return():
    good = (ug.K, Fr(0), 0)
    assert ug.ugd_closes(good)
    assert not all(x == 0 for x in good)


def test_full_rule_is_correct_on_both_witnesses():
    assert ug.ugd_closes((ug.K, Fr(0), 0))
    assert not ug.ugd_closes((0, Fr(0), ug.K))
    assert not ug.ugd_closes((0, Fr(1, 3), 0))
    assert not ug.ugd_closes((1, Fr(0), 0))


# ---------------- T2: the cocycle ----------------

def test_sheet_cocycle_is_invariant_under_every_chart_rechoice():
    n = (4, -1, 2)
    base = ug.sheet_cocycle(*n)
    for m_a, m_b, m_c in product(range(-4, 5), repeat=3):
        new = (n[0] + m_a - m_b, n[1] + m_b - m_c, n[2] + m_c - m_a)
        assert ug.sheet_cocycle(*new) == base


def test_the_algebraic_identity_behind_the_invariance():
    for m_a, m_b, m_c in ((7, -3, 11), (0, 0, 0), (-5, -5, 2),
                          (1, 2, 3)):
        assert (m_a - m_b) + (m_b - m_c) + (m_c - m_a) == 0


def test_nonzero_defect_cannot_be_erased_zero_stays_zero():
    assert ug.sheet_cocycle(4, -1, 2) == 5 != 0
    assert ug.sheet_cocycle(1, 1, -2) == 0


def test_phase_cocycle_also_invariant_mod_K():
    chi = (2, 3, 2)
    base = ug.phase_cocycle(*chi)
    for m_a, m_b, m_c in product(range(ug.K), repeat=3):
        new = (chi[0] + m_a - m_b, chi[1] + m_b - m_c,
               chi[2] + m_c - m_a)
        assert ug.phase_cocycle(*new) == base


def test_scale_cocycle_closes_with_the_constructed_chart():
    f_ab = ug.affine(Fr(2), Fr(1))
    f_bc = ug.affine(Fr(1, 2), Fr(-3))
    inner = ug.compose(f_bc, f_ab)
    a_i, b_i = inner
    f_ca = (Fr(1) / a_i, -b_i / a_i)
    assert ug.compose(f_ca, inner) == ug.IDENTITY
    assert ug.compose(ug.affine(Fr(3), Fr(0)), inner) != ug.IDENTITY


# ---------------- T3: flat but discrete memory ----------------

def test_each_constant_coefficient_form_is_flat():
    zero = ug.g3.bp({})
    for c in (Fr(3, 5), Fr(2, 25), Fr(0)):
        assert ug.g3.abelian_curvature(ug.g3.bp({(0, 0): c}), zero) == {}


def test_consumed_flat_with_period_from_pinned_emkg3():
    zero = ug.g3.bp({})
    a = ug.g3.bp({(0, 0): ug.g3.P_PHI})
    assert ug.g3.abelian_curvature(a, zero) == {}
    assert ug.g3.bp_int_du(a, Fr(0), ug.g3.L, Fr(0)) == ug.g3.ALPHA != 0


def test_relaxing_the_memory_sector_cannot_change_the_integer_class():
    zero = ug.g3.bp({})
    relaxed = ug.g3.bp({(0, 1): Fr(7)})
    assert ug.g3.abelian_curvature(relaxed, zero) != {}
    assert ug.sheet_cocycle(4, -1, 2) == 5


# ---------------- T4: non-implications ----------------

def test_base_closure_does_not_force_ugd_closure():
    w = {"base": Fr(0), "levi_civita": Fr(0), "ugd": Fr(3),
         "rtc_open": Fr(0)}
    c = ug.closed_sectors(w)
    assert c["base"] and not c["ugd"]


def test_ugd_closure_does_not_force_rtc_closure():
    w = {"base": Fr(0), "levi_civita": Fr(0), "ugd": Fr(0),
         "rtc_open": Fr(5, 2)}
    c = ug.closed_sectors(w)
    assert c["ugd"] and not c["rtc_open"]


def test_all_ordered_pairs_are_independent():
    count = 0
    for s in ug.SECTORS4:
        for t in ug.SECTORS4:
            if s == t:
                continue
            st = {x: (Fr(0) if x == s else Fr(1)) for x in ug.SECTORS4}
            c = ug.closed_sectors(st)
            assert c[s] and not c[t]
            count += 1
    assert count == 12


# ---------------- T5: faithful vs lossy ----------------

def test_faithful_projection_preserves_and_reflects():
    closed = {"phase": ug.K, "scale": Fr(0), "sheet": 0}
    op = {"phase": 0, "scale": Fr(0), "sheet": 3}
    assert ug.ugd_state_closed(ug.project_faithful(closed))
    assert not ug.ugd_state_closed(ug.project_faithful(op))


def test_lossy_projection_preserves_but_does_not_reflect():
    op = {"phase": 0, "scale": Fr(0), "sheet": 3}
    assert not ug.ugd_state_closed(op)
    assert ug.ugd_state_closed(ug.project_lossy(op))


def test_identifying_projection_discards_nothing_yet_fails():
    tricky = {"phase": 3, "scale": Fr(0), "sheet": -3}
    assert not ug.ugd_state_closed(tricky)
    m = ug.project_identifying(tricky)
    assert m["merged"] == 0 and m["scale"] == 0 and m["sheet"] == 0


def test_injectivity_is_what_separates_them():
    states = [{"phase": p, "scale": Fr(s), "sheet": k}
              for p in (0, 3) for s in (0, 1) for k in (-3, 0, 3)]
    faithful = {tuple(sorted(ug.project_faithful(x).items()))
                for x in states}
    lossy = {tuple(sorted(ug.project_lossy(x).items()))
             for x in states}
    assert len(faithful) == len(states)
    assert len(lossy) < len(states)


# ---------------- T6: the bridge ----------------

def test_numerals_share_visible_readings_but_differ_in_seam():
    Kn = ug.u1.K_DEFAULT
    clean = ug.u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    seamed = ug.u1.Numeral({0: (1, 1), 1: (1, -1)}, 2, Kn)
    assert clean.classical_projection() == seamed.classical_projection()
    assert clean.phase_content() == seamed.phase_content()
    assert clean.total_seam() != seamed.total_seam()


def test_transition_residue_is_visible_zero_sheet_nonzero():
    Kn = ug.u1.K_DEFAULT
    clean = ug.u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    seamed = ug.u1.Numeral({0: (1, 1), 1: (1, -1)}, 2, Kn)
    res = ug.ugd_residue(clean, seamed)
    assert res["phase"] == 0 and res["scale"] == 0
    assert res["sheet"] == seamed.total_seam() - clean.total_seam() == 2
    assert not ug.ugd_state_closed(res)


def test_lossy_projection_reads_the_memory_bearing_transition_closed():
    Kn = ug.u1.K_DEFAULT
    clean = ug.u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    seamed = ug.u1.Numeral({0: (1, 1), 1: (1, -1)}, 2, Kn)
    res = ug.ugd_residue(clean, seamed)
    assert ug.ugd_state_closed(ug.project_lossy(res))


def test_null_seam_transition_genuinely_closes():
    Kn = ug.u1.K_DEFAULT
    clean = ug.u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    assert ug.ugd_state_closed(ug.ugd_residue(clean, clean))


def test_visible_difference_stays_open_after_the_lossy_projection():
    Kn = ug.u1.K_DEFAULT
    clean = ug.u1.Numeral({0: (1, 0), 1: (1, 0)}, 0, Kn)
    other = ug.u1.Numeral({0: (2, 0), 1: (1, 0)}, 0, Kn)
    res = ug.ugd_residue(clean, other)
    assert res["scale"] != 0
    assert not ug.ugd_state_closed(ug.project_lossy(res))


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = ug.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_UGDG1.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_declares_dictionary_and_reduction():
    cert = json.load(open(os.path.join(CERT_DIR, "UGDG1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["declared_structure"].startswith("DECLARED")
    assert "never evaluated" in cb["no_logarithm_evaluated"]
    assert "exponent-mod-K" in cb["phase_convention"]
    assert cb["not_claimed"].startswith("NOT CLAIMED")
    assert "REDUCTION" in cb["reduction_not_replacement"]
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "ugdg1_phase_scale_seam_geometry.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "UGDG1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
