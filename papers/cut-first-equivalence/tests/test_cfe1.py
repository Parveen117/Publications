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


cfe1 = _load("cfe1", "cfe1_cut_first_equivalence.py")


# ---------------- memoryless limit ----------------

def test_curvature_vanishes_exactly_at_equilibrium():
    assert cfe1.curvature_density(Fr(1)) == 0


def test_curvature_nonzero_off_equilibrium():
    for chi in (Fr(3, 5), Fr(4, 5), Fr(6, 5), Fr(7, 5)):
        assert cfe1.curvature_density(chi) != 0


def test_every_maxwell_loop_closes_at_equilibrium():
    max_abs, stokes = cfe1.maxwell_loops_closed(Fr(1))
    assert max_abs == 0 and stokes


def test_per_cell_stokes_holds_off_equilibrium():
    # circulation of each cell must equal density * area, exactly
    _, stokes = cfe1.maxwell_loops_closed(Fr(7, 5))
    assert stokes


# ---------------- residue identity ----------------

def test_stokes_residue_identity_exact():
    p_lo, p_hi, v_lo, v_hi = cfe1.LOOP_RECT
    for chi in (Fr(3, 5), Fr(7, 5), Fr(11, 10), Fr(2)):
        circ = cfe1.circulation(
            cfe1.rect_boundary(p_lo, p_hi, v_lo, v_hi), chi)
        flux = cfe1.curvature_flux(p_lo, p_hi, v_lo, v_hi, chi)
        assert circ == flux


def test_residue_is_a_rational_not_a_bracket():
    p_lo, p_hi, v_lo, v_hi = cfe1.LOOP_RECT
    circ = cfe1.circulation(
        cfe1.rect_boundary(p_lo, p_hi, v_lo, v_hi), Fr(7, 5))
    assert isinstance(circ, Fr)


# ---------------- obstruction faithfulness ----------------

def test_residue_zero_iff_equilibrium():
    p_lo, p_hi, v_lo, v_hi = cfe1.LOOP_RECT
    loop = cfe1.rect_boundary(p_lo, p_hi, v_lo, v_hi)
    assert cfe1.circulation(loop, Fr(1)) == 0
    assert cfe1.circulation(loop, Fr(6, 5)) != 0


def test_residue_strictly_monotone_in_memory_dial():
    p_lo, p_hi, v_lo, v_hi = cfe1.LOOP_RECT
    loop = cfe1.rect_boundary(p_lo, p_hi, v_lo, v_hi)
    chis = [Fr(3, 5), Fr(4, 5), Fr(1), Fr(6, 5), Fr(7, 5)]
    vals = [cfe1.circulation(loop, c) for c in chis]
    assert all(vals[k] < vals[k + 1] for k in range(len(vals) - 1))


def test_residue_sign_tracks_dial_sign():
    p_lo, p_hi, v_lo, v_hi = cfe1.LOOP_RECT
    loop = cfe1.rect_boundary(p_lo, p_hi, v_lo, v_hi)
    assert cfe1.circulation(loop, Fr(3, 5)) < 0   # mu < 0
    assert cfe1.circulation(loop, Fr(7, 5)) > 0   # mu > 0


# ---------------- equilibrium invariant ----------------

def test_invariant_is_one_at_equilibrium():
    _, _, inv = cfe1.invariants(Fr(1))
    assert inv == 1


def test_invariant_departs_off_equilibrium():
    for chi in (Fr(4, 5), Fr(6, 5), Fr(7, 5)):
        _, _, inv = cfe1.invariants(chi)
        assert inv != 1


# ---------------- certificate integrity + honest boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = cfe1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_CFE1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_open_obligations_are_declared_open():
    cert = json.load(open(os.path.join(CERT_DIR, "CFE1_RESULT.json")))
    ob = cert["open_obligations_for_full_theorem"]
    assert ob["S_surjectivity_completeness"].startswith("OPEN")
    assert ob["U_uniqueness"].startswith("OPEN")
    cb = cert["claim_boundary"]
    assert cb["surjectivity_S"] == "OPEN"
    assert cb["uniqueness_U"] == "OPEN"
    assert cb["full_generality_beyond_witness_EOS"] == "OPEN"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["YM_continuum_gates"] == "not touched"


def test_target_theorem_stated_but_not_claimed_proven():
    cert = json.load(open(os.path.join(CERT_DIR, "CFE1_RESULT.json")))
    assert "Cut-First Equivalence" in cert["target_theorem"]["name"]
    # the claim_status must say 'core', never 'proven in full'
    assert "core" in cert["claim_status"]
    assert "OPEN" in cert["claim_status"]
