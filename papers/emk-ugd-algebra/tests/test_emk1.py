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


emk1 = _load("emk1", "emk1_determinant_seam_ladder.py")


# ---------------- primitive algebra ----------------

def test_primitive_relations():
    assert emk1.mmul(emk1.K2, emk1.K2) == emk1.I2
    assert emk1.mmul(emk1.R2, emk1.R2) == emk1.mscale(emk1.I2, -1)
    assert emk1.mmul(emk1.RK2, emk1.RK2) == emk1.I2
    assert emk1.mmul(emk1.R2, emk1.K2) == \
        emk1.mscale(emk1.mmul(emk1.K2, emk1.R2), -1)


def test_ugd_bracket_R_K_is_2RK():
    assert emk1.commutator(emk1.R2, emk1.K2) == emk1.mscale(emk1.RK2, 2)


def test_quarter_turn_squares_to_minus_identity():
    # the same quarter turn as LAM-2 (iota^2=-1) and CFE-Q (J^2=-I)
    assert emk1.mmul(emk1.R2, emk1.R2) == emk1.mscale(emk1.I2, -1)


# ---------------- the determinant identity ----------------

def test_determinant_identity_on_grid():
    g = [Fr(k) for k in range(-3, 4)]
    for a in g:
        for b in g:
            for c in g:
                for d in g:
                    M = emk1.emk_block(a, b, c, d)
                    assert emk1.det_exact(M) == \
                        emk1.delta_parallel(a, b) + emk1.delta_perp(c, d)


def test_channels_are_independent():
    # seam channel zero, rotational nonzero
    assert emk1.delta_parallel(Fr(2), Fr(2)) == 0
    assert emk1.delta_perp(Fr(3), Fr(1)) != 0
    # rotational zero, seam nonzero
    assert emk1.delta_perp(Fr(2), Fr(2)) == 0
    assert emk1.delta_parallel(Fr(3), Fr(1)) != 0


def test_both_channels_zero_but_operator_nonzero():
    M = emk1.emk_block(Fr(2), Fr(2), Fr(3), Fr(3))
    assert emk1.det_exact(M) == 0
    assert not emk1.is_zero(M)


def test_assembled_block_matches_appendix_form():
    a, b, c, d = Fr(5), Fr(2), Fr(3), Fr(1)
    M = emk1.emk_block(a, b, c, d)
    assert M == [[a - d, b - c], [b + c, a + d]]


# ---------------- projectors and mixing ----------------

def test_projector_identities():
    Pp, Pm = emk1.P_PLUS, emk1.P_MINUS
    assert emk1.mmul(Pp, Pp) == Pp
    assert emk1.mmul(Pm, Pm) == Pm
    assert emk1.is_zero(emk1.mmul(Pp, Pm))
    assert emk1.madd(Pp, Pm) == emk1.I2


def test_commutator_detects_seam_mixing():
    g = [Fr(k) for k in range(-2, 3)]
    for a in g:
        for b in g:
            for c in g:
                for d in g:
                    A = emk1.emk_block(a, b, c, d)
                    comm_zero = emk1.is_zero(emk1.commutator(A, emk1.K2))
                    mix_zero = (
                        emk1.is_zero(emk1.mmul(
                            emk1.mmul(emk1.P_PLUS, A), emk1.P_MINUS))
                        and emk1.is_zero(emk1.mmul(
                            emk1.mmul(emk1.P_MINUS, A), emk1.P_PLUS)))
                    assert comm_zero == mix_zero


# ---------------- Schur / second order ----------------

def test_schur_identity_exact():
    App = emk1.mat([[3, 1], [0, 2]])
    Amm = emk1.mat([[5, 0], [1, 4]])
    Epm = emk1.mat([[Fr(1, 10), 0], [0, Fr(1, 5)]])
    Emp = emk1.mat([[Fr(1, 5), 0], [0, Fr(1, 10)]])
    A4 = [App[0] + Epm[0], App[1] + Epm[1],
          Emp[0] + Amm[0], Emp[1] + Amm[1]]
    schur = emk1.msub(Amm, emk1.mmul(
        emk1.mmul(Emp, emk1.inv_exact(App)), Epm))
    assert emk1.det_exact(A4) == emk1.det_exact(App) * emk1.det_exact(schur)


def test_first_order_closure_is_incomplete():
    """Both diagonal blocks have determinant 1, yet the quadratic
    determinant residue is nonzero — the ladder's key warning."""
    Bpp = emk1.eye(2)
    Bmm = emk1.eye(2)
    Fpm = emk1.mat([[Fr(1, 4), 0], [0, Fr(1, 3)]])
    Fmp = emk1.mat([[Fr(1, 3), 0], [0, Fr(1, 4)]])
    assert emk1.det_exact(Bpp) == 1 and emk1.det_exact(Bmm) == 1
    quad = emk1.trace(emk1.mmul(
        emk1.mmul(emk1.inv_exact(Bpp), Fpm),
        emk1.mmul(emk1.inv_exact(Bmm), Fmp)))
    assert quad == Fr(1, 6)
    assert quad != 0


# ---------------- cut curvature ----------------

def test_cut_commutator_curvature_identity():
    Dx = emk1.mat([[0, 1], [0, 0]])
    Di = emk1.mat([[0, 0], [1, 0]])
    dx, di = Fr(1, 3), Fr(1, 7)
    Cx = emk1.madd(emk1.I2, emk1.mscale(Dx, dx))
    Ci = emk1.madd(emk1.I2, emk1.mscale(Di, di))
    lhs = emk1.msub(emk1.mmul(Ci, Cx), emk1.mmul(Cx, Ci))
    rhs = emk1.mscale(emk1.commutator(Dx, Di), -dx * di)
    assert lhs == rhs


def test_commuting_updates_close_the_path():
    Dx = emk1.mat([[0, 1], [0, 0]])
    Dc = emk1.mat([[2, 0], [0, 2]])
    assert emk1.is_zero(emk1.commutator(Dx, Dc))


# ---------------- winding ----------------

def test_winding_is_exact_integer_and_homotopy_invariant():
    assert emk1.winding_number(emk1.rational_circle(Fr(1))) == 1
    assert emk1.winding_number(emk1.rational_circle(Fr(2))) == 1
    assert emk1.winding_number(emk1.rational_circle(Fr(1, 3))) == 1


def test_loop_excluding_origin_has_zero_winding():
    assert emk1.winding_number(
        emk1.rational_circle(Fr(1), cx=Fr(3))) == 0


def test_winding_uses_no_floats():
    for (x, y) in emk1.rational_circle(Fr(1)):
        assert isinstance(x, Fr) and isinstance(y, Fr)


# ---------------- bridge to CFE ----------------

def test_memoryless_iff_rotation_free():
    g = [Fr(k) for k in range(-2, 3)]
    for a in g:
        for b in g:
            for c in g:
                for d in g:
                    M = emk1.emk_block(a, b, c, d)
                    memoryless = emk1.is_zero(emk1.commutator(M, emk1.K2))
                    assert memoryless == (c == 0 and d == 0)
                    if c == 0 and d == 0:
                        assert emk1.det_exact(M) == emk1.delta_parallel(a, b)


# ---------------- certificate integrity ----------------

def test_certificate_pin_matches_regeneration():
    cert = emk1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMK1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_infinite_dimensional_extensions_not_claimed():
    cert = json.load(open(os.path.join(CERT_DIR, "EMK1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["infinite_dimensional_extensions"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["YM_continuum_gates"] == "not touched"
    assert cb["quantum_gravity"] == "not touched"
    assert cb["CFE_uniqueness_U"].startswith("OPEN")


def test_provenance_records_no_prior_executable_version():
    cert = json.load(open(os.path.join(CERT_DIR, "EMK1_RESULT.json")))
    assert "NONE" in cert["provenance"]["prior_executable_version"]
