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


u = _load("cfeu", "cfeu_uniqueness.py")


# ---------------- the complex ----------------

def test_complex_ranks_and_composition():
    assert u.rank(u.D0) == 5
    assert u.rank(u.D1) == 1
    for j in range(6):
        col = [u.D0[i][j] for i in range(6)]
        assert u.mv(u.D1, col) == [Fr(0)]


def test_kernel_of_d0_is_constants():
    ker = u.nullspace(u.D0)
    assert len(ker) == 1
    assert ker[0][0] != 0 and all(x == 0 for x in ker[0][1:])


def test_closed_iff_gradient_constructively():
    for k in u.nullspace(u.D1):
        phi = u.potential_of_closed(k)
        assert u.mv(u.D0, phi) == k
    w = [Fr(1), Fr(2), Fr(-3, 7), Fr(4), Fr(-3, 7), Fr(0)]
    assert u.mv(u.D1, w) == [Fr(0)]
    assert u.mv(u.D0, u.potential_of_closed(w)) == w


def test_nonclosed_form_has_no_potential():
    w = [Fr(0), Fr(0), Fr(0), Fr(0), Fr(1), Fr(0)]     # p dv
    assert u.mv(u.D1, w) != [Fr(0)]
    ker1 = u.nullspace(u.D1)
    # w is independent of the closed space
    assert u.rank(ker1 + [w]) == 6


# ---------------- the uniqueness theorem ----------------

def test_annihilator_of_memoryless_space_is_one_dimensional():
    ker1 = u.nullspace(u.D1)
    ann = u.nullspace([k for k in ker1])
    assert len(ann) == 1


def test_residue_functional_spans_the_annihilator():
    loop = u.rect_loop(u.P_LO, u.P_HI, u.V_LO, u.V_HI)
    R = []
    for j in range(6):
        e = [Fr(0)] * 6
        e[j] = Fr(1)
        R.append(u.circulation_of(e, loop))
    for k in u.nullspace(u.D1):
        assert u.dot(R, k) == 0
    assert any(x != 0 for x in R)
    ann = u.nullspace([k for k in u.nullspace(u.D1)])[0]
    ratio = next(rj / gj for rj, gj in zip(R, ann) if gj != 0)
    assert all(rj == ratio * gj for rj, gj in zip(R, ann))


def test_residue_equals_area_times_curvature_map():
    loop = u.rect_loop(u.P_LO, u.P_HI, u.V_LO, u.V_HI)
    for j in range(6):
        e = [Fr(0)] * 6
        e[j] = Fr(1)
        assert u.circulation_of(e, loop) == u.AREA * u.D1[0][j]


# ---------------- the dial and the ladder ----------------

def test_witness_vector_matches_pinned_cfe1_fields():
    for chi in (Fr(1), Fr(7, 5), Fr(-2)):
        a, b, c, e, f, g = u.witness_vec(chi)
        for p in (Fr(0), Fr(2), Fr(-1, 3)):
            for v in (Fr(0), Fr(1), Fr(5, 2)):
                assert u.cfe1.lambda_p(p, v, chi) == a + b * p + c * v
                assert u.cfe1.lambda_v(p, v, chi) == e + f * p + g * v


def test_dial_direction_is_5_p_dv_and_excites_generator():
    d = [x - y for x, y in zip(u.witness_vec(Fr(3)), u.witness_vec(Fr(2)))]
    assert d == [Fr(0), Fr(0), Fr(0), Fr(0), u.cfe1.BETA, Fr(0)]
    assert u.mv(u.D1, d) == [u.cfe1.BETA]


def test_ladder_reproduces_cfe1_exactly_two_routes():
    loop = u.rect_loop(u.P_LO, u.P_HI, u.V_LO, u.V_HI)
    vals = []
    for chi in u.CHIS:
        mine = u.circulation_of(u.witness_vec(chi), loop)
        assert mine == u.cfe1.circulation(loop, chi)
        vals.append(mine)
    assert vals == [Fr(-32), Fr(-16), Fr(0), Fr(16), Fr(32)]


def test_equal_residue_implies_equal_obstruction_class():
    loop = u.rect_loop(u.P_LO, u.P_HI, u.V_LO, u.V_HI)
    base = u.witness_vec(Fr(6, 5))
    grad = u.mv(u.D0, [Fr(0), Fr(-2), Fr(1), Fr(0), Fr(4), Fr(1)])
    other = [x + y for x, y in zip(base, grad)]
    assert u.circulation_of(base, loop) == u.circulation_of(other, loop)
    diff = [x - y for x, y in zip(other, base)]
    assert u.mv(u.D1, diff) == [Fr(0)]
    assert u.mv(u.D0, u.potential_of_closed(diff)) == diff


# ---------------- discrete side ----------------

def test_cell_residue_matrix_rank_one_kernel_matches():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEU_RESULT.json")))
    t5 = cert["T5_discrete_residue_matrix_rank_1"]
    assert t5["residue_matrix_rank"] == 1
    assert t5["cells"] == 16
    assert t5["verdict"] == "PASS"


def test_discrete_stokes_per_cell_for_witness():
    p4, v4 = (u.P_HI - u.P_LO) / 4, (u.V_HI - u.V_LO) / 4
    for chi in (Fr(4, 5), Fr(6, 5)):
        wv = u.witness_vec(chi)
        for i in range(4):
            for j in range(4):
                loop = u.rect_loop(u.P_LO + i * p4, u.P_LO + (i + 1) * p4,
                                   u.V_LO + j * v4, u.V_LO + (j + 1) * v4)
                assert u.circulation_of(wv, loop) == \
                    u.cfe1.BETA * (chi - 1) * p4 * v4


# ---------------- boundary control ----------------

def test_enlarged_algebra_rank_three_declared_rank_one():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEU_RESULT.json")))
    t6 = cert["T6_boundary_control_rank_3"]
    assert t6["rank_enlarged"] == 3
    assert t6["rank_declared_inside"] == 1


# ---------------- integrity, flagship status, boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = u.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_CFEU.sha256")) as f:
        assert digest == f.read().strip()


def test_flagship_all_parts_certified_with_scope():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEU_RESULT.json")))
    fs = cert["flagship_status"]
    for part in ("part_1_memoryless_limit", "part_2_residue",
                 "part_S_surjectivity", "part_U_uniqueness"):
        assert fs[part].startswith("CERTIFIED")
    assert "DECLARED RESPONSE ALGEBRA" in fs["scope"]
    cb = cert["claim_boundary"]
    assert cb["generality"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(CERT_DIR, "cfeu_uniqueness.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEU_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
