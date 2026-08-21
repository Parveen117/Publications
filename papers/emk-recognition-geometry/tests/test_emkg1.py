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


g1 = _load("emkg1", "emkg1_rotational_seam_metric.py")

KAPPAS = (Fr(3), Fr(1, 5), Fr(0), Fr(-2), Fr(7, 4))
VS = (Fr(0), Fr(1, 2), Fr(2), Fr(-3, 4), Fr(1, 3))


# ---------------- involution ----------------

def test_involution_swaps_coordinates():
    for x in range(-4, 5):
        for y in range(-4, 5):
            assert g1.J_real(x, y) == (y, x)


def test_J_is_an_involution():
    for x in range(-3, 4):
        for y in range(-3, 4):
            assert g1.J_real(*g1.J_real(x, y)) == (x, y)


def test_seam_coordinate_invariant_normal_negated():
    for x in range(-4, 5):
        for y in range(-4, 5):
            u0, v0 = g1.seam_coords(x, y)
            u1, v1 = g1.seam_coords(*g1.J_real(x, y))
            assert u1 == u0 and v1 == -v0


def test_fixed_set_is_the_45_degree_line():
    fixed = [(x, y) for x in range(-4, 5) for y in range(-4, 5)
             if g1.J_real(x, y) == (x, y)]
    assert len(fixed) == 9
    assert all(x == y for (x, y) in fixed)
    assert all(g1.seam_coords(x, y)[1] == 0 for (x, y) in fixed)


# ---------------- metric ----------------

def test_determinant_is_W_and_positivity():
    for k in KAPPAS:
        for v in VS:
            W = g1.W_quadratic(v, k)
            if W > 0:
                assert W > 0            # eigenvalues W and 1
            else:
                assert W <= 0


def test_degeneracy_locus_is_exact():
    k = Fr(-4)
    assert g1.W_quadratic(Fr(1, 2), k) == 0
    assert g1.W_quadratic(Fr(-1, 2), k) == 0
    assert g1.W_quadratic(Fr(1, 4), k) != 0


# ---------------- connection and curvature ----------------

def test_christoffels_match_levi_civita():
    for k in KAPPAS:
        for v in VS:
            W = g1.W_quadratic(v, k)
            if W == 0:
                continue
            lc_u_uv = Fr(1, 2) * (Fr(1) / W) * g1.dW_quadratic(v, k)
            lc_v_uu = -Fr(1, 2) * g1.dW_quadratic(v, k)
            assert g1.christoffel_u_uv(v, k) == lc_u_uv
            assert g1.christoffel_v_uu(v, k) == lc_v_uu


def test_curvature_two_routes_agree():
    for k in KAPPAS:
        for v in VS:
            W = g1.W_quadratic(v, k)
            if W == 0:
                continue
            assert g1.riemann_uvuv(v, k) / W == \
                g1.gaussian_curvature_closed(v, k)
            assert g1.gaussian_curvature(v, k) == \
                g1.gaussian_curvature_closed(v, k)


def test_seam_curvature_is_minus_kappa():
    for k in KAPPAS:
        assert g1.gaussian_curvature(Fr(0), k) == -k


def test_curvature_is_rational_no_floats():
    K = g1.gaussian_curvature(Fr(1, 3), Fr(7, 4))
    assert isinstance(K, Fr)


# ---------------- geodesic seam ----------------

def test_reflection_symmetric_seam_is_geodesic():
    for k in KAPPAS:
        assert g1.dW_quadratic(Fr(0), k) == 0
        assert g1.christoffel_v_uu(Fr(0), k) == 0


def test_odd_perturbation_breaks_the_geodesic_seam():
    mu, k = Fr(3, 5), Fr(2)
    dW0 = g1.dW_odd(Fr(0), k, mu)
    assert dW0 == mu and mu != 0
    assert -dW0 / 2 != 0


def test_even_warp_has_vanishing_first_derivative_at_zero():
    for k in KAPPAS:
        assert g1.W_quadratic(Fr(1, 3), k) == g1.W_quadratic(Fr(-1, 3), k)


# ---------------- sign trichotomy ----------------

def test_positive_kappa_gives_negative_curvature():
    for v in (Fr(0), Fr(1, 2), Fr(3), Fr(-5, 2)):
        assert g1.W_quadratic(v, Fr(3)) > 0
        assert g1.gaussian_curvature(v, Fr(3)) < 0


def test_zero_kappa_is_euclidean():
    for v in (Fr(0), Fr(1, 2), Fr(3), Fr(-5, 2)):
        assert g1.gaussian_curvature(v, Fr(0)) == 0


def test_negative_kappa_positive_curvature_inside_strip():
    k = Fr(-4)
    for v in (Fr(0), Fr(1, 5), Fr(-1, 4)):
        assert v * v < Fr(-1) / k
        assert g1.gaussian_curvature(v, k) > 0
    assert g1.W_quadratic(Fr(1, 2), k) == 0        # strip boundary


# ---------------- the separation ----------------

def test_metric_and_recognition_curvature_are_independent():
    # Euclidean metric, nonzero recognition curvature
    assert g1.gaussian_curvature(Fr(0), Fr(0)) == 0
    assert Fr(5) * Fr(2, 5) != 0
    # curved metric, zero recognition curvature
    assert g1.gaussian_curvature(Fr(0), Fr(3)) != 0
    assert Fr(5) * Fr(0) == 0


# ---------------- certificate integrity ----------------

def test_certificate_pin_matches_regeneration():
    cert = g1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKG1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_metric_family_declared_and_identification_not_claimed():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert "DECLARED" in cb["metric_family_is_declared"]
    assert cb["identification_with_recognition_two_form"].startswith(
        "NOT CLAIMED")
    assert cb["global_and_completeness_statements"] == "NOT CLAIMED"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["quantum_gravity"] == "not touched"


def test_provenance_records_no_prior_executable_version():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG1_RESULT.json")))
    assert "NONE" in cert["provenance"]["prior_executable_version"]
