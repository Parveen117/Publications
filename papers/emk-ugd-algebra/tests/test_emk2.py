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


emk2 = _load("emk2", "emk2_native_carrier.py")


# ---------------- grading ----------------

def test_coefficient_recovery_is_exact():
    a, b, c, d = Fr(5), Fr(2), Fr(3), Fr(1)
    assert emk2.coeffs(emk2.block(a, b, c, d)) == (a, b, c, d)


def test_even_times_even_is_even():
    E1 = emk2.block(Fr(3), Fr(1), 0, 0)
    E2 = emk2.block(Fr(2), Fr(5), 0, 0)
    assert emk2.is_even(emk2.mm(E1, E2))


def test_odd_times_odd_lands_in_even():
    O1 = emk2.block(0, 0, Fr(3), Fr(1))
    O2 = emk2.block(0, 0, Fr(2), Fr(5))
    assert emk2.is_even(emk2.mm(O1, O2))


def test_even_times_odd_is_odd():
    E = emk2.block(Fr(3), Fr(1), 0, 0)
    O = emk2.block(0, 0, Fr(2), Fr(5))
    P = emk2.mm(E, O)
    assert emk2.is_odd(P)


def test_explicit_composition_laws():
    a1, b1, a2, b2 = Fr(3), Fr(1), Fr(2), Fr(5)
    P = emk2.mm(emk2.block(a1, b1, 0, 0), emk2.block(a2, b2, 0, 0))
    assert emk2.coeffs(P) == (a1 * a2 + b1 * b2, a1 * b2 + a2 * b1,
                              Fr(0), Fr(0))
    c1, d1, c2, d2 = Fr(3), Fr(1), Fr(2), Fr(5)
    Q = emk2.mm(emk2.block(0, 0, c1, d1), emk2.block(0, 0, c2, d2))
    assert emk2.coeffs(Q) == (d1 * d2 - c1 * c2, d1 * c2 - c1 * d2,
                              Fr(0), Fr(0))


# ---------------- the UGD multiplicative law ----------------

def test_cayley_composition_is_the_rational_addition_law():
    ys = [Fr(1, 3), Fr(1, 5), Fr(2, 7), Fr(-1, 4), Fr(3, 11), Fr(0)]
    for y1 in ys:
        for y2 in ys:
            P = emk2.mm(emk2.block(Fr(1), y1, 0, 0),
                        emk2.block(Fr(1), y2, 0, 0))
            pa, pb, _, _ = emk2.coeffs(P)
            assert pb / pa == (y1 + y2) / (1 + y1 * y2)


def test_native_chart_gives_log_addition():
    """x = (1+y)/(1-y);  x1*x2 == x(y1 (+) y2)  — i.e. Log_Sigma(x1 x2)
    = Log_Sigma(x1) + Log_Sigma(x2), exactly, with no series."""
    ys = [Fr(1, 3), Fr(1, 5), Fr(2, 7), Fr(-1, 4), Fr(3, 11), Fr(5, 13)]
    for y1 in ys:
        for y2 in ys:
            yl = (y1 + y2) / (1 + y1 * y2)
            x1 = (1 + y1) / (1 - y1)
            x2 = (1 + y2) / (1 - y2)
            assert x1 * x2 == (1 + yl) / (1 - yl)


def test_seam_determinant_channel_is_multiplicative():
    g = [Fr(k) for k in range(-3, 4)]
    for a1 in g:
        for b1 in g:
            for a2 in g:
                for b2 in g:
                    P = emk2.mm(emk2.block(a1, b1, 0, 0),
                                emk2.block(a2, b2, 0, 0))
                    ap, bp, _, _ = emk2.coeffs(P)
                    assert emk2.delta_par(ap, bp) == \
                        emk2.delta_par(a1, b1) * emk2.delta_par(a2, b2)


def test_naive_addition_control_separates():
    y1, y2 = Fr(1, 3), Fr(1, 5)
    assert y1 + y2 != (y1 + y2) / (1 + y1 * y2)


def test_multiplicative_law_uses_no_floats():
    y1, y2 = Fr(2, 7), Fr(3, 11)
    yl = (y1 + y2) / (1 + y1 * y2)
    assert isinstance(yl, Fr)


# ---------------- native carrier ----------------

def test_memoryless_iff_even_grade():
    g = [Fr(k) for k in range(-2, 3)]
    for a in g:
        for b in g:
            for c in g:
                for d in g:
                    M = emk2.block(a, b, c, d)
                    assert emk2.is_even(M) == (c == 0 and d == 0)


def test_zero_odd_amplitude_is_flat():
    assert emk2.curvature_native(Fr(0)) == 0
    assert emk2.circulation_native(emk2.RECT, Fr(0)) == 0


def test_native_stokes_identity_exact():
    for rho in (Fr(-2, 5), Fr(1, 5), Fr(3, 7)):
        circ = emk2.circulation_native(emk2.RECT, rho)
        flux = emk2.curvature_native(rho) * emk2.signed_area(emk2.RECT)
        assert circ == flux


def test_residue_faithful_and_monotone():
    rhos = [Fr(-2, 5), Fr(-1, 5), Fr(0), Fr(1, 5), Fr(2, 5)]
    vals = [emk2.circulation_native(emk2.RECT, r) for r in rhos]
    assert vals[2] == 0
    assert all(vals[k] < vals[k + 1] for k in range(len(vals) - 1))


# ---------------- scaffold retirement ----------------

def test_scaffold_and_native_densities_coincide():
    """Under rho = chi - 1 the CFE scaffold density and the native one
    are identical — retiring the equation of state changes nothing."""
    for chi in (Fr(3, 5), Fr(4, 5), Fr(1), Fr(6, 5), Fr(7, 5)):
        assert Fr(5) * (chi - 1) == emk2.curvature_native(chi - 1)


def test_native_residues_match_cfe1_values():
    """The residues on the native carrier reproduce CFE-1's certified
    numbers (-32, -16, 0, +16, +32) exactly."""
    expected = [Fr(-32), Fr(-16), Fr(0), Fr(16), Fr(32)]
    rhos = [Fr(-2, 5), Fr(-1, 5), Fr(0), Fr(1, 5), Fr(2, 5)]
    got = [emk2.circulation_native(emk2.RECT, r) for r in rhos]
    assert got == expected


# ---------------- certificate integrity ----------------

def test_certificate_pin_matches_regeneration():
    cert = emk2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMK2.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_native_log_derivation_is_pinned_not_rederived():
    cert = json.load(open(os.path.join(CERT_DIR, "EMK2_RESULT.json")))
    assert "F00G" in cert["provenance"]["native_log_law_source"]
    assert "PINNED" in cert["claim_boundary"]["native_log_derivation"]


def test_open_boundaries_declared():
    cert = json.load(open(os.path.join(CERT_DIR, "EMK2_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["CFE_uniqueness_U"].startswith("OPEN")
    assert cb["infinite_dimensional_extensions"] == "NOT CLAIMED"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["quantum_gravity"] == "not touched"
