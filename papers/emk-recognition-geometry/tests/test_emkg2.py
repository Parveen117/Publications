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


g2 = _load("emkg2", "emkg2_global_quotient_holonomy.py")


# ---------------- T1: closed geodesic seam ----------------

def test_geodesic_criterion_equivalent_in_A_and_W():
    for c in g2.CS:
        A = g2.A_rational(c)
        W = g2.pmul(A, A)
        assert g2.peval(A, Fr(0)) > 0
        assert (g2.peval(g2.pderiv(A), Fr(0)) == 0) == \
            (g2.peval(g2.pderiv(W), Fr(0)) == 0)


def test_seam_length_squared_is_exact():
    for kappa in g2.KAPPAS:
        W = g2.W_family(kappa)
        assert g2.L_PERIOD ** 2 * g2.peval(W, Fr(0)) == g2.L_PERIOD ** 2


def test_scaled_warp_changes_the_length():
    W = g2.pmul([Fr(9)], g2.W_family(Fr(1)))
    assert g2.L_PERIOD ** 2 * g2.peval(W, Fr(0)) == 9 * g2.L_PERIOD ** 2


def test_odd_perturbation_breaks_geodesy_not_closure():
    for mu in (Fr(1), Fr(-3, 7), Fr(5, 2)):
        W = g2.padd(g2.W_family(Fr(2)), [Fr(0), mu])
        assert g2.peval(g2.pderiv(W), Fr(0)) == mu != 0
        assert g2.peval(W, Fr(0)) > 0            # still closes


# ---------------- T2: quotient is not doubling ----------------

def test_involution_fixed_set_is_exactly_the_seam():
    pts = [(u, v) for u in g2.GRID_U for v in g2.GRID_V]
    fixed = [p for p in pts if g2.j_involution(p) == p]
    assert all(v == 0 for (_, v) in fixed)
    assert len(fixed) == len(g2.GRID_U)


def test_orbit_sizes_two_off_seam_one_on():
    pts = [(u, v) for u in g2.GRID_U for v in g2.GRID_V]
    orbits = {}
    for p in pts:
        orbits.setdefault((p[0], abs(p[1])), set()).add(p)
    for (u, av), orb in orbits.items():
        assert len(orb) == (1 if av == 0 else 2)


def test_side_target_does_not_factor_through_the_quotient():
    pts = [(u, v) for u in g2.GRID_U for v in g2.GRID_V if v != 0]
    for p in pts:
        q = g2.j_involution(p)
        sp = 1 if p[1] > 0 else -1
        sq = 1 if q[1] > 0 else -1
        assert sp == -sq


def test_double_has_twice_the_points_off_the_seam():
    pts = [(u, v) for u in g2.GRID_U for v in g2.GRID_V]
    for av in {abs(v) for v in g2.GRID_V}:
        double = [p for p in pts if abs(p[1]) == av]
        quotient = {(p[0], abs(p[1])) for p in double}
        assert len(double) == (1 if av == 0 else 2) * len(quotient)


# ---------------- T3: gluing criterion ----------------

def test_even_warps_glue_smoothly_to_themselves():
    for kappa in g2.KAPPAS:
        W = g2.W_family(kappa)
        assert g2.reflect_poly(W) == W
        assert g2.glues_smoothly(W, W)


def test_odd_term_fails_at_first_jet_but_glues_to_its_reflection():
    W = g2.padd(g2.W_family(Fr(2)), [Fr(0), Fr(3)])
    assert not g2.glues_smoothly(W, W)
    assert g2.jet(W, 1) == 3
    assert g2.glues_smoothly(W, g2.reflect_poly(W))


def test_mismatched_pair_fails_at_the_identified_jet():
    Wa, Wb = g2.W_family(Fr(2)), g2.W_family(Fr(5))
    assert not g2.glues_smoothly(Wa, Wb)
    first = next(n for n in range(7)
                 if g2.jet(Wa, n) != (-1) ** n * g2.jet(Wb, n))
    assert first == 2


def test_jets_of_polynomial_warps_are_exact():
    W = g2.W_family(Fr(7, 3))
    assert g2.jet(W, 0) == 1
    assert g2.jet(W, 1) == 0
    assert g2.jet(W, 2) == 2 * Fr(7, 3)
    assert g2.jet(W, 3) == 0


# ---------------- T4: global classification ----------------

def test_positive_kappa_is_bounded_below_by_one():
    for kappa in (Fr(3), Fr(1, 5), Fr(10 ** 3)):
        W = g2.W_family(kappa)
        for v in (Fr(0), Fr(17), Fr(-10 ** 4), Fr(1, 9)):
            assert g2.peval(W, v) >= 1


def test_flat_case_is_identically_one():
    W = g2.W_family(Fr(0))
    for v in (Fr(0), Fr(6), Fr(-5, 3)):
        assert g2.peval(W, v) == 1


def test_negative_kappa_degeneracy_at_exact_finite_distance():
    W = g2.W_family(Fr(-4))
    assert Fr(-1) / Fr(-4) == Fr(1, 4)
    assert g2.peval(W, Fr(1, 2)) == 0
    assert g2.peval(W, Fr(-1, 2)) == 0
    assert g2.peval(W, Fr(1, 4)) > 0


def test_curvature_growth_witness_is_strictly_increasing():
    W = g2.W_family(Fr(-4))
    seq = [Fr(1, 2) - Fr(1, 10 ** k) for k in range(1, 6)]
    mags = [abs(Fr(-4) / (g2.peval(W, v) ** 2)) for v in seq]
    assert all(a < b for a, b in zip(mags, mags[1:]))
    assert mags[-1] > 10 ** 8


# ---------------- T5: holonomy two routes ----------------

def test_two_routes_agree_on_every_rectangle_and_warp():
    rects = ((Fr(0), Fr(2), Fr(-1), Fr(1)),
             (Fr(-1), Fr(3), Fr(0), Fr(1, 2)),
             (Fr(1, 2), Fr(5, 2), Fr(-2), Fr(1)))
    for c in g2.CS:
        A = g2.A_rational(c)
        for r in rects:
            assert g2.holonomy_boundary(A, *r) == g2.holonomy_area(A, *r)


def test_holonomy_matches_the_closed_form():
    A = g2.A_rational(Fr(5, 3))
    u0, u1, v0, v1 = Fr(0), Fr(3), Fr(-1), Fr(2)
    closed = -(u1 - u0) * (g2.peval(g2.pderiv(A), v1)
                           - g2.peval(g2.pderiv(A), v0))
    assert g2.holonomy_boundary(A, u0, u1, v0, v1) == closed


def test_symmetric_strip_formula():
    for c in g2.CS:
        A = g2.A_rational(c)
        b, u0, u1 = Fr(3, 4), Fr(0), Fr(2)
        assert g2.holonomy_boundary(A, u0, u1, -b, b) == \
            -2 * (u1 - u0) * g2.peval(g2.pderiv(A), b)


def test_seam_loop_holonomy_is_exactly_zero():
    for c in g2.CS:
        A = g2.A_rational(c)
        assert g2.holonomy_boundary(
            A, Fr(0), g2.L_PERIOD, Fr(0), Fr(0)) == 0


def test_polynomial_division_certifies_K_times_A():
    A = g2.A_rational(Fr(3))
    App = g2.pderiv(g2.pderiv(A))
    num = g2.pscale(Fr(-1), App)
    quot, rem = g2.pdivmod(g2.pmul(num, A), A)
    assert rem == [Fr(0)] and quot == num


# ---------------- T6: the two holonomies ----------------

def test_levi_civita_trivial_with_nontrivial_monodromy():
    A = g2.A_rational(Fr(3))
    assert g2.holonomy_boundary(
        A, Fr(0), g2.L_PERIOD, Fr(0), Fr(0)) == 0
    H = g2.mmul(g2.NEG_I, g2.meye())
    assert H == g2.NEG_I != g2.meye()


def test_levi_civita_nontrivial_with_trivial_monodromy():
    A = g2.A_rational(Fr(3))
    assert g2.holonomy_boundary(A, Fr(0), Fr(2), Fr(0), Fr(1)) != 0
    N = [[Fr(0), Fr(1)], [Fr(0), Fr(0)]]
    H = g2.mmul(g2.nilpotent_exp(Fr(2), N),
                g2.nilpotent_exp(Fr(-2), N))
    assert H == g2.meye()


def test_monodromy_ledger_closes_exactly():
    assert g2.mmul(g2.NEG_I, g2.NEG_I) == g2.meye()


def test_global_closure_needs_every_active_sector():
    r = {"endpoint": Fr(0), "levi_civita": Fr(0),
         "recognition": Fr(1), "sheet": Fr(2, 3)}
    allact = {s: True for s in g2.SECTORS}
    assert not g2.survives(r, allact)
    visible = {s: (s in ("endpoint", "levi_civita"))
               for s in g2.SECTORS}
    assert g2.survives(r, visible)
    assert g2.survives(r, allact,
                       {"recognition": Fr(1), "sheet": Fr(2, 3)})
    assert not g2.survives(r, allact, {"recognition": Fr(1)})


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = g2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKG2.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_declares_families_and_refuses_identification():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG2_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["declared_families"].startswith("DECLARED")
    assert cb["blowup_is_a_growth_witness"].startswith("NOT CLAIMED")
    assert cb["identification_of_holonomies"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "emkg2_global_quotient_holonomy.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src
    assert "** 0.5" not in src and "**0.5" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKG2_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
