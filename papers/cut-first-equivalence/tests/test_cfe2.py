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


cfe2 = _load("cfe2", "cfe2_surjectivity.py")


def test_potential_is_convex_on_lattice():
    # Hessian determinant strictly positive and T>0 everywhere used
    for s in range(1, 5):
        for v in range(1, 5):
            S, V = Fr(s), Fr(v)
            det = (cfe2.U_SS(S, V) * cfe2.U_VV(S, V)
                   - cfe2.U_SV(S, V) ** 2)
            assert det > 0
            assert cfe2.U_S(S, V) > 0


def test_exact_rank_helper():
    assert cfe2.mat_rank([[1, 0], [0, 1]]) == 2
    assert cfe2.mat_rank([[1, 2], [2, 4]]) == 1
    assert cfe2.mat_rank([[Fr(1, 3), Fr(2, 3)], [Fr(1), Fr(2)]]) == 1


def test_all_maxwell_relations_hold_on_potential():
    for s in range(1, 5):
        for v in range(1, 5):
            for res in cfe2.maxwell_residuals(Fr(s), Fr(v)):
                assert res == 0


def test_nonpotential_field_breaks_maxwell():
    assert cfe2.maxwell_from_nonpotential(Fr(2), Fr(3)) != 0


def test_ratio_law_and_mayer_exact():
    for s in range(1, 5):
        for v in range(1, 5):
            S, V = Fr(s), Fr(v)
            T = cfe2.U_S(S, V)
            Cv, Cp, kT, kS, det = cfe2.response_coeffs(S, V, T)
            assert kT / kS == Cp / Cv                     # gamma ratio law
            assert Cp - Cv == T * cfe2.U_SV(S, V) ** 2 \
                / (cfe2.U_SS(S, V) * det)                 # Mayer shape


def test_surjectivity_zero_cokernel():
    cert = json.load(open(os.path.join(CERT_DIR, "CFE2_RESULT.json")))
    t3 = cert["T3_surjectivity_rank_equality"]
    assert t3["surjective"] is True
    assert t3["cokernel_dimension"] == 0
    assert t3["rank_cut_generated_space"] == t3["rank_classical_response_space"]
    assert t3["rank_union"] == t3["rank_classical_response_space"]


def test_cut_basis_is_finite_difference_not_analytic():
    # the cut rank must be reached by centered finite differences of U,
    # a procedure independent of the analytic partials
    h = Fr(1)
    LAT = [(Fr(s), Fr(v)) for s in range(1, 5) for v in range(1, 5)]

    def fd_S(fn, S, V):
        return (fn(S + h, V) - fn(S - h, V)) / (2 * h)

    def fd_V(fn, S, V):
        return (fn(S, V + h) - fn(S, V - h)) / (2 * h)

    def vv(fn):
        return [fn(S, V) for (S, V) in LAT[:8]]

    rows = [
        vv(lambda S, V: fd_S(cfe2.U, S, V)),
        vv(lambda S, V: fd_V(cfe2.U, S, V)),
        vv(lambda S, V: fd_S(lambda s, v: fd_S(cfe2.U, s, v), S, V)),
        vv(lambda S, V: fd_V(lambda s, v: fd_V(cfe2.U, s, v), S, V)),
        vv(lambda S, V: fd_V(lambda s, v: fd_S(cfe2.U, s, v), S, V)),
    ]
    assert cfe2.mat_rank(rows) == 5


def test_renaming_test_random_potentials_also_rank5():
    # structural, not point-luck: same monomial shape, random rational
    # coefficients -> still full rank 5 by cut finite differences
    import random
    random.seed(11)
    h = Fr(1)
    LAT = [(Fr(s), Fr(v)) for s in range(1, 5) for v in range(1, 5)]

    def fd_S(fn, S, V):
        return (fn(S + h, V) - fn(S - h, V)) / (2 * h)

    def fd_V(fn, S, V):
        return (fn(S, V + h) - fn(S, V - h)) / (2 * h)

    for _ in range(4):
        a, b, c, d, e = [Fr(random.randint(6, 12)) for _ in range(5)]

        def Ur(S, V, a=a, b=b, c=c, d=d, e=e):
            return (a * S * S / 2 + b * S * V + c * V * V / 2
                    + d * S * S * V / 2 + e * S * V * V / 2)

        def vv(fn):
            return [fn(S, V) for (S, V) in LAT[:8]]

        rows = [
            vv(lambda S, V: fd_S(Ur, S, V)),
            vv(lambda S, V: fd_V(Ur, S, V)),
            vv(lambda S, V: fd_S(lambda s, v: fd_S(Ur, s, v), S, V)),
            vv(lambda S, V: fd_V(lambda s, v: fd_V(Ur, s, v), S, V)),
            vv(lambda S, V: fd_V(lambda s, v: fd_S(Ur, s, v), S, V)),
        ]
        assert cfe2.mat_rank(rows) == 5


def test_control_separates_recorded():
    cert = json.load(open(os.path.join(CERT_DIR, "CFE2_RESULT.json")))
    assert cert["T3_surjectivity_rank_equality"]["nonclosed_control_separates"] is True


def test_certificate_pin_matches_regeneration():
    cert = cfe2.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_CFE2.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_uniqueness_remains_open():
    cert = json.load(open(os.path.join(CERT_DIR, "CFE2_RESULT.json")))
    ro = cert["remaining_open"]
    assert ro["U_uniqueness"].startswith("OPEN")
    assert ro["generality_beyond_witness_EOS"].startswith("OPEN")
    cb = cert["claim_boundary"]
    assert cb["uniqueness_U"] == "OPEN"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["YM_continuum_gates"] == "not touched"
