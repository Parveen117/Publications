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


cfeq = _load("cfeq", "cfeq_quantum_witness.py")


def test_psd_helper_exact():
    assert cfeq.psd_2x2(Fr(1), Fr(0), Fr(1))
    assert not cfeq.psd_2x2(Fr(1), Fr(2), Fr(1))
    assert cfeq.psd_2x2(Fr(0), Fr(0), Fr(0))


def test_markovian_limit_is_cp_divisible():
    div, step = cfeq.cp_divisible(Fr(1))
    assert div is True and step is None


def test_every_off_equilibrium_breaks_cp_divisibility():
    for chi in (Fr(3, 5), Fr(4, 5), Fr(6, 5), Fr(7, 5)):
        div, step = cfeq.cp_divisible(chi)
        assert div is False and step is not None


def test_geometric_residue_zero_iff_markovian():
    assert cfeq.circulation_loop(cfeq.DIAMOND, Fr(1)) == 0
    assert cfeq.circulation_loop(cfeq.DIAMOND, Fr(6, 5)) != 0


def test_quantum_stokes_identity_exact():
    for chi in (Fr(3, 5), Fr(7, 5), Fr(2)):
        circ = cfeq.circulation_loop(cfeq.DIAMOND, chi)
        flux = cfeq.curvature_density_q(chi) * cfeq.signed_area(cfeq.DIAMOND)
        assert circ == flux


def test_residue_strictly_monotone_in_memory():
    chis = [Fr(3, 5), Fr(4, 5), Fr(1), Fr(6, 5), Fr(7, 5)]
    vals = [cfeq.circulation_loop(cfeq.DIAMOND, c) for c in chis]
    assert all(vals[k] < vals[k + 1] for k in range(len(vals) - 1))
    assert vals[0] < 0 and vals[-1] > 0


def test_quarter_turn_squares_to_minus_identity():
    J = [[Fr(0), Fr(-1)], [Fr(1), Fr(0)]]

    def mm(P, Q):
        return [[P[0][0] * Q[0][0] + P[0][1] * Q[1][0],
                 P[0][0] * Q[0][1] + P[0][1] * Q[1][1]],
                [P[1][0] * Q[0][0] + P[1][1] * Q[1][0],
                 P[1][0] * Q[0][1] + P[1][1] * Q[1][1]]]
    assert mm(J, J) == [[Fr(-1), Fr(0)], [Fr(0), Fr(-1)]]


def test_signed_area_diamond():
    assert cfeq.signed_area(cfeq.DIAMOND) == Fr(1, 2)


def test_certificate_pin_matches_regeneration():
    cert = cfeq.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_CFEQ.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_quantum_gravity_is_not_claimed():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEQ_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["quantum_gravity"].startswith("NOT TOUCHED")
    assert cb["generality_beyond_qubit_witness"] == "OPEN"
    assert cb["full_open_system_theorem"] == "OPEN"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["YM_continuum_gates"] == "not touched"


def test_sequence_context_recorded():
    cert = json.load(open(os.path.join(CERT_DIR, "CFEQ_RESULT.json")))
    sc = cert["sequence_context"]
    assert "quantum thermodynamics" in sc["step"]
    assert "horizon" in sc["next"].lower()
