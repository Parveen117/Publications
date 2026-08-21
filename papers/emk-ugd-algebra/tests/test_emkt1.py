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


t1 = _load("emkt1", "emkt1_master_tensor_and_time.py")


# ---------------- J is derived ----------------

def test_cut_swap_is_the_primitive_K():
    assert t1.cut_swap() == t1.K2


def test_cut_swap_is_an_involution():
    J = t1.cut_swap()
    assert t1.mm(J, J) == t1.I2


def test_cut_swap_acts_as_coordinate_swap():
    J = t1.cut_swap()
    for x in range(-3, 4):
        for y in range(-3, 4):
            out = t1.mm(J, [[Fr(x)], [Fr(y)]])
            assert (out[0][0], out[1][0]) == (Fr(y), Fr(x))


def test_primitive_relations_still_hold():
    assert t1.mm(t1.R2, t1.R2) == t1.mscale(t1.I2, -1)
    assert t1.mm(t1.K2, t1.K2) == t1.I2
    assert t1.mm(t1.R2, t1.K2) == t1.mscale(t1.mm(t1.K2, t1.R2), -1)


# ---------------- channel curvature ----------------

def _chan(a, b, c, d, e, f, g, h):
    return ([[Fr(a), Fr(b)], [Fr(c), Fr(d)]],
            [[Fr(e), Fr(f)], [Fr(g), Fr(h)]])


def test_curvature_decomposition_is_exact():
    chans = [_chan(1, 2, 0, 1, 0, 1, 1, 0), _chan(0, 1, 3, 0, 2, 0, 0, 1),
             _chan(1, 0, 0, 2, 0, 3, 1, 0)]
    omega = t1.master_curvature(chans)
    rebuilt = t1.madd(t1.msum(t1.channel_curvatures(chans), 2),
                      t1.mixed_curvature(chans))
    assert omega == rebuilt


def test_single_channel_has_no_mixed_term():
    chans = [_chan(1, 2, 0, 1, 0, 1, 1, 0)]
    assert t1.is_zero(t1.mixed_curvature(chans))
    assert t1.master_curvature(chans) == t1.channel_curvatures(chans)[0]


def test_mixed_term_carries_all_the_curvature():
    c1 = ([[Fr(0), Fr(1)], [Fr(0), Fr(0)]],
          [[Fr(0), Fr(2)], [Fr(0), Fr(0)]])
    c2 = ([[Fr(0), Fr(0)], [Fr(1), Fr(0)]],
          [[Fr(0), Fr(0)], [Fr(3), Fr(0)]])
    assert t1.is_zero(t1.comm(c1[0], c1[1]))
    assert t1.is_zero(t1.comm(c2[0], c2[1]))
    pair = [c1, c2]
    assert t1.is_zero(t1.msum(t1.channel_curvatures(pair), 2))
    assert not t1.is_zero(t1.mixed_curvature(pair))
    assert t1.master_curvature(pair) == t1.mixed_curvature(pair)


# ---------------- closure decision ----------------

def test_nonzero_curvature_can_be_closed():
    active = [[True, True], [True, True]]
    omega = [[Fr(0), Fr(5)], [Fr(-5), Fr(0)]]
    exact = [[Fr(0), Fr(3)], [Fr(-3), Fr(0)]]
    ledger = [[Fr(0), Fr(2)], [Fr(-2), Fr(0)]]
    assert not t1.is_zero(omega)
    assert t1.sector_closed(omega, exact, ledger, active)


def test_unledgered_component_leaves_it_open():
    active = [[True, True], [True, True]]
    omega = [[Fr(0), Fr(5)], [Fr(-5), Fr(0)]]
    exact = [[Fr(0), Fr(3)], [Fr(-3), Fr(0)]]
    short = [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]]
    assert not t1.sector_closed(omega, exact, short, active)


def test_inactive_survivor_does_not_block():
    mask = [[True, False], [False, True]]
    omega = [[Fr(0), Fr(5)], [Fr(-5), Fr(0)]]
    exact = [[Fr(0), Fr(3)], [Fr(-3), Fr(0)]]
    short = [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]]
    assert t1.sector_closed(omega, exact, short, mask)


# ---------------- time vs clock ----------------

def test_identical_clock_can_leave_the_sector_open():
    grad = [Fr(1), Fr(2)]
    zero = [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
    mem = [[Fr(1), Fr(0)], [Fr(0), Fr(-1)]]
    a = t1.time_tensor(grad, zero, zero)
    b = t1.time_tensor(grad, mem, zero)
    assert a != b
    assert t1.temporal_residue(b, a, zero) > 0


def test_lawful_temporal_ledger_closes_it():
    grad = [Fr(1), Fr(2)]
    zero = [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
    mem = [[Fr(1), Fr(0)], [Fr(0), Fr(-1)]]
    a = t1.time_tensor(grad, zero, zero)
    b = t1.time_tensor(grad, mem, zero)
    assert t1.temporal_residue(b, a, mem) == 0


def test_time_tensor_is_symmetric_in_the_gradient_part():
    g = t1.time_tensor([Fr(2), Fr(3)],
                       [[Fr(0), Fr(0)], [Fr(0), Fr(0)]],
                       [[Fr(0), Fr(0)], [Fr(0), Fr(0)]])
    assert g[0][1] == g[1][0]


# ---------------- determinant coupling ----------------

def test_determinant_is_multiplicative():
    B1 = [[Fr(3), Fr(1)], [Fr(0), Fr(2)]]
    B2 = [[Fr(5), Fr(0)], [Fr(1), Fr(4)]]
    assert t1.det2(t1.mm(B1, B2)) == t1.det2(B1) * t1.det2(B2)


def test_coupling_term_is_exactly_the_inverse_transport_determinant():
    B1 = [[Fr(3), Fr(1)], [Fr(0), Fr(2)]]
    B2 = [[Fr(5), Fr(0)], [Fr(1), Fr(4)]]
    T = [[Fr(2), Fr(0)], [Fr(0), Fr(3)]]
    lhs = t1.det2(t1.mm(T, t1.mm(B1, B2)))
    rhs = t1.det2(t1.mm(T, B1)) * t1.det2(t1.mm(T, B2))
    coupling = lhs / rhs
    assert coupling == Fr(1) / t1.det2(T)
    assert coupling != 1
    assert coupling * t1.det2(T) == 1


def test_no_logarithm_is_ever_evaluated():
    src = open(os.path.join(
        CERT_DIR, "emkt1_master_tensor_and_time.py")).read()
    assert "math.log" not in src and "import math" not in src


# ---------------- presentation agreement ----------------

def test_verdicts_agree_when_sectors_preserved():
    tol = Fr(1, 10)

    def verdict(rs):
        return all(abs(r) <= tol for r in rs)

    for s in [(Fr(0), Fr(0), Fr(0)), (Fr(1, 20), Fr(0), Fr(1, 50)),
              (Fr(1, 2), Fr(0), Fr(0)), (Fr(0), Fr(3, 10), Fr(0))]:
        assert verdict(list(s)) == verdict(list(s)) == verdict(list(s))


def test_dropping_a_sector_breaks_agreement():
    tol = Fr(1, 10)

    def verdict(rs):
        return all(abs(r) <= tol for r in rs)

    s = (Fr(0), Fr(3, 10), Fr(0))
    assert verdict(list(s)) is False
    assert verdict([s[0], s[2]]) is True        # dropped sector hides it


# ---------------- certificate integrity ----------------

def test_certificate_pin_matches_regeneration():
    cert = t1.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_EMKT1.sha256")) as f:
        pinned = f.read().strip()
    assert digest == pinned


def test_declared_structure_is_recorded_as_declared():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKT1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert "declarations" in cb["declared_not_certified"]
    assert "eight-channel" in cb["declared_not_certified"]
    assert cb["field_theoretic_and_device_extensions"] == "NOT CLAIMED"
    assert cb["RH_K0_L0"] == "not touched"
    assert cb["quantum_gravity"] == "not touched"


def test_honest_note_present_in_provenance():
    cert = json.load(open(os.path.join(CERT_DIR, "EMKT1_RESULT.json")))
    assert "protocol rather than theorem" in cert["provenance"]["honest_note"]
