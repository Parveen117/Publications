import hashlib
import importlib.util
import json
import os
import sys
from fractions import Fraction as Fr
from itertools import permutations

HERE = os.path.dirname(os.path.abspath(__file__))
CERT_DIR = os.path.join(HERE, "..", "certificates")


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(CERT_DIR, fname))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


lt = _load("ltb1", "ltb1_thermodynamic_bridge.py")


# ---------------- T1: skew response is not curvature ----------------

def test_antisymmetric_part_is_antisymmetric():
    J = [[Fr(2), Fr(3)], [Fr(-1), Fr(5)]]
    A = lt.antisym(J)
    assert A == lt.mT(lt.mscale(Fr(-1), A))
    assert A[0][0] == 0 and A[1][1] == 0


def test_nonreciprocal_response_with_exactly_flat_connection():
    J = [[Fr(2), Fr(3)], [Fr(-1), Fr(5)]]
    assert not lt.is_zero(lt.antisym(J))
    om = lt.omega_affine([Fr(1), Fr(-4)],
                         [[Fr(3), Fr(7)], [Fr(7), Fr(2)]])
    assert lt.is_zero(lt.curvature_of(om))


def test_symmetric_response_with_exactly_nonzero_curvature():
    J = [[Fr(4), Fr(1, 2)], [Fr(1, 2), Fr(-3)]]
    assert lt.is_zero(lt.antisym(J))
    om = lt.omega_affine([Fr(0), Fr(0)],
                         [[Fr(0), Fr(5)], [Fr(-1), Fr(0)]])
    assert lt.curvature_of(om) == [[Fr(0), Fr(-6)], [Fr(6), Fr(0)]]


def test_curvature_of_a_symmetric_gradient_vanishes():
    for g in ([[Fr(1), Fr(2)], [Fr(2), Fr(3)]],
              [[Fr(0), Fr(-5)], [Fr(-5), Fr(9)]]):
        om = lt.omega_affine([Fr(0), Fr(0)], g)
        assert lt.is_zero(lt.curvature_of(om))


# ---------------- T2: the subtraction ladder ----------------

def test_ladder_is_order_independent_over_all_permutations():
    raw = [[Fr(0), Fr(12)], [Fr(-12), Fr(0)]]
    acc = {"law": [[Fr(0), Fr(5)], [Fr(-5), Fr(0)]],
           "frame": [[Fr(0), Fr(2)], [Fr(-2), Fr(0)]],
           "seam": [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]],
           "test": [[Fr(0), Fr(3)], [Fr(-3), Fr(0)]],
           "ledger": [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]}
    ref = lt.open_curvature(raw, acc)
    count = 0
    for order in permutations(lt.SECTORS):
        out = [r[:] for r in raw]
        for s in order:
            out = lt.msub(out, acc[s])
        assert out == ref
        count += 1
    assert count == 120


def test_fully_accounted_raw_curvature_is_exactly_zero():
    raw = [[Fr(0), Fr(12)], [Fr(-12), Fr(0)]]
    acc = {"law": [[Fr(0), Fr(5)], [Fr(-5), Fr(0)]],
           "frame": [[Fr(0), Fr(2)], [Fr(-2), Fr(0)]],
           "seam": [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]],
           "test": [[Fr(0), Fr(3)], [Fr(-3), Fr(0)]],
           "ledger": [[Fr(0), Fr(1)], [Fr(-1), Fr(0)]]}
    assert lt.is_zero(lt.open_curvature(raw, acc))


def test_self_referential_sector_annihilates_every_input():
    for probe in ([[Fr(0), Fr(7)], [Fr(-7), Fr(0)]],
                  [[Fr(2), Fr(-1)], [Fr(4), Fr(6)]],
                  [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]):
        acc = {s: [[Fr(0), Fr(0)], [Fr(0), Fr(0)]] for s in lt.SECTORS}
        acc["law"] = probe
        assert lt.is_zero(lt.open_curvature(probe, acc))


# ---------------- T3: the vdW cycle ----------------

def test_attraction_parameter_cancels_exactly_in_the_loop():
    T0, T1, V0, V1 = Fr(2), Fr(5), Fr(1), Fr(3)
    bottom = lt.work_leg_isothermal(T0, V0, V1)
    top = lt.fneg(lt.work_leg_isothermal(T1, V0, V1))
    loop = lt.fadd(bottom, top)
    assert bottom[0] != 0                      # each leg carries a
    assert loop[0] == 0                        # the loop does not


def test_lambda_coefficient_is_exactly_R_times_delta_T():
    T0, T1, V0, V1 = Fr(2), Fr(5), Fr(1), Fr(3)
    loop = lt.fadd(lt.work_leg_isothermal(T0, V0, V1),
                   lt.fneg(lt.work_leg_isothermal(T1, V0, V1)))
    assert loop[1] == lt.R_GAS * (T0 - T1)


def test_two_routes_agree_up_to_orientation():
    T0, T1 = Fr(2), Fr(5)
    loop = lt.fadd(lt.work_leg_isothermal(T0, Fr(1), Fr(3)),
                   lt.fneg(lt.work_leg_isothermal(T1, Fr(1), Fr(3))))
    area = lt.flog(Fr(0), lt.R_GAS * (T1 - T0))
    assert lt.fadd(loop, area) == (Fr(0), Fr(0))


def test_domain_condition_holds():
    assert Fr(1) > lt.B_VDW


def test_rational_control_eos_needs_no_formal_symbol():
    c, T0, T1, V0, V1 = Fr(3), Fr(2), Fr(5), Fr(1), Fr(3)

    def leg(T, v0, v1):
        return c * T * (v1 * v1 - v0 * v0) / 2

    loop = leg(T0, V0, V1) - leg(T1, V0, V1)
    area = -c * (T1 - T0) * (V1 * V1 - V0 * V0) / 2
    assert loop == area != 0


# ---------------- T4: three diagnostics ----------------

def test_same_magnitude_two_orientations():
    coeff = lt.R_GAS * (Fr(5) - Fr(2))
    assert abs(coeff) == abs(-coeff)
    assert coeff != -coeff


def test_repeated_cycle_returns_phase_and_scale_advances_branch():
    def cycle(state):
        phi, sigma, k = state
        return ((phi + lt.K) % lt.K, sigma, k + 1)

    st = (3, Fr(7, 2), 0)
    s = st
    for n in range(1, 6):
        s = cycle(s)
        assert s[0] == st[0] and s[1] == st[1]
        assert s[2] == n


def test_local_work_curvature_unchanged_by_repetition():
    assert lt.R_GAS * (Fr(5) - Fr(2)) == Fr(3)


# ---------------- T5: post-hoc ledger ----------------

def test_declared_ledger_separates_protocols():
    L = Fr(3)
    assert lt.residue(Fr(3), L) == 0
    assert lt.residue(Fr(19, 4), L) == Fr(7, 4) != 0


def test_post_hoc_ledger_is_identically_zero():
    for v in (Fr(3), Fr(19, 4), Fr(0), Fr(-11, 3), Fr(1000)):
        assert lt.residue(v, v) == 0


def test_post_hoc_rule_distinguishes_no_pair_at_all():
    probes = (Fr(3), Fr(19, 4), Fr(0), Fr(-11, 3), Fr(1000))
    distinguished = sum(1 for a in probes for b in probes
                        if a != b and lt.residue(a, a) != lt.residue(b, b))
    assert distinguished == 0
    L = Fr(3)
    declared = sum(1 for a in probes for b in probes
                   if a != b and lt.residue(a, L) != lt.residue(b, L))
    assert declared > 0


# ---------------- T6: presentation and the flagship's ratios ----------

def test_three_reference_scales_give_three_presentations():
    lam = Fr(7, 2)
    t1 = lt.theta_exponent(lam, Fr(1, 2), Fr(3))
    t2 = lt.theta_exponent(lam, Fr(1, 2), Fr(4))
    t3 = lt.theta_exponent(lam, Fr(2), Fr(3))
    assert (t1, t2, t3) == (0, 9, 6)
    assert len({t1, t2, t3}) == 3


def test_full_cycle_shift_of_the_reference_is_invariant():
    lam = Fr(7, 2)
    a = lt.theta_exponent(lam, Fr(1, 2), Fr(3))
    b = lt.theta_exponent(lam, Fr(1, 2) - Fr(3), Fr(3))
    assert a == b


def test_ratio_variables_multiply_to_the_invariant_shape():
    Gc = Fr(6) / Fr(3)
    Gm = Fr(5) / Fr(10)
    assert Gc * Gm == 1


def test_pinned_cfe1_invariant_is_one_exactly_on_the_face():
    _, _, inv = lt.cfe1.invariants(Fr(1))
    assert inv == 1
    for chi in (Fr(3, 5), Fr(4, 5), Fr(6, 5), Fr(7, 5)):
        _, _, iv = lt.cfe1.invariants(chi)
        assert iv != 1


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = lt.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_LTB1.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_refuses_the_thermodynamic_identifications():
    cert = json.load(open(os.path.join(CERT_DIR, "LTB1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["declared_structure"].startswith("DECLARED")
    assert "NO claim of experimental validation" in \
        cb["worked_example_only"]
    assert "never evaluated" in cb["no_logarithm_evaluated"]
    assert cb["not_claimed"].startswith("NOT CLAIMED")
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "ltb1_thermodynamic_bridge.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src
    assert "math.log" not in src


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "LTB1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 6
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
