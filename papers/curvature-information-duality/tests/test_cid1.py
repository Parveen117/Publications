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


c = _load("cid1", "cid1_curvature_information_duality.py")


# ---------------- T1: the information metric ----------------

def test_metric_is_symmetric_and_positive_definite():
    g = c.covariance(c.P_POINT, c.STATS)
    assert g == c.mT(g)
    assert c.is_pd(g)


def test_metric_is_the_covariance_by_an_independent_route():
    p, stats = c.P_POINT, c.STATS
    m = [c.expect(p, T) for T in stats]
    for i in range(2):
        for j in range(2):
            centered = sum(
                p[x] * (stats[i][x] - m[i]) * (stats[j][x] - m[j])
                for x in range(len(p)))
            assert c.covariance(p, stats)[i][j] == centered


def test_affine_dependence_makes_metric_degenerate_with_exact_kernel():
    dep = tuple(c.STATS[0][x] + c.STATS[1][x] for x in range(4))
    g3 = c.covariance(c.P_POINT, c.STATS + (dep,))
    assert c.is_psd(g3) and not c.is_pd(g3)
    assert c.mv(g3, [Fr(1), Fr(1), Fr(-1)]) == [Fr(0)] * 3


def test_constants_lie_in_the_kernel():
    const = tuple(Fr(1) for _ in c.OUTCOMES)
    g = c.covariance(c.P_POINT, c.STATS + (const,))
    assert c.mv(g, [Fr(0), Fr(0), Fr(1)]) == [Fr(0)] * 3


# ---------------- T2: the recognition ledger ----------------

def test_total_covariance_identity_exact_for_every_partition():
    g = c.covariance(c.P_POINT, c.STATS)
    for blocks in ([(0, 1), (2, 3)], [(0,), (1,), (2,), (3,)],
                   [(0, 1, 2, 3)], [(0,), (1, 2), (3,)],
                   [(0, 3), (1,), (2,)]):
        _, _, between, within = c.coarse_grain(c.P_POINT, c.STATS, blocks)
        assert c.madd(between, within) == g
        assert c.is_psd(between) and c.is_psd(within)


def test_monotonicity_recognition_never_manufactures_information():
    g = c.covariance(c.P_POINT, c.STATS)
    for blocks in ([(0, 1), (2, 3)], [(0,), (1, 2, 3)],
                   [(0, 1, 2, 3)]):
        _, _, between, _ = c.coarse_grain(c.P_POINT, c.STATS, blocks)
        assert c.is_psd(c.msub(g, between))


def test_sufficient_partition_discards_nothing():
    _, _, between, within = c.coarse_grain(
        c.P_POINT, c.STATS, [(0,), (1,), (2,), (3,)])
    assert within == [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
    assert between == c.covariance(c.P_POINT, c.STATS)


def test_total_collapse_recognizes_nothing():
    _, _, between, within = c.coarse_grain(
        c.P_POINT, c.STATS, [(0, 1, 2, 3)])
    assert between == [[Fr(0), Fr(0)], [Fr(0), Fr(0)]]
    assert within == c.covariance(c.P_POINT, c.STATS)


# ---------------- T3: the obstruction is metric-free ----------------

def test_all_declared_metrics_are_positive_definite():
    for name, g in c.METRICS.items():
        assert c.is_pd(g), name


def test_residue_ladder_independent_of_metric():
    vals = []
    for chi in c.CHIS:
        w = c.cfeu.witness_vec(chi)
        r = c.cfeu.circulation_of(w, c.LOOP)
        assert r == c.cfe1.circulation(c.LOOP, chi)
        assert r == c.cfe1.BETA * (chi - 1) * c.cfeu.AREA
        vals.append(r)
    assert vals == [Fr(-32), Fr(-16), Fr(0), Fr(16), Fr(32)]


# ---------------- T4: the split is metric-dependent ----------------

def test_split_is_genuine_for_every_metric():
    basis = c.cfeu.nullspace(c.cfeu.D1)
    w = c.cfeu.witness_vec(Fr(7, 5))
    for name, g in c.METRICS.items():
        _, closed, rem = c.hodge_split(w, g, basis)
        assert c.cfeu.mv(c.cfeu.D1, closed) == [Fr(0)]
        assert c.cfeu.circulation_of(closed, c.LOOP) == 0
        for k in basis:
            assert c.inner(rem, k, g) == 0
        assert [x + y for x, y in zip(closed, rem)] == w


def test_closed_part_integrates_to_explicit_potential():
    basis = c.cfeu.nullspace(c.cfeu.D1)
    w = c.cfeu.witness_vec(Fr(3, 5))
    for g in c.METRICS.values():
        _, closed, _ = c.hodge_split(w, g, basis)
        phi = c.cfeu.potential_of_closed(closed)
        assert c.cfeu.mv(c.cfeu.D0, phi) == closed


def test_splits_genuinely_differ_but_residue_does_not():
    basis = c.cfeu.nullspace(c.cfeu.D1)
    w = c.cfeu.witness_vec(Fr(7, 5))
    seen = set()
    for g in c.METRICS.values():
        _, closed, rem = c.hodge_split(w, g, basis)
        assert c.cfeu.circulation_of(rem, c.LOOP) == Fr(32)
        seen.add(tuple(rem))
    assert len(seen) > 1


# ---------------- T5: one obstruction, two presentations ----------------

def test_complement_rank_is_one_for_every_metric():
    basis = c.cfeu.nullspace(c.cfeu.D1)
    for name, g in c.METRICS.items():
        comp = []
        for j in range(6):
            e = [Fr(0)] * 6
            e[j] = Fr(1)
            _, _, rem = c.hodge_split(e, g, basis)
            comp.append(rem)
        assert c.cfeu.rank(comp) == 1, name


def test_generators_differ_but_all_carry_nonzero_residue():
    basis = c.cfeu.nullspace(c.cfeu.D1)
    gens = set()
    for g in c.METRICS.values():
        comp = []
        for j in range(6):
            e = [Fr(0)] * 6
            e[j] = Fr(1)
            _, _, rem = c.hodge_split(e, g, basis)
            comp.append(rem)
        gen = next(x for x in comp if any(y != 0 for y in x))
        assert c.cfeu.circulation_of(gen, c.LOOP) != 0
        gens.add(tuple(gen))
    assert len(gens) > 1


def test_degenerate_form_is_not_invertible_control():
    d = [[Fr(1), Fr(1)], [Fr(1), Fr(1)]]
    assert c.is_psd(d) and not c.is_pd(d)
    assert c.det(d) == 0


# ---------------- integrity and boundary ----------------

def test_certificate_pin_matches_regeneration():
    cert = c.build_certificate()
    payload = json.dumps(cert, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(payload.encode()).hexdigest()
    with open(os.path.join(CERT_DIR, "EXPECTED_CID1.sha256")) as f:
        assert digest == f.read().strip()


def test_claim_boundary_declares_family_form_and_source():
    cert = json.load(open(os.path.join(CERT_DIR, "CID1_RESULT.json")))
    cb = cert["claim_boundary"]
    assert cb["family_form_declared"].startswith("DECLARED")
    assert cb["not_claimed"].startswith("NOT CLAIMED")
    assert "has NOT been read" in cb["source"]
    assert cb["RH_K0_L0"] == "not touched"


def test_no_transcendental_machinery_in_certificate_source():
    src = open(os.path.join(
        CERT_DIR, "cid1_curvature_information_duality.py")).read()
    assert "import math" not in src
    assert "numpy" not in src
    assert "float(" not in src
    assert "exp(" not in src.replace("expect(", "").replace(
        "exponential", "")


def test_all_blocks_pass():
    cert = json.load(open(os.path.join(CERT_DIR, "CID1_RESULT.json")))
    blocks = [k for k in cert if k.startswith("T")]
    assert len(blocks) == 5
    for k in blocks:
        assert cert[k]["verdict"] == "PASS"
