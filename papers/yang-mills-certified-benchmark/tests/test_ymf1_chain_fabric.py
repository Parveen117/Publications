import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "certificates"))

import ymf1_chain_fabric as f1  # noqa: E402


def test_verdict_and_pin():
    cert = f1.run()
    assert cert["verdict"] == "PASS"
    pin = open(os.path.join(HERE, "..", "certificates",
                            "EXPECTED_YMF1.sha256")).read().strip()
    assert f1.canonical_sha(cert) == pin


def test_quaternion_identification_is_homomorphism():
    q1 = f1.rational_unit(F(1, 2), F(-1, 3), F(2))
    q2 = f1.rational_unit(F(3), F(0), F(-1, 5))
    assert f1.meq(f1.to_matrix(q1 * q2),
                  f1.mmul(f1.to_matrix(q1), f1.to_matrix(q2)))
    assert f1.meq(f1.to_matrix(f1.QI * f1.QJ), f1.to_matrix(f1.QK))


def test_det_split_twisted_vs_untwisted():
    a, b, c, d = F(2), F(2), F(3), F(1)          # EMK-1's channel witness
    assert f1.mdet(f1.emk_block(a, b, c, d, twist=False)) == f1.GQ(8)   # 0 + 8
    assert f1.mdet(f1.emk_block(a, b, c, d, twist=True)) == f1.GQ(18)   # norm


def test_stokes_telescopes_and_tamper_breaks():
    rungs = [f1.rational_unit(F(k, 3), F(1, k + 1), F(-k, 7)) for k in range(1, 6)]
    faces = f1.ladder_faces(rungs)
    assert f1.stokes_product(faces) == rungs[0].inv() * rungs[-1]
    faces[2] = rungs[2].inv() * rungs[3].inv()
    assert f1.stokes_product(faces) != rungs[0].inv() * rungs[-1]


def test_face_residue_endpoints():
    assert f1.QONE.residue() == 0
    assert f1.Quat(-1, 0, 0, 0).residue() == 2
    h = f1.rational_unit(F(1), F(1), F(1))
    assert 0 < h.residue() < 2
