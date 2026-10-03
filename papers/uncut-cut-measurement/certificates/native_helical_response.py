"""Exact controls for NH-1--NH-4; native cut scalars, no complex or float.

The integer sheet is never reduced modulo a finite sheet count. Finite
coefficient checks accompany the written bounded-operator/infinite-chain proof.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


native = load("nh_native_scalar", HERE / "native_sources/native_summability.py")
g3 = load("nh_emkg3", ROOT.parent / "emk-recognition-geometry/certificates/emkg3_helical_sheet_memory.py")
C = native.NativeCutScalar
ZERO, ONE, IOTA = C.zero(), C.one(), C.iota()


def padd(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, ZERO) + v
        if out[k].is_zero():
            del out[k]
    return out


def pmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out = padd(out, {i + j: x * y})
    return out


def pscale(a, c):
    return {k: c * v for k, v in a.items() if not (c * v).is_zero()}


def pdagger(a):
    return {-k: v.dagger() for k, v in a.items()}


def apply_poly(a, x):
    # (S x)_n = x_(n-1); fields here have finite support.
    return pmul(a, x)


def inner(a, b):
    return sum((v.dagger() * b.get(k, ZERO) for k, v in a.items()), ZERO)


def inverse_real(a):
    n = len(a)
    aug = [[F(v) for v in row] + [F(i == j) for j in range(n)]
           for i, row in enumerate(a)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if aug[i][j])
        aug[j], aug[pivot] = aug[pivot], aug[j]
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(n):
            if i != j:
                factor = aug[i][j]
                aug[i] = [x - factor*y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def h_cycle(n, w, eta):
    h = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        h[i][i] += eta + 2*w
        h[i][(i-1) % n] -= w
        h[i][(i+1) % n] -= w
    return h


def run():
    pins = json.loads((ROOT / "NATIVE_HELICAL_SOURCE_PINS.json").read_text())
    consumed = 0
    for rec in pins["records"]:
        if rec["local_consumed_path"] is None:
            continue
        data = (ROOT / rec["local_consumed_path"]).read_bytes()
        assert hashlib.sha256(data).hexdigest() == rec["sha256"]
        assert hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == rec["git_blob_sha1"]
        consumed += 1
    assert IOTA * IOTA == -ONE and IOTA.dagger() == -IOTA
    assert IOTA.norm_square() == 1

    # Explicit faithful scalar-to-pair intertwiner, rather than an unrelated R.
    for a in (C(2, 3), C(F(1, 2), F(-2, 3)), IOTA):
        for z in (C(1, 4), C(-3, 2)):
            az = a*z
            paired = (a.radial*z.radial-a.turn*z.turn,
                      a.turn*z.radial+a.radial*z.turn)
            assert paired == (az.radial, az.turn)
            assert (IOTA*z).radial == -z.turn
            assert (IOTA*z).turn == z.radial

    # Original G3 phase, scale and integer-sheet advance, all retained.
    nbase = 3
    seed = (0, F(0), 0)

    def label(s):
        turns, j = divmod(s, nbase)
        return j, g3.rho_power(seed, turns)

    labels = [label(s) for s in range(-60, 61)]
    assert len(set(labels)) == len(labels)
    assert label(6) == (0, (0, F(4, 5), 4))
    assert label(6) != label(0)
    assert label(6)[1][0] == label(0)[1][0]

    w, eta = F(1), F(1, 2)
    h = {-1: C(-w), 0: C(eta+2*w), 1: C(-w)}
    d = {0: ONE, 1: -ONE}
    assert h == padd(pscale(pmul(pdagger(d), d), C(w)), {0: C(eta)})
    assert pdagger(h) == h
    psi = {-3: C(1, 2), 0: C(2, -1), 1: C(1, 1), 4: C(F(1, 3), 1)}
    lhs = inner(psi, apply_poly(h, psi))
    rhs = C(w)*inner(apply_poly(d, psi), apply_poly(d, psi)) + C(eta)*inner(psi, psi)
    assert lhs == rhs and lhs.turn == 0 and lhs.radial > 0
    wrong_h = dict(h); wrong_h[1] = C(w)
    assert wrong_h != pdagger(wrong_h)
    assert inner(psi, apply_poly(wrong_h, psi)) != rhs

    # Generic recurrence coefficients plus the source jump determine the
    # infinite Green function in the written proof. No truncated chain inverse.
    r, amp = F(1, 2), F(2, 3)
    assert (eta+2*w)*r-w*(1+r*r) == 0
    assert amp*(eta+2*w-2*w*r) == 1

    def green(s):
        return amp*r**abs(s)

    for s in range(-100, 101):
        assert (eta+2*w)*green(s)-w*green(s-1)-w*green(s+1) == F(s == 0)
    norm_square = amp*amp*(1+2*r*r/(1-r*r))
    total_mass = amp*(1+2*r/(1-r))
    assert norm_square == F(20, 27)
    assert total_mass == 1/eta == 2
    assert amp*(eta+2*w-2*w*F(1, 4)) != 1  # wrong decay fails

    # Exact cyclic aggregation: sum all integer-sheet copies, not finite wrap.
    periodic_origin = amp*(1+r**nbase)/(1-r**nbase)
    cyc = h_cycle(nbase, w, eta)
    inv = inverse_real(cyc)
    assert periodic_origin == inv[0][0] == F(6, 7)
    assert periodic_origin-green(0) == F(4, 21)
    assert periodic_origin != green(0)  # aggregate is not sheet-local readout
    for n in (3, 4, 5, 7):
        assert inverse_real(h_cycle(n, w, eta))[0][0] == amp*(1+r**n)/(1-r**n)

    # On finite-support fields P H = H_cycle P; this is a lawful summed
    # observable, but erases the integer sheet. No false nonclosure assertion.
    def project(x):
        return [sum((v for s, v in x.items() if s % nbase == j), ZERO)
                for j in range(nbase)]
    px = project(psi)
    assert project(apply_poly(h, psi)) == [sum((C(cyc[i][j])*px[j] for j in range(nbase)), ZERO) for i in range(nbase)]
    assert project({0: ONE}) == project({3: ONE})
    assert {0: ONE}.get(0, ZERO) != {3: ONE}.get(0, ZERO)

    # Exact bounded-operator factorial coefficients over Laurent shift algebra.
    # The infinite convergence/unitarity argument is written in NH-3.
    generator = pscale(h, -IOTA)
    assert pdagger(generator) == pscale(generator, C(-1))
    coeff = [{0: ONE}]
    for degree in range(1, 13):
        coeff.append(pscale(pmul(coeff[-1], generator), C(F(1, degree))))
    for degree in range(13):
        conv = {}
        for j in range(degree+1):
            conv = padd(conv, pmul(pdagger(coeff[j]), coeff[degree-j]))
        assert conv == ({0: ONE} if degree == 0 else {})
    bad_generator = h
    assert padd(pdagger(bad_generator), bad_generator) != {}
    wrong_square = (C(1, 1)*C(1, 1)).radial
    assert wrong_square != C(1, 1).norm_square()

    return {
        "status": "PASS_NATIVE_HELICAL_EXACT_CONTROLS",
        "consumed_source_files_checked": consumed,
        "native_scalar": "pinned NativeCutScalar, no Python complex or float",
        "retained_G3_parameters": {"alpha_mod_6": 3, "beta": "2/5", "q_integer": 2},
        "two_circuit_fibre": [0, "4/5", 4],
        "sheet_index_reduced_modulo_finite_size": False,
        "w": str(w), "eta": str(eta),
        "green_kernel": "(2/3)*(1/2)^abs(s)",
        "recurrence_points_checked": 201,
        "green_norm_square": str(norm_square),
        "green_l1_mass": str(total_mass),
        "sheet_local_response": str(green(0)),
        "three_site_sheet_summed_response": str(periodic_origin),
        "exact_hidden_sheet_contribution": str(periodic_origin-green(0)),
        "unitary_series_coefficients_checked_through_degree": 12,
        "negative_controls": ["wrong edge sign", "wrong decay", "sheet-local equals aggregate", "missing iota in generator", "square substituted for dagger norm"],
        "infinite_operator_claim_evidence": "written proof; finite checks are not a universal proof",
        "energy_law_weights_and_pinning_selected_from_iota_alone": False,
        "physical_gravity_or_hbar_derived": False,
    }


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
