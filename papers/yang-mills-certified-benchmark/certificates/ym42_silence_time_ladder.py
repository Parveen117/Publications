"""YM-42: structural silence contraction at t=2,3, with exact scope controls.

Reuses YM-37's carrier and T53 elimination, YM-41's native character
projectors, and YM-38's quadratic arithmetic. No new operator engine.
See ../YM42_SILENCE_TIME_LADDER.md for the general proofs and the distinction
between spatial Cauchy convergence and an operator bound in time.

Default execution verifies committed evidence without writing it. --write
is the explicit evidence-generation mode; --check is equivalent to default.
All verdicts use Fraction arithmetic. Decimal upper bounds are directed.
"""

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import sys
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))

from ym1_certified_gap import canonical_sha
from ym31_uniform_floor import fj, lam, rnd_down
import ym37_space_transfer as y37
import ym38_native_uniformity as y38
import ym41_tools_dock as y41

GRID = (F(1, 8), F(1, 4), F(1, 2))
CONTENTS = (F(0), F(1, 2), F(1))
RESULT = HERE / "YM42_RESULT.json"
PIN = HERE / "EXPECTED_YM42.sha256"
SOURCES = HERE / "YM42_SOURCE_PINS.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def upper_decimal(x, places=12):
    scale = 10 ** places
    q = -((-x.numerator * scale) // x.denominator)
    sign = "-" if q < 0 else ""
    q = abs(q)
    return f"{sign}{q // scale}.{q % scale:0{places}d}"


def sym(A):
    return bool(A) and all(len(row) == len(A) for row in A) and all(
        A[i][j] == A[j][i] for i in range(len(A)) for j in range(len(A)))


def ceiling_checks(M, B, sigma):
    """Two cut-square sign checks equivalent to a doubled norm ceiling.

    The sum/difference congruence and Schur substitution are proved in
    YM42-T1. Neither eigenvectors nor a spectral theorem are used.
    """
    if not sym(M) or not sym(B) or sigma <= 0:
        return False, {}
    n = len(M)
    if len(B) != n:
        return False, {}
    pos_b = y37.inertia(B)
    signs = [y37.inertia([[sigma * B[i][j] + sign * M[i][j]
                         for j in range(n)] for i in range(n)])
             for sign in (-1, 1)]
    ok = pos_b == (n, 0, 0) and all(v[1] == 0 for v in signs)
    return ok, {"B": list(pos_b), "sigma_B_minus_M": list(signs[0]),
                "sigma_B_plus_M": list(signs[1])}


def silence_restriction(M, B, e):
    """Exact congruence by u_i=e_i-(B_ie/B_ee)e_e, no fitted vector."""
    W = B[e][e]
    inds = [i for i in range(len(B)) if i != e]
    c = [B[i][e] / W for i in inds]
    Mp = [[M[i][j] - c[a] * M[e][j] - c[b] * M[i][e]
           + c[a] * c[b] * M[e][e]
           for b, j in enumerate(inds)] for a, i in enumerate(inds)]
    Bp = [[B[i][j] - B[i][e] * B[e][j] / W for j in inds] for i in inds]
    return Mp, Bp


def structural_budget(theta, W, D, sigma):
    """Theorem YM42-T2, all inequalities rational, including ball radius^2."""
    if theta <= 0 or W <= 0 or D < 0 or sigma < 0:
        return {"accepted": False}
    d = D / (W * theta * theta)
    s = sigma / theta
    if 1 - 2 * d <= 0:
        return {"accepted": False}
    beta = (s + 2 * d) / (1 - 2 * d)
    ball_left = 2 * s + 4 * d
    floor = theta * (1 - 2 * d)
    return {"accepted": ball_left < 1 and beta < 1,
            "d": d, "s": s, "beta": beta, "ball_left": ball_left,
            "floor": floor, "radius_squared": 4 * D / (theta * theta),
            "tail_prefactor_squared": D / (theta * theta * (1 - beta) ** 2)
            if beta != 1 else None}


def factorial_moment(e):
    """Independent factorial form, not YM-37's double-factorial algorithm."""
    if any(v % 2 for v in e):
        return F(0)
    b = [v // 2 for v in e]
    out = F(1, math.factorial(sum(b) + 1))
    for v in b:
        out *= F(math.factorial(2 * v), 4 ** v * math.factorial(v))
    return out


def grading_control():
    """The degree<=2 formula is consumed only after native projector checks."""
    mons = [e for e in itertools.product(range(3), repeat=4) if sum(e) <= 2]
    G = [[y37.moment4(tuple(a+b for a,b in zip(e, f))) for f in mons] for e in mons]
    assert y37.inertia(G) == (14, 0, 1)
    relation = [F(-1) if sum(e) == 0 else F(1) if 2 in e else F(0) for e in mons]
    assert y37.matvec(G, relation) == [F(0)] * len(mons)
    # The sole radical is sum x_i^2 - 1. Norm-zero degree<=2 differences
    # below therefore vanish by the declared sphere relation itself.
    coef = {F(0): F(2), F(1, 2): F(1, 3), F(1): F(1, 7)}
    checks = 2
    def add(p, q, scale=F(1)):
        return y37.padd(p, q, scale)
    def project(p, c):
        out = {}
        for e, v in p.items():
            out = add(out, y41.char_component(e, c), v)
        return out
    def same_on_sphere(p, q):
        diff = add(p, q, F(-1))
        return y37.integrate(y37.pmul(diff, diff), 1) == 0
    for e in mons:
        parts = {c: y41.char_component(e, c) for c in CONTENTS}
        total, native_inv = {}, {}
        for c in CONTENTS:
            total = add(total, parts[c])
            native_inv = add(native_inv, parts[c], 1 / coef[c])
            for other in CONTENTS:
                assert same_on_sphere(project(parts[c], other),
                                      parts[c] if other == c else {})
                checks += 1
        assert same_on_sphere(total, {e: F(1)})
        assert same_on_sphere(native_inv, y37.cinv_monomial(e, coef))
        checks += 2
    for e in itertools.product(range(7), repeat=4):
        if sum(e) <= 12:
            assert y37.moment4(e) == factorial_moment(e)
            checks += 1
    return {"monomials": len(mons), "Gram_inertia": [14,0,1],
            "radical": "sum x_i^2 - 1", "exact_equalities": checks}


def rail_gram(basis, f, t):
    cb = [y37.cinv(b, f, t) for b in basis]
    n = len(basis)
    # Compute both triangles independently: symmetry is an actual check.
    return [[y37.integrate(y37.pmul(basis[i], cb[j]), t)
             for j in range(n)] for i in range(n)]


def finite_iterate_controls(M, B, e, theta, r, D, W, budget, horizon=3):
    """Tests the derived recurrence/tail on exact iterates; not its proof."""
    n = len(M)
    w = [F(i == e) for i in range(n)]
    z, b, x = [F(0)] * n, F(1), w
    states = [z]
    steps = 0
    for p in range(horizon):
        xnext = y37.tau_apply(M, B, x)
        bnext = y38.bq(B, w, xnext) / W
        floor_actual = bnext / b
        denominator = theta + y38.bq(B, r, z) / W
        assert floor_actual == denominator >= budget["floor"] > 0
        tz = y37.tau_apply(M, B, z)
        wt = y38.bq(B, w, tz) / W
        az = [tz[i] - wt * w[i] for i in range(n)]
        zn = [xnext[i] / bnext - w[i] for i in range(n)]
        assert zn == [(r[i] + az[i]) / denominator for i in range(n)]
        assert y38.bq(B, w, zn) == 0
        assert y37.bnorm2(B, zn) <= budget["radius_squared"]
        step = [zn[i] - z[i] for i in range(n)]
        assert y37.bnorm2(B, step) <= D / theta**2 * budget["beta"]**(2*p)
        x, b, z = xnext, bnext, zn
        states.append(z)
        steps += 5
    for p in range(horizon):
        for q in range(p + 1, horizon + 1):
            diff = [states[q][i] - states[p][i] for i in range(n)]
            assert y37.bnorm2(B, diff) <= budget["tail_prefactor_squared"] * budget["beta"]**(2*p)
            steps += 1
    return steps


def time_blindness_control():
    """Positive finite time chains with identical perfect spatial contraction.

    Independent spatial columns => tau_t=1*mu_t, rank one, beta_space=0.
    In time P_delta has mean-zero ratio 1-delta, arbitrarily close to one.
    This is not a counterexample to the declared YM model; it refutes the
    inference from spatial contraction alone to a uniform time gap.
    """
    rows = []
    for delta in (F(1, 2), F(1, 8), F(1, 64), F(1, 1024)):
        P = [[1 - delta/2, delta/2], [delta/2, 1 - delta/2]]
        assert y37.inertia(P) == (2, 0, 0)
        assert y37.matvec(P, [F(1), F(1)]) == [F(1), F(1)]
        f = [F(1), F(-1)]
        assert y37.matvec(P, f) == [(1-delta)*v for v in f]
        for length in (2, 3, 4):
            paths = list(itertools.product(range(2), repeat=length))
            weights = []
            for path in paths:
                v = F(1, 2)
                for a, b in zip(path, path[1:]):
                    v *= P[a][b]
                weights.append(v)
            assert all(v > 0 for v in weights) and sum(weights) == 1
            corr = sum(w * f[path[0]] * f[path[-1]]
                       for w, path in zip(weights, paths))
            assert corr == (1-delta)**(length-1)
            # Gram B=diag(mu), M=mu mu^T: tau=1 mu^T exactly.
            # On the B-orthogonal zero-mean space it is zero, all t.
            probes = ([F(i == j) for i in range(len(paths))]
                      for j in (0, len(paths)-1))
            for probe in probes:
                mean = sum(a*b for a,b in zip(weights, probe))
                centered = [v-mean for v in probe]
                assert sum(a*b for a,b in zip(weights, centered)) == 0
            assert corr != 0  # false identification beta_space=ratio_time rejected
        rows.append({"delta": str(delta), "beta_space": "0",
                     "time_ratio": str(1-delta)})
    return {"finite_positive_models": rows,
            "universal_argument": "0<delta<1; beta_space=0 for all t, but rho_time=1-delta approaches 1",
            "scope": "refutes implication from spatial contraction alone; not a counterexample to NG or the YM family"}


def negative_controls():
    I = [[F(1), F(0)], [F(0), F(1)]]
    heavy = [[F(5), F(0)], [F(0), F(0)]]
    negative = [[F(-5), F(0)], [F(0), F(0)]]
    assert not ceiling_checks(heavy, I, F(1))[0]
    assert not ceiling_checks(negative, I, F(1))[0]
    assert not ceiling_checks([[F(0), F(1)], [F(0), F(0)]], I, F(1))[0]
    assert not ceiling_checks(I, [[F(1), F(0)], [F(0), F(-1)]], F(1))[0]
    assert not structural_budget(F(1), F(1), F(1, 3), F(1, 10))["accepted"]
    assert not structural_budget(F(1), F(1), F(0), F(5))["accepted"]
    assert not structural_budget(F(0), F(1), F(0), F(0))["accepted"]
    assert not structural_budget(F(1), F(0), F(0), F(0))["accepted"]
    return {"heavy_positive_complement": True, "heavy_negative_complement": True,
            "nonsymmetric_pencil": True, "indefinite_Gram": True,
            "noninvariant_ball": True, "noncontracting_complement": True,
            "zero_floor": True, "zero_silence_norm": True}


def verify_sources():
    pins = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in pins["upstream_sha256"].items():
        actual = digest(ROOT / path)
        if actual != expected:
            raise ValueError(f"upstream source pin changed: {path}")
        checked[path] = actual
    for path in pins["local_inputs"]:
        checked[path] = digest(ROOT / path)
    checked[str(SOURCES.relative_to(ROOT))] = digest(SOURCES)
    return checked


def run():
    if not __debug__:
        raise RuntimeError("YM42 verification requires assertions; optimized Python is refused")
    inputs = verify_sources()
    grading = grading_control()
    negatives = negative_controls()
    blind = time_blindness_control()
    kcoef = {F(0): F(1), F(1, 2): rnd_down(lam(1).lo), F(1): rnd_down(lam(2).lo)}
    seed = {c: rnd_down(fj(int(2*c), GRID[-1]).lo) for c in CONTENTS}
    grid = {}
    cross_checks = 0
    for t in (2, 3):
        basis, M, _, seed_B, e = y37.build_pencil(t, seed, kcoef, with_ins=False)
        n = len(basis)
        assert n == {2: 8, 3: 52}[t]
        assert sym(M) and M[e][e] == 1
        for kap in GRID:
            fcoef = {c: rnd_down(fj(int(2*c), kap).lo) for c in CONTENTS}
            B = rail_gram(basis, fcoef, t)
            if kap == GRID[-1]:
                assert B == seed_B
            assert sym(B) and y37.inertia(B) == (n, 0, 0)
            W = B[e][e]
            assert W == 1 / fcoef[F(0)]**t
            theta = M[e][e] / W
            tw = y37.solve(B, [row[e] for row in M])
            assert y37.matvec(B, tw) == [row[e] for row in M]
            w = [F(i == e) for i in range(n)]
            r = [v - theta*w[i] for i,v in enumerate(tw)]
            assert y38.bq(B, w, r) == 0
            wrong = [v - theta*F(101,100)*w[i] for i,v in enumerate(tw)]
            assert y38.bq(B, w, wrong) != 0
            D = y37.bnorm2(B, r)
            assert D > 0
            assert D == sum(M[e][i]*tw[i] for i in range(n)) - M[e][e]**2/W
            Mp, Bp = silence_restriction(M, B, e)
            sigma = theta * kap**2 / 4  # declared coarse rational ceiling, not fitted eigen-data
            valid, signs = ceiling_checks(Mp, Bp, sigma)
            assert valid
            budget = structural_budget(theta, W, D, sigma)
            assert budget["accepted"]
            # Independent prior construction and doubled norm route, t=2.
            if t == 2:
                _, old_Mp, old_Bp = y38.restricted(M, B, w, W)
                assert old_Mp == Mp and old_Bp == Bp
                cols = [y37.solve(Bp, [Mp[i][j] for i in range(n-1)]) for j in range(n-1)]
                doubled = [[sum(Mp[i][k]*cols[j][k] for k in range(n-1))
                            - sigma**2*Bp[i][j] for j in range(n-1)] for i in range(n-1)]
                assert y37.inertia(doubled)[0] == 0
                cross_checks += 3
            steps = finite_iterate_controls(M, B, e, theta, r, D, W, budget)
            grid.setdefault(str(kap), {})[f"t{t}"] = {
                "dimension": n, "complement_dimension": n-1,
                "face_coefficients": {str(c):str(v) for c,v in fcoef.items()},
                "W": str(W), "theta": str(theta), "defect_squared": str(D),
                "sigma_over_theta": str(kap**2/4),
                "budget_exact": {k:str(v) for k,v in budget.items() if k != "accepted"},
                "beta_upper_decimal": upper_decimal(budget["beta"]),
                "ball_left_upper_decimal": upper_decimal(budget["ball_left"]),
                "sign_checks": signs, "iterate_checks": steps,
                "wrong_theta_rejected": True, "accepted": True}
    return {"certificate_type": "YM42_STRUCTURAL_SILENCE_TIME_LADDER",
            "verdict": "PASS", "input_sha256": inputs,
            "rung_coefficients": {str(c):str(v) for c,v in kcoef.items()},
            "grid": grid, "native_grading_and_independent_moments": grading,
            "prior_route_cross_checks": cross_checks,
            "negative_controls": negatives, "time_blindness_control": blind,
            "claims": {"structural_ball_and_Cauchy_tail_all_spatial_steps": "PROVED_UNDER_EXPLICIT_QUADRATIC_HYPOTHESES",
                       "t2_t3_declared_rational_grid": "CERTIFIED",
                       "all_time_rails": "OPEN",
                       "t3_infinite_content_release": "OPEN",
                       "operator_time_gap": "OPEN",
                       "native_measure_dictionary": "OPEN",
                       "continuum_or_Clay": "NOT_CLAIMED"}}


def check_evidence(cert, result=RESULT, pin=PIN):
    sha = canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError("fresh YM42 certificate differs from committed evidence")
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    cert = run()
    sha = canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, indent=2, sort_keys=True)+"\n")
        PIN.write_text(sha+"\n")
    else:
        check_evidence(cert)
    print("YM42 PASS", sha)
    print(json.dumps({k:{t:v["beta_upper_decimal"] for t,v in row.items()}
                      for k,row in cert["grid"].items()}, sort_keys=True))


if __name__ == "__main__":
    main()
