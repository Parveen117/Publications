"""YM64: vacuum-adapted spatial cuts and a uniformly bounded response tower.

Python 3.12 only. Native controls and explicitly labelled rational operator
fixtures check the written proof; fixtures do not replace the YM62 vacuum.
Default/--check are read-only. --write changes only this result and hash.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym63_spatial_blocking_scale as previous

native, p, y, a = previous.native, previous.p, previous.y, previous.a
model, exact = previous.previous, previous.exact
RESULT = HERE / 'YM64_RESULT.json'
PIN = HERE / 'EXPECTED_YM64.sha256'
SOURCES = HERE / 'YM64_SOURCE_PINS.json'


def matrix(rows):
    out = [[exact(v) for v in row] for row in rows]
    if not out or not out[0] or any(len(row) != len(out[0]) for row in out):
        raise ValueError('Nonempty rectangular exact matrix required')
    return out


def sub(A, rows, columns=None):
    columns = rows if columns is None else columns
    return [[A[i][j] for j in columns] for i in rows]


def schur(A, retained):
    A = matrix(A)
    if len(A) != len(A[0]):
        raise ValueError('Square matrix required')
    retained = tuple(retained)
    if not retained or len(set(retained)) != len(retained) or any(
            isinstance(i, bool) or not isinstance(i, int) or not 0 <= i < len(A)
            for i in retained):
        raise ValueError('Distinct valid retained indices required')
    hidden = tuple(i for i in range(len(A)) if i not in retained)
    C = sub(A, retained)
    if not hidden:
        return C
    B, D, Bt = sub(A, retained, hidden), sub(A, hidden), sub(A, hidden, retained)
    return a.add(C, a.multiply(a.multiply(B, a.inverse(D)), Bt), -1)


def effective(H, retained, damping):
    damping = previous.positive(damping)
    return a.add(schur(a.add(H, a.identity(len(H)), damping), retained),
                 a.identity(len(retained)), -damping)


def moment(B, D, order):
    native.width(order, 0)
    B, D = matrix(B), matrix(D)
    if len(D) != len(D[0]) or len(B[0]) != len(D) or not a.psd(D):
        raise ValueError('Compatible positive hidden matrix required')
    return a.multiply(a.multiply(B, a.power(a.inverse(D), order+1)), a.transpose(B))


def self_energy(B, D, damping):
    damping = exact(damping)
    shifted = a.add(D, a.identity(len(D)), damping)
    if not a.psd(shifted):
        raise ValueError('Positive hidden resolvent required')
    return a.multiply(a.multiply(B, a.inverse(shifted)), a.transpose(B))


def tower_budget(gap, damping, order):
    gap, damping = previous.positive(gap), exact(damping)
    native.width(order, 0)
    ratio = abs(damping)/gap
    if ratio >= 1:
        raise ValueError('The convergent tower requires abs(damping) < gap')
    return dict(ratio=ratio,remainder_factor=ratio**(order+1)/(1-ratio),
                coefficient_factor=gap**(-order))


def hankel(moments, size, shift=0):
    native.width(size)
    native.width(shift, 0)
    width = len(moments[0])
    if 2*(size-1)+shift >= len(moments):
        raise ValueError('Insufficient moments')
    return [[moments[i+j+shift][r][c]
             for j in range(size) for c in range(width)]
            for i in range(size) for r in range(width)]


def fixture():
    """Six-dimensional operator fixture, not a Yang--Mills truncation."""
    C = matrix(((0, 0, 0), (0, 3, F(1, 5)), (0, F(1, 5), 4)))
    B = matrix(((0, 0, 0), (F(1, 2), F(1, 3), 0), (F(1, 4), 0, F(1, 5))))
    D = a.diag((2, 3, 5))
    H = [C[i]+B[i] for i in range(3)]
    H += [a.transpose(B)[i]+D[i] for i in range(3)]
    assert a.psd(a.add(H, a.diag((0, 1, 1, 1, 1, 1)), -1))
    return H, C, B, D


def partition(weights, groups):
    weights = tuple(previous.positive(w) for w in weights)
    groups = tuple(tuple(group) for group in groups)
    n = len(weights)
    flat = [i for group in groups for i in group]
    if not groups or any(not group for group in groups) or any(
            isinstance(i, bool) or not isinstance(i, int) for i in flat
            ) or sorted(flat) != list(range(n)):
        raise ValueError('A disjoint exhaustive nonempty partition is required')
    I = [[F(i in group) for group in groups] for i in range(n)]
    masses = tuple(sum((weights[i] for i in group), F(0)) for group in groups)
    J = [[weights[i]/masses[k] if i in group else F(0) for i in range(n)]
         for k, group in enumerate(groups)]
    return I, J, masses


def weighted_effective(H, I, J, damping):
    damping = previous.positive(damping)
    G = a.multiply(a.multiply(J, a.inverse(a.add(H, a.identity(len(H)), damping))), I)
    return a.add(a.inverse(G), a.identity(len(J)), -damping)


def conditional_controls():
    weights = tuple(F(v, 59) for v in (1, 2, 3, 5, 7, 11, 13, 17))
    assert sum(weights) == 1
    I, J, coarse = partition(weights, ((0, 1), (2, 3), (4, 5), (6, 7)))
    U, Z, twice = partition(coarse, ((0, 1), (2, 3)))
    P = a.multiply(I, J)
    R = a.multiply(a.multiply(I, U), a.multiply(Z, J))
    assert a.multiply(J, I) == a.identity(4)
    assert a.multiply(P, P) == P and a.multiply(R, P) == R and a.multiply(P, R) == R
    assert a.multiply(a.diag(weights), P) == a.multiply(a.transpose(P), a.diag(weights))
    assert all(sum(row) == 1 for row in P)
    H = [[F(i == j)-weights[j] for j in range(8)] for i in range(8)]
    for i, j, rate in ((0, 2, F(1, 17)), (1, 7, F(1, 23)), (4, 6, F(1, 19))):
        H[i][i] += rate/weights[i]
        H[j][j] += rate/weights[j]
        H[i][j] -= rate/weights[i]
        H[j][i] -= rate/weights[j]
    assert all(sum(row) == 0 for row in H)
    assert a.psd(a.multiply(a.diag(weights), H))
    rows = []
    for s in (F(1, 7), F(1), F(3)):
        first = weighted_effective(H, I, J, s)
        nested = weighted_effective(first, U, Z, s)
        direct = weighted_effective(H, a.multiply(I, U), a.multiply(Z, J), s)
        assert nested == direct and all(sum(row) == 0 for row in direct)
        rows.append(dict(damping=s,coarse_operator=direct))
    return dict(fixture='Weighted finite reversible chain; not the YM vacuum',
                fine_masses=weights,coarse_masses=coarse,twice_masses=twice,
                weighted_orthogonal_projection=True,conditional_tower=True,
                exact_nested_resolvents=rows,unweighted_averaging_differs=P[0][0] != F(1, 2))


def operator_controls():
    H, C, B, D = fixture()
    moments = [moment(B, D, n) for n in range(8)]
    M0 = moments[0]
    assert a.psd(a.add(a.add(C, M0, -1), a.diag((0, 1, 1)), -1))
    for n, value in enumerate(moments):
        assert a.psd(value) and a.psd(a.add(M0, value, -1))
    block, shifted = hankel(moments, 3), hankel(moments, 3, 1)
    assert a.psd(block) and a.psd(shifted) and a.psd(a.add(block, shifted, -1))
    commutator = a.add(a.multiply(moments[0], moments[1]), a.multiply(moments[1], moments[0]), -1)
    assert any(v for row in commutator for v in row)
    rows = []
    for s in (F(1, 7), F(1, 2), F(2)):
        sigma = self_energy(B, D, s)
        diff = a.add(M0, sigma, -1)
        assert a.psd(diff) and a.psd(a.add(a.scale(M0, s/(1+s)), diff, -1))
        nested = effective(effective(H, (0, 1, 2), s), (0, 1), s)
        direct = effective(H, (0, 1), s)
        assert nested == direct
        G = sub(a.inverse(a.add(H, a.identity(6), s)), (0, 1, 2))
        assert a.inverse(G) == a.add(a.add(C, a.identity(3), s), sigma, -1)
        trunc = [[F(0)]*3 for _ in range(3)]
        for n in range(6):
            trunc = a.add(trunc, moments[n], (-s)**n)
            signed_error = a.scale(a.add(sigma, trunc, -1), (-1)**(n+1))
            factor = s**(n+1)/(1+s)
            assert a.psd(signed_error) and a.psd(a.add(a.scale(M0, factor), signed_error, -1))
        reset = effective(a.add(C, M0, -1), (0, 1), s)
        assert reset != direct
        rows.append(dict(damping=s,direct=direct,static_reset=reset))
    kappa = F(7, 3)
    for n in range(6):
        assert moment(a.scale(B, kappa), a.scale(D, kappa), n) == a.scale(moments[n], kappa**(1-n))
    assert effective(a.scale(H, kappa), (0, 1), F(1, 2)) == a.scale(
        effective(H, (0, 1), F(1, 2)/kappa), kappa)
    tails = []
    for ticks in (1, 2, 5):
        decay = a.diag([F(1, 2)**(ticks*int(D[i][i]))/D[i][i] for i in range(3)])
        tail = a.multiply(a.multiply(B, decay), a.transpose(B))
        assert a.psd(tail) and a.psd(a.add(a.scale(M0, F(1, 2)**(2*ticks)), tail, -1))
        tails.append(dict(time_in_log2_units=ticks,tail=tail,relative_bound=F(1, 2)**(2*ticks)))
    return dict(fixture='Rational positive operator with an exact vacuum; not a YM spin truncation',
                full_gap_lower=1,hidden_gap_lower=2,first_four_moments=moments[:4],
                hankel_inertia=p.y37.inertia(block),shifted_hankel_positive=True,
                gap_localizing_hankel_positive=True,noncommuting_moments=commutator,
                nested_blocks_and_reset_counterexamples=rows,
                positive_remainder_checks=18,clock_covariance_checks=7,
                integrated_memory_tails=tails)


def obstruction_controls():
    overlaps = []
    for t in (F(1, 2), F(1, 4), F(1, 16), F(1, 64)):
        h = (2*t/(1+t*t), (1-t*t)/(1+t*t))
        H = [[F(i == j)-h[i]*h[j] for j in range(2)] for i in range(2)]
        assert sum(v*v for v in h) == 1 and H[1][1] == h[0]**2
        overlaps.append(dict(parameter=t,free_cut_hidden_gap=H[1][1],adapted_hidden_gap=1))
    long_range = []
    for n in (8, 16, 32):
        H = [[F(i == j)-F(1, n) for j in range(n)] for i in range(n)]
        distant = n//2
        H[0][0] += 1
        H[distant][distant] += 1
        H[0][distant] -= 1
        H[distant][0] -= 1
        base = [[F(i == j)-F(1, n) for j in range(n)] for i in range(n)]
        assert a.psd(a.add(H, base, -1)) and all(sum(row) == 0 for row in H)
        long_range.append(dict(sites=n,separation=distant,off_diagonal=H[0][distant],gap_lower=1))
    return dict(free_overlap_can_vanish=overlaps,
                gap_alone_does_not_imply_spatial_locality=long_range,
                examples_are_not_claims_about_the_YM_Hamiltonian=True)


def gamma_poly(f, g, width):
    return p.add(*(p.mul(native.derivative(f, edge, axis), native.derivative(g, edge, axis))
                   for edge in range(width) for axis in (1, 2, 3)))


def word_derivative(values, word, edge, axis):
    if not isinstance(axis, int) or isinstance(axis, bool) or axis not in (1, 2, 3):
        raise ValueError('A native derivative axis is required')
    out = (F(0),)*4
    unit = model.UNITS[axis]
    for index, (current, orientation) in enumerate(word):
        if current != edge:
            continue
        product = y.ONE
        for position, (link, sign) in enumerate(word):
            q = values.get(link, y.ONE)
            factor = q if sign == 1 else y.qdagger(q)
            if position == index:
                factor = (tuple(-v/2 for v in y.qmul(unit, q)) if sign == 1
                          else tuple(v/2 for v in y.qmul(y.qdagger(q), unit)))
            product = y.qmul(product, factor)
        out = tuple(x+v for x, v in zip(out, product))
    return out[0]


@lru_cache(maxsize=None)
def native_controls():
    checks = []
    quadratic = p.mul(p.COORD[0], p.COORD[1])
    for length, f, g in ((2, p.COORD[0], p.COORD[1]),
                          (2, quadratic, p.COORD[0]), (3, p.COORD[0], p.COORD[0])):
        lhs = gamma_poly(previous.pullback(f, length), previous.pullback(g, length), length)
        rhs = p.scale(previous.pullback(gamma_poly(f, g, 1), length), length)
        assert lhs == rhs
        checks.append(dict(factor=length,terms=len(lhs)))
    # Independently validate both word-derivative signs against native polynomials.
    word = model.ORIENTED_WORD
    poly = model.loop_polynomial()
    values = {i: q for i, q in enumerate(y.turns()[:4])}
    for edge in range(4):
        for axis in (1, 2, 3):
            lhs = word_derivative(values, word, edge, axis)
            rhs = model.evaluate(native.derivative(poly, edge, axis), tuple(values.values()))
            assert lhs == rhs
    # Gauge-invariant physical coarse plaquette, b=3. Its first-order leakage
    # is evaluated in the actual native link algebra, without a spin truncation.
    geometry = previous.blocked_geometry(4, 3)
    coarse_word = geometry['coarse']['faces'][0]
    fine_word = tuple((edge, sign) for coarse_edge, sign in coarse_word
                      for edge in (geometry['paths'][coarse_edge] if sign == 1
                                   else reversed(geometry['paths'][coarse_edge])))
    active = {edge for edge, sign in fine_word}
    exceptional = geometry['paths'][coarse_word[0][0]][1]
    values = {exceptional: (F(3, 5), F(4, 5), F(0), F(0))}
    adjacent = tuple(face for face in geometry['fine']['faces']
                     if active.intersection(edge for edge, sign in face))
    def response(point):
        return sum((word_derivative(point, face, edge, axis)*word_derivative(point, fine_word, edge, axis)
                    for face in adjacent for edge in active.intersection(e for e, sign in face)
                    for axis in (1, 2, 3)), F(0))
    gamma = response(values)
    assert gamma == F(16, 25)
    vertex = geometry['fine']['ends'][exceptional][0]
    gauge = (F(3, 5), F(0), F(4, 5), F(0))
    changed = {}
    for edge, (start, end) in enumerate(geometry['fine']['ends']):
        q = values.get(edge, y.ONE)
        if start == vertex:
            q = y.qmul(gauge, q)
        if end == vertex:
            q = y.qmul(q, y.qdagger(gauge))
        if q != y.ONE:
            changed[edge] = q
    assert response(changed) == gamma
    far_word = next(face for face in geometry['coarse']['faces']
                    if geometry['coarse']['edges'][face[0][0]] == ((2, 2, 2), 0))
    far_active = {e for link, sign in far_word for e in geometry['paths'][link]}
    far_faces = [face for face in geometry['fine']['faces']
                 if far_active.intersection(e for e, sign in face)]
    support = {e for face in adjacent for e, sign in face}
    far_support = {e for face in far_faces for e, sign in face}
    assert support.isdisjoint(far_support)
    return dict(pointwise_carre_du_champ_pullbacks=checks,word_derivative_sign_checks=12,
                physical_probe='One coarse plaquette represented by a 3-by-3 fine boundary',
                fine_boundary_links=len(fine_word),adjacent_plaquettes=len(adjacent),
                gamma_V_probe_at_rational_configuration=gamma,
                first_order_hidden_leakage=-F(2, 3)*gamma,
                local_gauge_transformation_preserves_response=True,
                distant_leading_response_supports_disjoint=True,
                full_vacuum_not_replaced_by_this_first_order_coefficient=True)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = native.digest(ROOT/path)
        if actual != expected:
            raise ValueError('upstream source changed: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = native.digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))] = native.digest(SOURCES)
    return checked


def run():
    native.previous.runtime_check()
    return native.previous.encode(dict(
        certificate_type='YM64_VACUUM_ADAPTED_NESTED_BLOCKS_AND_POSITIVE_RESPONSE_TOWER',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        native=native_controls(),conditional=conditional_controls(),operators=operator_controls(),
        obstructions=obstruction_controls(),
        series_budget=tower_budget(F(1, 10), F(1, 40), 4),
        claim_status='FINITE_VOLUME_EXACT_BLOCK_COMPOSITION_AND_UNIFORM_RELATIVE_SPECTRAL_BOUNDS',
        model=dict(YM62_operator_retained=True,actual_vacuum_weighted_cut=True,
                   full_constant_vacuum_retained=True,uniform_hidden_gap_from_YM62=True,
                   nested_resolvent_composition=True,positive_inverse_moment_tower=True,
                   uniform_energy_relative_memory_bounds=True,
                   local_first_order_physical_leakage=True,
                   full_spatial_locality_bound=False,uniform_absolute_memory_norm=False,
                   explicit_interacting_vacuum_computed=False,exact_scalar_RG=False,
                   continuum_R4_QFT=False,physical_mass_prediction=False),
        proof_scope='Written finite-volume form proof; finite rational fixtures are separate algebraic controls',
        open='All-order spatial locality, constructive coarse vacuum control and a nontrivial continuum trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM64 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM64 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],hidden_gap='gamma(theta)',
                          physical_first_order_leakage='-32/75',continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
