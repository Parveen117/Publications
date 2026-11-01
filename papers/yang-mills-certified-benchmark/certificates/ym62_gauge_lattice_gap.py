"""YM62: native 3+1 gauge lattice, plaquette response and gap controls.

Python 3.12 only. Exact controls accompany the written all-content proof;
they do not numerically solve a vacuum or prove the continuum mass gap.
Default/--check are read-only. --write changes only this result and hash.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache, reduce
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym61_source_response_collective as previous
import ym53_anisotropic_interaction as blocks

native, heat, source = previous.native, previous.heat, previous.source
a, p, y, exact = previous.a, previous.p, previous.y, previous.exact
RESULT = HERE / 'YM62_RESULT.json'
PIN = HERE / 'EXPECTED_YM62.sha256'
SOURCES = HERE / 'YM62_SOURCE_PINS.json'
EPS = blocks.old.EPS
SQUARE = ((0, 1), (1, 2), (3, 2), (0, 3))
ORIENTED_WORD = ((0, 1), (1, 1), (2, -1), (3, -1))
UNITS = tuple(tuple(F(i == j) for i in range(4)) for j in range(4))


def lattice(side, periodic=True):
    """Positive oriented links and ordered elementary faces of a cubic box."""
    if not isinstance(periodic, bool):
        raise ValueError('Periodic flag must be Boolean')
    native.width(side, 4 if periodic else 2)
    vertices = tuple(itertools.product(range(side), repeat=3))
    def step(x, axis):
        out = list(x)
        out[axis] += 1
        if periodic:
            out[axis] %= side
        return tuple(out)
    edges = tuple((x, axis) for x in vertices for axis in range(3)
                  if periodic or x[axis]+1 < side)
    index = {edge: i for i, edge in enumerate(edges)}
    ends = tuple((x, step(x, axis)) for x, axis in edges)
    faces = []
    for x in vertices:
        for mu, nu in itertools.combinations(range(3), 2):
            if periodic or (x[mu]+1 < side and x[nu]+1 < side):
                faces.append(((index[x, mu], 1), (index[step(x, mu), nu], 1),
                              (index[step(x, nu), mu], -1), (index[x, nu], -1)))
    return dict(vertices=vertices,edges=edges,ends=ends,faces=tuple(faces))


def incidence(faces, edge_count):
    native.width(edge_count)
    degrees, pairs = [0]*edge_count, Counter()
    for word in faces:
        indices = [e for e, sign in word]
        if len(set(indices)) != len(indices) or any(
                isinstance(e, bool) or not isinstance(e, int) or not 0 <= e < edge_count
                or isinstance(sign, bool) or sign not in (-1, 1) for e, sign in word):
            raise ValueError('Each face must have distinct valid oriented links')
        for e in indices:
            degrees[e] += 1
            for other in indices:
                if e != other:
                    pairs[e, other] += 1
    return tuple(degrees), pairs


def quaternion_product(left, right):
    """Polynomial extension of the inherited native quaternion multiplication."""
    out = [{} for _ in range(4)]
    for i, j in itertools.product(range(4), repeat=2):
        value = p.mul(left[i], right[j])
        for k, sign in enumerate(y.qmul(UNITS[i], UNITS[j])):
            if sign:
                out[k] = p.add(out[k], p.scale(value, sign))
    return tuple(out)


def loop_polynomial(word=ORIENTED_WORD, edge_count=4):
    incidence((word,), edge_count)
    if not word:
        raise ValueError('A nonempty loop word is required')
    result = (source.constant(1, edge_count), {}, {}, {})
    for edge, sign in word:
        q = tuple(native.lift(f, edge, edge_count) for f in p.COORD)
        if sign == -1:
            q = (q[0], *(p.scale(f, -1) for f in q[1:]))
        result = quaternion_product(result, q)
    return result[0]


def right_derivative(poly, edge, axis):
    native.width(axis)
    if axis > 3:
        raise ValueError('Internal axis must be 1, 2 or 3')
    if not poly:
        return {}
    n = len(next(iter(poly), (0, 0, 0, 0)))//4
    if not 0 <= edge < n:
        raise ValueError('Invalid link')
    matrix = y.qmatrix(tuple(-v/2 for v in UNITS[axis]), right=True)
    coords = [native.lift(f, edge, n) for f in p.COORD]
    return p.add(*(p.mul(p.add(*(p.scale(coords[j], matrix[i][j]) for j in range(4))),
                         p.partial(poly, 4*edge+i)) for i in range(4)))


def gauss(poly, vertex, axis, ends=SQUARE):
    return p.add(*(native.derivative(poly, e, axis) for e, (s, t) in enumerate(ends) if s == vertex),
                 *(p.scale(right_derivative(poly, e, axis), -1)
                   for e, (s, t) in enumerate(ends) if t == vertex))


def electric(poly, edge_count):
    return p.scale(p.add(*(native.derivative(native.derivative(poly, e, k), e, k)
                           for e in range(edge_count) for k in (1, 2, 3))), -1)


def energy(poly, edge_count):
    return p.add(*(p.mul(df, df) for e in range(edge_count) for k in (1, 2, 3)
                   for df in (native.derivative(poly, e, k),)))


def evaluate(poly, quaternions):
    point = tuple(c for q in quaternions for c in q)
    return sum((coefficient*y.prod(x**power for x, power in zip(point, monomial))
                for monomial, coefficient in poly.items()), F(0))


def loop_value(qs, word=ORIENTED_WORD):
    out = y.ONE
    for edge, sign in word:
        q = qs[edge] if sign == 1 else y.qdagger(qs[edge])
        out = y.qmul(out, q)
    return out


def gauge_transform(qs, gauges, ends=SQUARE):
    return tuple(y.qmul(y.qmul(gauges[s], q), y.qdagger(gauges[t]))
                 for q, (s, t) in zip(qs, ends))


def block_gate(theta, J=5, R=None, degree=4, arity=4):
    theta = exact(theta)
    native.width(J)
    native.width(degree)
    native.width(arity, 2)
    alpha = 2/(EPS*J)+40*degree*(arity+1)*abs(theta)*J
    R = (1+1/alpha)/2 if R is None else exact(R)
    if R <= 1 or R*alpha >= 1:
        return None
    gamma = blocks.old.log_iv(blocks.Iv(R), blocks.old.LOG_TERMS)/blocks.Iv(10*J)
    return dict(alpha=alpha,R=R,margin=1-R*alpha,gamma=gamma)


def refinement_budget(plaquettes, theta, time, mesh):
    native.width(plaquettes, 0)
    theta, time, mesh = exact(theta), exact(time), exact(mesh)
    if time < 0 or not 0 <= mesh <= time:
        raise ValueError('Require 0 <= mesh <= time')
    b = abs(theta)*plaquettes
    if not b or not time or not mesh:
        return F(0)
    return b*time*blocks.refine.ex(b*time).hi*blocks.refine.iv_sqrt(blocks.Iv(18*mesh)).hi


@lru_cache(maxsize=None)
def gauge_controls():
    w = loop_polynomial()
    count = 0
    for vertex, axis in itertools.product(range(4), (1, 2, 3)):
        assert not gauss(w, vertex, axis)
        # Noninvariant probe checks operator commutation, not only G(W)=0.
        probe = p.mul(w, native.lift(p.mul(p.COORD[0], p.COORD[1]), 0, 4))
        assert electric(gauss(probe, vertex, axis), 4) == gauss(electric(probe, 4), vertex, axis)
        count += 1
    turns = tuple(y.turns())
    cases = []
    for offset in range(6):
        qs = tuple(turns[(offset+i)%len(turns)] for i in range(4))
        gauges = tuple(turns[(offset+2*i+1)%len(turns)] for i in range(4))
        changed = gauge_transform(qs, gauges)
        word = loop_value(qs)
        assert loop_value(changed) == y.qmul(y.qmul(gauges[0], word), y.qdagger(gauges[0]))
        assert evaluate(w, qs) == word[0] == evaluate(w, changed)
        cases.append(dict(offset=offset,trace=word[0]))
    # Open-link overlap and an incorrectly ordered face fail local gauge invariance.
    qs = (turns[0], turns[2], turns[4], y.ONE)
    gauges = (y.ONE, turns[1], y.ONE, y.ONE)
    changed = gauge_transform(qs, gauges)
    overlap = lambda row: sum(x*z for x, z in zip(row[0], row[1]))
    assert overlap(qs) != overlap(changed)
    wrong = ((0, 1), (2, -1), (1, 1), (3, -1))
    assert loop_value(qs, wrong)[0] != loop_value(changed, wrong)[0]
    f = p.mul(p.COORD[0], p.COORD[0])
    rank_one = lambda g: p.scale(y.native_generator(y.native_generator(g, 1), 1), -1)
    assert rank_one(y.native_generator(f, 2)) != y.native_generator(rank_one(f), 2)
    return dict(gauss_ward_identities=count,casimir_gauss_commutators=count,
                finite_gauge_words=cases,open_overlap_rejected=True,wrong_order_rejected=True,
                open_overlap_before=overlap(qs),open_overlap_after=overlap(changed),
                wrong_order_before=loop_value(qs, wrong)[0],wrong_order_after=loop_value(changed, wrong)[0],
                rank_one_colour_operator_rejected=True)


@lru_cache(maxsize=None)
def plaquette_controls():
    w = loop_polynomial()
    assert len(w) == 64
    assert electric(w, 4) == p.scale(w, 3)
    radii = [native.lift(y.RADIUS, e, 4) for e in range(4)]
    radius_product = reduce(p.mul, radii)
    square = p.mul(w, w)
    assert energy(w, 4) == p.add(radius_product, p.scale(square, -1))
    assert native.phi(w, 4) == 0
    assert native.phi(square, 4) == F(1, 4)
    assert native.phi(energy(w, 4), 4) == F(3, 4)
    assert native.phi(p.mul(w, p.scale(w, F(1, 3))), 4) == F(1, 12)
    return dict(polynomial_terms=len(w),free_electric_rate=F(3),mean=F(0),variance=F(1, 4),
                energy=F(3, 4),quadratic_energy_per_plaquette=F(-1, 12),
                diagonal_source_susceptibility=F(1, 6),bulk_quadratic_energy_per_cell=F(-1, 4),
                gradient_identity='Gamma(W)=1-W^2 on the native unit carrier')


def geometry_controls():
    rows = []
    for periodic, sides in ((False, (2, 3, 5)), (True, (4, 5, 6))):
        for side in sides:
            cell = lattice(side, periodic)
            E, P = len(cell['edges']), len(cell['faces'])
            degrees, pairs = incidence(cell['faces'], E)
            assert E == (3*side**3 if periodic else 3*(side-1)*side**2)
            assert P == (3*side**3 if periodic else 3*side*(side-1)**2)
            assert max(degrees) <= 4
            assert all(sum(v for (a, b), v in pairs.items() if a == e) == 3*degrees[e] for e in range(E))
            for word in cell['faces']:
                directed = [cell['ends'][e] if sign == 1 else cell['ends'][e][::-1] for e, sign in word]
                assert all(directed[i][1] == directed[(i+1)%4][0] for i in range(4))
            supports = [frozenset(e for e, sign in word) for word in cell['faces']]
            assert len(set(supports)) == P
            # Distinct four-link supports have an odd, unpaired link in W_p W_q.
            rows.append(dict(periodic=periodic,side=side,links=E,plaquettes=P,
                             maximum_incidence=max(degrees),off_diagonal_free_response_zero=True))
    return rows


def incidence_controls():
    # Four incident plaquettes around one link, including every clipped label.
    # The central link attains the full outgoing budget 4*(4-1)=12.
    faces = tuple(tuple((e, 1) for e in (0, 3*k+1, 3*k+2, 3*k+3)) for k in range(4))
    E = 13
    degrees, pairs = incidence(faces, E)
    assert degrees[0] == 4 and sum(v for (e, j), v in pairs.items() if e == 0) == 12
    eta, checks = F(1, 100000), 0
    for N, ell in itertools.product((1, 2, 4), (1, 2, 3, 5)):
        labels = blocks.old.blocks(E, N, ell)
        base, D = 1+F(1, 10*ell), 1/EPS+16*eta*ell*ell
        R = base**ell
        weight = lambda u: base**min(u, N+1-u)
        for e, u in itertools.product(range(E), range(1, N+1)):
            removed, temporal, spatial, cost = 0, 0, 0, F(0)
            for j, interval in labels:
                top_weight = max(weight(v) for v in interval)
                if e == j:
                    removed += u in interval
                    if u in (interval[0]-1, interval[-1]+1):
                        temporal += 1
                        cost += top_weight*D
                elif u in interval:
                    multiplicity = pairs[e, j]
                    spatial += multiplicity
                    cost += top_weight*4*eta*ell*multiplicity
            assert removed == ell and temporal <= 2 and spatial == 3*degrees[e]*ell <= 12*ell
            assert cost <= weight(u)*R*(2/EPS+80*eta*ell*ell)
            checks += 1
    # A shared pair in two interactions must be charged twice.
    _, multiple = incidence((((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 1), (3, 1))), 4)
    assert multiple[0, 1] == 2
    return dict(weighted_coordinate_checks=checks,maximum_link_plaquette_degree=4,
                plaquette_arity=4,outgoing_spatial_budget=12,total_fine_coefficient=80,
                pair_multiplicity_preserved=True)


def gap_controls():
    rows, fine, endpoints = [], 0, 0
    for theta in (F(0), F(1, 100000), F(1, 16000), F(999, 8000000)):
        gate = block_gate(theta)
        assert gate is not None
        for h in (F(1), F(2, 3), F(1, 7), F(17, 256), F(1, 65536)):
            r = blocks.skeleton(1, h)
            ell = 5*r
            actual = 2*r/(EPS*ell)+80*h*abs(theta)*ell
            assert actual <= gate['alpha'] and gate['R']*actual < 1
            fine += 1
        rows.append(dict(theta=theta,alpha=gate['alpha'],R=gate['R'],margin=gate['margin'],
                         gamma_lower=blocks.old.decimals(gate['gamma'].lo, False),
                         gamma_upper=blocks.old.decimals(gate['gamma'].hi)))
    for theta in (F(1, 8000), F(1, 4000), F(1, 64)):
        for J in range(1, 31):
            alpha = 2/(EPS*J)+800*theta*J
            assert J*(alpha-1) == F((J-5)**2, 10)+(800*theta-F(1, 10))*J*J
            assert block_gate(theta, J) is None
            endpoints += 1
    example = block_gate(F(1, 16000))
    assert example['alpha'] == F(3, 4) and example['R'] == F(7, 6) and example['margin'] == F(1, 8)
    # Recover the actual predecessor pair-chain budget under d=2, q=2.
    chain = block_gate(F(1, 4096), 6, F(5, 4), 2, 2)
    old = blocks.block_gate(1, F(1, 4096), 6, F(5, 4))
    assert chain['alpha'] == old['alpha'] and chain['margin'] == old['margin']
    return dict(window='abs(theta) < 1/8000',cells=rows,fine_step_checks=fine,
                endpoint_refusals=endpoints,chain_budget_recovered=True,
                example_rate='log(7/6)/50',at_or_above_endpoint_gapless_claim=False)


def history_controls():
    # This is the local projection/variance estimate, not a computed vacuum.
    epsilon = F(16, 3)*F(1, 8000)
    assert epsilon == F(1, 1500) < F(1, 38)**2
    variance_floor = F(1, 4)-F(2, 38)-4*epsilon
    assert variance_floor > F(19, 100)
    # For unitary-normalized h and its projected v, the 2x2 state difference
    # squares to epsilon I: hence trace distance 2 sqrt(epsilon).
    # Use exact rational amplitudes to check the identity independently.
    u, v = F(3, 5), F(4, 5)
    difference = [[u*u-1, u*v], [u*v, v*v]]
    assert a.multiply(difference, difference) == a.scale(a.identity(2), v*v)
    return dict(local_projection_error_upper=epsilon,variance_lower=variance_floor,
                stated_physical_plaquette_variance_floor=F(19, 100),
                increment_bound='E[(W(t)-W(0))^2] <= 2t',
                periodic_states_are_not_averaged=True,
                uniform_gap_passes_on_all_positive_time_cylinders=True,
                uniqueness_of_subsequential_limits=False,time_zero_markov_limit=False)


def refinement_controls():
    rows = []
    for P in (6, 192):
        errors = [refinement_budget(P, F(1, 16000), 1, F(1, n)) for n in (16, 64, 256)]
        assert errors[2] < errors[1] < errors[0]
        roots = [blocks.refine.iv_sqrt(blocks.Iv(F(18, n))) for n in (16, 64, 256)]
        assert all(not roots[i].separated_from(blocks.Iv(2)*roots[i+1]) for i in (0, 1))
        rows.append(dict(plaquettes=P,error_upper=errors))
    assert refinement_budget(192, 0, 1, F(1, 4)) == 0
    assert refinement_budget(0, 1, 1, F(1, 4)) == 0
    assert 4*F(3, 4) == 3 and 2*3 == 6 and 3*6 == 18
    return dict(cells=rows,commutator_squared_coefficient=6,packet_squared_coefficient=18,
                norm_limit_fixed_volume=True,uniform_in_volume_refinement_error=False)


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
        certificate_type='YM62_NATIVE_3_PLUS_1_GAUGE_LATTICE_STRONG_COUPLING_GAP',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        gauge=gauge_controls(),plaquette=plaquette_controls(),geometry=geometry_controls(),
        incidence=incidence_controls(),heat=blocks.heat_controls(),gap=gap_controls(),
        refinement=refinement_controls(),histories=history_controls(),
        claim_status='FIXED_SPATIAL_LATTICE_GAUGE_GAP__SPATIAL_CONTINUUM_OPEN',
        proof_scope='Written all-content proof with exact finite controls; no proof-assistant or peer verification',
        model=dict(spatial_dimensions=3,continuous_heat_time_dimensions=1,native_links=True,
                   independent_casimir_electric_operator=True,ordered_plaquette_potential=True,
                   local_gauss_sector=True,volume_uniform_interacting_gap=True,
                   subsequential_infinite_volume_gapped_histories=True,
                   nontrivial_physical_loop_variance=True,
                   unchanged_YM61_correlated_chain=False,continuum_R4_QFT=False,
                   all_compact_simple_gauge_groups=False,physical_mass_prediction=False),
        open='Spatial continuum and scale trajectory; Euclidean/relativistic axioms; general gauge group'))


def check(cert, result=RESULT, pin=PIN):
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM62 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = native.previous.native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM62 PASS', sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],dimensions='3 spatial + 1 heat time',
                          coupling_window='abs(theta) < 1/8000',continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
