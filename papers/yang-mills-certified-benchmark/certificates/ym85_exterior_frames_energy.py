"""YM85: exact exterior-frame and physical-energy controls, Python 3.12.

The finite Z2 model and rational configurations are comparison controls.
Actual vacuum norm equalities and energy domains are proved in the paper.
Default/--check are read-only; --write changes only this result and pin.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import sys

import ym84_gauge_resolved_scaling as previous

p = previous.prev
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE / 'YM85_RESULT.json'
PIN = HERE / 'EXPECTED_YM85.sha256'
SOURCES = HERE / 'YM85_SOURCE_MANIFEST.json'
BASE_COMMIT = '51d3de6ba1e9fb37c2e836453361e911b6d9f003'
BASE_TREE = '4de17bca837e3ab1a40efc3c871cb6fae7c5da56'
PARENT_HASH = '51396ef9db120b46555f7430775379734286b38b150bbc92c5851ac47731ec9a'
MODEL = dict(
    actual_positive_vacuum_exterior_frame_construction=True,
    physical_conditional_norm_equals_unrestricted_with_complete_exterior=True,
    physical_partial_floor_equals_conditional_with_complete_exterior=True,
    cubic_single_link_physical_row_test_improved_by_Gauss_restriction=False,
    actual_dressed_plaquette_partial_energy_matrix=True,
    partial_marker_stiffness_is_full_physical_gap=False,
    first_frame_marker_has_finite_full_energy=True,
    second_and_third_normalized_frame_markers_have_finite_full_energy=False,
    finite_group_comparison_is_SU2_vacuum=False,
    finite_configuration_checks_are_universal_proofs=False,
    actual_density_floor_uniform_in_volume_or_coupling=False,
    unconditioned_response_used_as_conditioned_response=False,
    full_lattice_weak_coupling_gap_proved=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def vector(values, length):
    out = tuple(p.exact(x) for x in values)
    if len(out) != length:
        raise ValueError('Wrong vector dimension')
    return out


def dot(a, b):
    if len(a) != len(b):
        raise ValueError('Dot dimensions disagree')
    return sum((x * y for x, y in zip(a, b)), F(0))


def cross(a, b):
    a, b = vector(a, 3), vector(b, 3)
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def plus(a, b, factor=1):
    if len(a) != len(b):
        raise ValueError('Vector dimensions disagree')
    return tuple(x + factor * y for x, y in zip(a, b))


def scale(a, factor):
    return tuple(factor * x for x in a)


def qmul(a, b):
    a, b = vector(a, 4), vector(b, 4)
    return (a[0]*b[0]-dot(a[1:], b[1:]),
            *plus(plus(scale(b[1:], a[0]), scale(a[1:], b[0])), cross(a[1:], b[1:])))


def qdagger(a):
    a = vector(a, 4)
    return (a[0], *scale(a[1:], -1))


def unit(a):
    a = vector(a, 4)
    if dot(a, a) != 1:
        raise ValueError('An exact unit quaternion is required')
    return a


ONE = unit((1, 0, 0, 0))
UNITS = tuple(tuple(F(i == j) for i in range(4)) for j in range(4))
QA = unit((F(3, 5), F(4, 5), 0, 0))
QB = unit((0, F(3, 5), F(4, 5), 0))
QX = unit((F(1, 2), F(1, 2), F(1, 2), F(1, 2)))


def rational_sqrt(value):
    """Exact rational-frame fixtures only; general analytic frames may be irrational."""
    value = p.nonnegative(value)
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    if a*a != value.numerator or b*b != value.denominator:
        raise ValueError('This exact fixture needs a rational square root')
    return F(a, b)


def body_frame(a, b):
    a, b = unit(a), unit(b)
    r = rational_sqrt(dot(a[1:], a[1:]))
    if not r:
        raise ValueError('First reference holonomy is central')
    n = scale(a[1:], 1/r)
    perpendicular = plus(b[1:], n, -dot(n, b[1:]))
    s = rational_sqrt(dot(perpendicular, perpendicular))
    if not s:
        raise ValueError('Reference holonomies have parallel imaginary parts')
    second = scale(perpendicular, 1/s)
    return dict(axes=(n, second, cross(n, second)), r=r, s=s)


def canonical(q, frame):
    q = unit(q)
    return (q[0], *(dot(axis, q[1:]) for axis in frame['axes']))


def frame_differential(a, b, da, db):
    a, b, da, db = unit(a), unit(b), vector(da, 4), vector(db, 4)
    frame = body_frame(a, b)
    n, second, _ = frame['axes']
    dn = scale(plus(da[1:], n, -dot(n, da[1:])), 1/frame['r'])
    dperp = plus(plus(db[1:], dn, -dot(n, b[1:])), n,
                 -(dot(dn, b[1:]) + dot(n, db[1:])))
    dsecond = scale(plus(dperp, second, -dot(second, dperp)), 1/frame['s'])
    return (dn, dsecond, plus(cross(dn, second), cross(n, dsecond)))


def native_direction(q, axis):
    p.size(axis)
    if axis > 3:
        raise ValueError('SU(2) has three tangent axes')
    return scale(qmul(UNITS[axis], unit(q)), F(1, 2))


def pointwise_energy(a, b, x):
    frame = body_frame(a, b)
    x = unit(x)
    marker = canonical(x, frame)[1:]
    target_derivatives = [[dot(axis, native_direction(x, k)[1:]) for axis in frame['axes']]
                          for k in (1, 2, 3)]
    # Four distinct plaquette links each contribute the same native Gram matrix.
    gram = [[4 * sum(row[i]*row[j] for row in target_derivatives)
             for j in range(3)] for i in range(3)]
    assert gram == [[F(i == j)-marker[i]*marker[j] for j in range(3)] for i in range(3)]
    b_derivatives = [[dot(axis, x[1:]) for axis in frame_differential(a, b, (0,)*4, native_direction(b, k))]
                     for k in (1, 2, 3)]
    cost = [sum(row[i]**2 for row in b_derivatives) for i in range(3)]
    assert cost == [0, marker[2]**2/(4*frame['s']**2), marker[1]**2/(4*frame['s']**2)]
    return dict(markers=marker, frame_r=frame['r'], frame_s=frame['s'],
                plaquette_partial_gram=gram, exterior_chord_energy=cost,
                finite_configuration_control_not_vacuum_expectation=True)


def lattice(side):
    p.size(side, 4)
    vertices = tuple(product(range(side), repeat=3))
    edges = tuple((v, axis) for v in vertices for axis in range(3))
    def shift(v, axis, step=1):
        return tuple((x + (step if i == axis else 0)) % side for i, x in enumerate(v))
    ends = {e: (e[0], shift(*e)) for e in edges}
    return dict(side=side, vertices=vertices, edges=edges, ends=ends)


def shifted(v, axis, step, side):
    return tuple((x + (step if i == axis else 0)) % side for i, x in enumerate(v))


def detours(edge, side):
    v, axis = edge
    end = shifted(v, axis, 1, side)
    paths = []
    for transverse in range(3):
        if transverse == axis:
            continue
        paths.append((((v, transverse), 1), ((shifted(v, transverse, 1, side), axis), 1),
                      ((end, transverse), -1)))
        down = shifted(v, transverse, -1, side)
        paths.append((((down, transverse), -1), ((down, axis), 1),
                      ((shifted(end, transverse, -1, side), transverse), 1)))
    return paths


def exterior_tree(graph, removed):
    removed = set(removed)
    if not removed <= set(graph['edges']):
        raise ValueError('Removed set contains nonedges')
    adjacency = {v: [] for v in graph['vertices']}
    for edge in graph['edges']:
        if edge not in removed:
            left, right = graph['ends'][edge]
            adjacency[left].append((right, edge, 1))
            adjacency[right].append((left, edge, -1))
    root = graph['vertices'][0]
    parent, queue = {root: None}, deque([root])
    while queue:
        v = queue.popleft()
        for w, edge, sign in adjacency[v]:
            if w not in parent:
                parent[w] = (v, edge, sign)
                queue.append(w)
    if len(parent) != len(graph['vertices']):
        raise ValueError('Exterior graph is disconnected')
    tree = {item[1] for item in parent.values() if item is not None}
    chords = tuple(e for e in graph['edges'] if e not in tree and e not in removed)
    if len(chords) < 2:
        raise ValueError('Exterior has fewer than two independent cycles')
    return dict(root=root, parent=parent, tree=tree, chords=chords)


def tree_transports(tree, field):
    transports = {tree['root']: ONE}
    for v, entry in tree['parent'].items():
        if entry is not None:
            parent, edge, sign = entry
            q = unit(field[edge])
            transports[v] = qmul(transports[parent], q if sign == 1 else qdagger(q))
    return transports


def tree_holonomy(edge, graph, transports, field):
    left, right = graph['ends'][edge]
    return qmul(qmul(transports[left], field[edge]), qdagger(transports[right]))


def gauge_transform(graph, field, gauges):
    return {edge: qmul(qmul(unit(gauges[left]), unit(field[edge])), qdagger(unit(gauges[right])))
            for edge, (left, right) in graph['ends'].items()}


def lattice_frame_control(side=4):
    graph = lattice(side)
    origin = (0, 0, 0)
    removed = ((origin, 0), ((1, 0, 0), 1), ((0, 1, 0), 0), (origin, 1))
    tree = exterior_tree(graph, removed)
    ea, eb = tree['chords'][:2]
    field = {edge: ONE for edge in graph['edges']}
    field[ea], field[eb], field[removed[0]] = QA, QB, QX
    palette = (QX, QA, QB, unit((F(5, 13), 0, F(12, 13), 0)))
    gauges = {v: palette[(sum(v)+2*v[0]) % len(palette)] for v in graph['vertices']}
    changed = gauge_transform(graph, field, gauges)
    original_t, changed_t = tree_transports(tree, field), tree_transports(tree, changed)
    frame = body_frame(tree_holonomy(ea, graph, original_t, field), tree_holonomy(eb, graph, original_t, field))
    new_frame = body_frame(tree_holonomy(ea, graph, changed_t, changed), tree_holonomy(eb, graph, changed_t, changed))
    for edge in graph['edges']:
        before = canonical(tree_holonomy(edge, graph, original_t, field), frame)
        after = canonical(tree_holonomy(edge, graph, changed_t, changed), new_frame)
        assert before == after
    for edge in removed:
        paths = detours(edge, side)
        assert any(not set(item[0] for item in path) & set(removed) for path in paths)
    return dict(side=side, vertices=len(graph['vertices']), links=len(graph['edges']),
                removed_links=4, exterior_tree_edges=len(tree['tree']),
                exterior_cycle_rank=len(tree['chords']),
                all_canonical_links_gauge_invariant=True,
                full_measure_frame_proof_is_analytic=True)


def finite_reference_model(t_minus=F(1, 4), t_plus=F(7, 8), rare_probability=F(1, 100)):
    """Z2 acts by (x,y,r,z)->(-x,-y,-r,z); r supplies a complete frame."""
    ts = {-1: p.exact(t_minus), 1: p.exact(t_plus)}
    prob = p.positive(rare_probability)
    if prob >= 1 or any(abs(t) >= 1 for t in ts.values()):
        raise ValueError('All sixteen states must have positive weight')
    states = list(product((-1, 1), repeat=4))
    weights = [(prob if z == 1 else 1-prob) * (1+ts[z]*x*y)/8 for x,y,r,z in states]
    index = {s: i for i, s in enumerate(states)}
    gauge = [[F(0) for _ in states] for _ in states]
    for i, (x,y,r,z) in enumerate(states):
        gauge[i][i] += F(1, 2)
        gauge[i][index[-x,-y,-r,z]] += F(1, 2)
    px, py = (p.conditional_projection(states, weights, [i]) for i in (0, 1))
    pxy = p.conditional_projection(states, weights, [0, 1])
    response = p.add(p.multiply(px, py), pxy, -1)
    maximum = max(abs(t) for t in ts.values())
    # Weighted adjoint, so the PSD test certifies the full operator norm upper bound.
    adjoint = [[weights[j]*response[j][i]/weights[i] for j in range(16)] for i in range(16)]
    bound = p.add([[maximum**2*x for x in row] for row in p.identity(16)], p.multiply(adjoint, response), -1)
    assert p.psd([[weights[i]*x for x in row] for i,row in enumerate(bound)])
    for z, t in ts.items():
        bare = [F(x) if zz == z else F(0) for x,y,r,zz in states]
        dressed = [F(x*r) if zz == z else F(0) for x,y,r,zz in states]
        partner = [F(y*r) if zz == z else F(0) for x,y,r,zz in states]
        assert not any(previous.apply(gauge, bare))
        assert previous.apply(gauge, dressed) == dressed
        assert previous.apply(response, dressed) == [t*x for x in partner]
    return dict(states=states, weights=weights, gauge=gauge, response=response,
                ts=ts, full_norm=maximum, physical_norm=maximum,
                rare_probability=prob, kind='finite_Z2_reference_comparison_not_SU2_vacuum')


def finite_reference_control(prob):
    data = finite_reference_model(rare_probability=prob)
    return {key: data[key] for key in ('full_norm', 'physical_norm', 'rare_probability', 'kind')}


def sphere_moment(powers):
    """Even coordinate powers on Haar S^3; powers are the half-exponents."""
    if len(powers) != 4:
        raise ValueError('Four half-exponents are required')
    for n in powers:
        p.size(n, 0)
    numerator = F(1)
    for n in powers:
        for k in range(n):
            numerator *= F(1, 2)+k
    denominator = F(1)
    for k in range(sum(powers)):
        denominator *= 2+k
    return numerator/denominator


def log_two_interval(terms=12):
    p.size(terms)
    z = F(1, 3)
    lower = sum((2*z**(2*k+1)/(2*k+1) for k in range(terms)), F(0))
    tail = 2*z**(2*terms+1)/((2*terms+1)*(1-z*z))
    return dict(lower=lower, upper=lower+tail)


def singularity_control(level):
    p.size(level)
    logarithm = log_two_interval()
    return dict(epsilon=F(1, 2**level),
                Haar_exterior_energy_interval={key: level*value/8 for key,value in logarithm.items()},
                actual_cutoff_energy_coefficient_of_minimum_density={key: level*value/16 for key,value in logarithm.items()},
                positive_series_Haar_energy_lower=sum((F(1, n) for n in range(1, level+1)), F(0))/16,
                actual_density_floor_numerically_supplied=False,
                cutoff_Haar_region_probability_lower=F(1, 2))


def stiffness_control(sqrt_theta):
    s = p.positive(sqrt_theta)
    if s < 5:
        raise ValueError('The actual concentration input requires theta >= 25')
    delta = F(9, 4)/s
    lower = (1-2*delta)/(2*delta)
    assert lower == 2*s/9-1 and lower > 0
    return dict(theta=s*s, actual_delta_upper=delta,
                specified_marker_partial_quotient_lower=lower,
                simpler_sqrt_theta_over_nine_valid=s >= 9,
                all_physical_states_bound=False, full_gap_bound=False)


def controls():
    rotations = (ONE, QX, QA, QB)
    points = []
    for g in rotations:
        conjugate = lambda q: qmul(qmul(g, q), qdagger(g))
        points.append(pointwise_energy(conjugate(QA), conjugate(QB), conjugate(QX)))
    singular_points = []
    for denominator in (2, 4, 16, 64):
        t = F(1, denominator)
        b = (0, (1-t*t)/(1+t*t), 2*t/(1+t*t), 0)
        singular_points.append(pointwise_energy(QA, b, QX))
    return dict(lattice_frames=[lattice_frame_control(n) for n in (4, 5)],
                exact_native_grams=points, approaching_reference_degeneracy=singular_points,
                finite_reference_norms=[finite_reference_control(prob) for prob in (F(1, 2), F(1, 100), F(1, 1000000))],
                logarithmic_energy_budgets=[singularity_control(n) for n in (1, 4, 16, 64)],
                actual_marker_stiffness=[stiffness_control(s) for s in (5, 9, 16, 100, 1000)],
                first_axis_Haar_inverse_squared_radius=F(2),
                second_axis_Haar_inverse_squared_radius_finite=False)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage'] != 'YM85' or manifest['base_commit'] != BASE_COMMIT or manifest['base_tree'] != BASE_TREE:
        raise ValueError('YM85 source base mismatch')
    if len(manifest['upstream_sha256']) != 368 or len(manifest['local_inputs']) != 4:
        raise ValueError('YM85 source manifest shape mismatch')
    if previous.PIN.read_text().strip() != PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text())) != PARENT_HASH:
        raise ValueError('YM84 result/pin mismatch')
    checked = {}
    for path, expected in manifest['upstream_sha256'].items():
        actual = sha256((ROOT/path).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError('Changed upstream input: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2] != (3, 12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return p.encode(dict(certificate_type='YM85_EXTERIOR_FRAMES_AND_PHYSICAL_ENERGY', verdict='PASS',
                         runtime_policy='Python 3.12 only', model=MODEL,
                         input_sha256=source_checks(), controls=controls(),
                         claim_status='ACTUAL_CONDITIONAL_NORM_EQUALITY__MARKER_ENERGY_DOMAIN__CONTINUUM_OPEN'))


def check(cert, result=RESULT, pin=PIN):
    if cert.get('model') != MODEL or cert.get('certificate_type') != 'YM85_EXTERIOR_FRAMES_AND_PHYSICAL_ENERGY':
        raise ValueError('YM85 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != digest:
        raise ValueError('YM85 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument('--write', action='store_true')
    choice.add_argument('--check', action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert, sort_keys=True, indent=2)+'\n')
        PIN.write_text(digest+'\n')
    else:
        check(cert)
    print('YM85 PASS', digest)
    print('Actual exterior-frame norm identities and marker energy costs: proved; 4D continuum gap: open')


if __name__ == '__main__':
    main()
