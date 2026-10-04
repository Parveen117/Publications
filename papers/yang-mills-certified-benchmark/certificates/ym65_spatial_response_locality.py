"""YM65: spatial propagation, connected inverse responses and bulk control.

Python 3.12 only. Exact incidence and interval controls accompany the
written all-content proof. No finite-spin approximation of YM62 is used.
Default/--check are read-only; --write changes only this result and hash.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from math import factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym64_vacuum_adapted_tower as previous

native, p, y, a = previous.native, previous.p, previous.y, previous.a
model, exact = previous.model, previous.exact
Iv, refine = model.blocks.Iv, model.blocks.refine
RESULT = HERE / 'YM65_RESULT.json'
PIN = HERE / 'EXPECTED_YM65.sha256'
SOURCES = HERE / 'YM65_SOURCE_PINS.json'


def nonnegative(value):
    value = exact(value)
    if value < 0:
        raise ValueError('A nonnegative exact scalar is required')
    return value


def exp_negative(x):
    if not isinstance(x, Iv):
        x = Iv(nonnegative(x))
    if x.lo < 0:
        raise ValueError('Nonnegative exponent interval required')
    return Iv(p.exp_negative(x.hi).lo, p.exp_negative(x.lo).hi)


def iv_power(x, order):
    native.width(order, 0)
    out = Iv(1)
    for _ in range(order):
        out = p._r(out*x)
    return out


def parameters(theta):
    theta = exact(theta)
    gate = model.block_gate(theta)
    if gate is None:
        raise ValueError('Outside the inherited strict YM62 gap window')
    gamma = p._r(gate['gamma'])
    if gamma.lo <= 0:
        raise ArithmeticError('Positive gap is below the declared interval rounding precision')
    log_two = p._r(model.blocks.old.log_iv(Iv(2), 80))
    velocity = 64*abs(theta)
    rho = p._r(log_two*gamma/(Iv(2*velocity)+gamma))
    return dict(theta=theta,gamma=gamma,log_two=log_two,velocity=velocity,rho=rho)


def lr_upper(theta, distance, support, time):
    theta, time = exact(theta), abs(exact(time))
    native.width(distance)
    native.width(support)
    value = 64*abs(theta)*time
    if not value:
        return F(0)
    if value > 8:
        return F(2)  # The universal commutator bound is sufficient here.
    return min(F(2), F(2*support, 2**distance)*(refine.ex(value).hi-1))


def response_budget(theta, distance, support, order=0, clock=F(1)):
    native.width(distance)
    native.width(support)
    native.width(order, 0)
    clock = previous.previous.positive(clock)
    cell = parameters(theta)
    r = p._r(cell['rho']*Iv(distance))
    if r.lo <= 0:
        raise ArithmeticError('Spatial rate below interval rounding precision')
    decay = exp_negative(r)
    prefactor = p._r(Iv(1+F(2*support, 3))+Iv(1)/(Iv(3)*r))
    if order == 0:
        envelope = p._r(prefactor*decay)
    else:
        x = p._r(Iv(2)*r)
        partial = Iv(0)
        for k in range(order):
            partial = p._r(partial+iv_power(x, k)/Iv(factorial(k)))
        envelope = p._r(prefactor*decay*iv_power(x, order)/Iv(factorial(order))
                         +exp_negative(x)*partial)
    normalized_upper = min(F(1), envelope.hi)
    if not cell['theta']:
        normalized_upper = F(0)  # Exact free product-vacuum factorization.
    absolute_upper = normalized_upper/(clock*cell['gamma'].lo)**order
    return dict(distance=distance,smaller_support=support,inverse_order=order,
                rate=cell['rho'],normalized_upper=normalized_upper,
                absolute_upper=absolute_upper,clock=clock,
                cutoff_time=p._r(Iv(2)*r/(Iv(clock)*cell['gamma'])))


def geometric_moments(q):
    q = exact(q)
    if not 0 <= q < 1:
        raise ValueError('A geometric ratio in [0,1) is required')
    return (q/(1-q), q/(1-q)**2, q*(1+q)/(1-q)**3,
            q*(1+4*q+q*q)/(1-q)**4,
            q*(1+11*q+11*q*q+q**3)/(1-q)**5)


def susceptibility_sum_upper(theta):
    cell = parameters(theta)
    rho = cell['rho']
    q = exp_negative(rho).hi
    if q >= 1:
        raise ArithmeticError('Decay too small for the declared interval precision')
    S = geometric_moments(q)
    u = rho.hi
    total = 16+176*u*S[4]+(1320*u+40)*S[3]+(3300*u+300)*S[2]
    total += (2750*u+750)*S[1]+625*S[0]
    return p._r(Iv(2*total/cell['gamma'].lo)).hi


def bulk_mean_error(theta_max, side, difference_step, theta=F(0)):
    theta_max = previous.previous.positive(theta_max)
    native.width(side, 4)
    theta, difference_step = exact(theta), previous.previous.positive(difference_step)
    if abs(theta)+difference_step > theta_max:
        raise ValueError('The centered difference must stay inside the controlled interval')
    bound = susceptibility_sum_upper(theta_max)
    epsilon = 4*theta_max/side
    return min(F(2), epsilon/difference_step+bound*difference_step)


@lru_cache(maxsize=None, typed=True)
def graph(side, periodic=True):
    lattice = model.lattice(side, periodic)
    rows = [Counter() for _ in lattice['edges']]
    faces = tuple(frozenset(e for e, sign in word) for word in lattice['faces'])
    for face in faces:
        for e in face:
            rows[e].update(face)
    return lattice, tuple(rows), faces


def distances(rows, start):
    out = {e: 0 for e in start}
    queue = deque(start)
    while queue:
        e = queue.popleft()
        for f in rows[e]:
            if f not in out:
                out[f] = out[e]+1
                queue.append(f)
    return out


def advance(rows, state):
    out = Counter()
    for e, weight in state.items():
        for f, multiplicity in rows[e].items():
            out[f] += weight*multiplicity
    return out


def incidence_controls():
    geometries = []
    for side, periodic in ((4, True), (5, True), (6, True), (3, False), (5, False)):
        lattice, rows, faces = graph(side, periodic)
        maximum = max(sum(row.values()) for row in rows)
        assert maximum <= 16
        if periodic:
            assert {sum(row.values()) for row in rows} == {16}
            assert {row[e] for e, row in enumerate(rows)} == {4}
        midpoint = [tuple(2*x[k]+int(axis == k) for k in range(3)) for x, axis in lattice['edges']]
        for face in faces:
            for e, f in itertools.product(face, repeat=2):
                delta = [abs(u-v) for u, v in zip(midpoint[e], midpoint[f])]
                if periodic:
                    delta = [min(v, 2*side-v) for v in delta]
                assert sum(delta) <= 2
        origin = faces[0]
        d = distances(rows, origin)
        state = Counter({e: 1 for e in origin})
        for n in range(7):
            assert sum(state.values()) <= len(origin)*16**n
            assert all(d[e] <= n for e in state)
            if periodic:
                assert sum(state.values()) == len(origin)*16**n
            state = advance(rows, state)
        geometries.append(dict(side=side,periodic=periodic,links=len(rows),faces=len(faces),
                               max_incidence_row=maximum,graph_diameter_from_probe=max(d.values())))
    lattice, rows, faces = graph(4)
    X, Y = faces[0], faces[126]
    d = min(distances(rows, X)[e] for e in Y)
    neighbours = [tuple(j for j, other in enumerate(faces) if face & other) for face in faces]
    face_state = {i: 1 for i, face in enumerate(faces) if face & X}
    link_state = Counter({e: 1 for e in X})
    comparisons = []
    for n in range(1, 7):
        link_state = advance(rows, link_state)
        chains = sum(value for i, value in face_state.items() if faces[i] & Y)
        walks = sum(link_state[e] for e in Y)
        assert chains <= walks <= len(X)*16**n
        if n < d:
            assert chains == walks == 0
        if n == d:
            assert walks > 0
        comparisons.append(dict(order=n,plaquette_chains=chains,link_walks=walks))
        fresh = Counter()
        for i, value in face_state.items():
            for j in neighbours[i]:
                fresh[j] += value
        face_state = fresh
    return dict(geometries=geometries,independent_chain_distance=d,
                independent_chain_vs_walk_counts=comparisons,
                midpoint_distance_changes_by_at_most_one=True)


def atan_small(x, terms=32):
    x = exact(x)
    native.width(terms)
    if not 0 < x < 1:
        raise ValueError('Alternating arctangent requires 0 < x < 1')
    value = sum(((-1)**k*x**(2*k+1)/(2*k+1) for k in range(terms)), F(0))
    next_term = (-1)**terms*x**(2*terms+1)/(2*terms+1)
    return Iv(min(value, value+next_term), max(value, value+next_term))


@lru_cache(maxsize=None)
def pi_interval():
    tangent = F(1, 5)
    for _ in range(2):
        tangent = 2*tangent/(1-tangent*tangent)
    assert (tangent-F(1, 239))/(1+tangent/F(239)) == 1
    return p._r(Iv(16)*atan_small(F(1, 5))-Iv(4)*atan_small(F(1, 239)))


@lru_cache(maxsize=None, typed=True)
def erfc_interval(x):
    x = exact(x)
    if abs(x) > 4:
        raise ValueError('Gaussian finite-control argument exceeds the declared range')
    if x < 0:
        return Iv(2)-erfc_interval(-x)
    if x == 0:
        return Iv(1)
    terms = 96
    value = sum(((-1)**k*x**(2*k+1)/(factorial(k)*(2*k+1)) for k in range(terms)), F(0))
    next_term = (-1)**terms*x**(2*terms+1)/(factorial(terms)*(2*terms+1))
    integral = Iv(min(value, value+next_term), max(value, value+next_term))
    root_pi = refine.iv_sqrt(pi_interval())
    result = p._r(Iv(1)-Iv(2)*integral/root_pi)
    return Iv(max(F(0), result.lo), min(F(2), result.hi))


def scalar_filter(energy, heat_time):
    """The spectral image of the Gaussian commutator filter, alpha=1/4."""
    energy, heat_time = previous.previous.positive(energy), nonnegative(heat_time)
    if heat_time > 2*energy or energy+heat_time/2 > 4 or energy*heat_time > 8:
        raise ValueError('Outside the declared scalar filter test range')
    plus = refine.ex(-heat_time*energy)*erfc_interval(heat_time/2-energy)/Iv(2)
    minus = refine.ex(heat_time*energy)*erfc_interval(heat_time/2+energy)/Iv(2)
    return p._r(plus-minus)


def filter_controls():
    rows = []
    for energy in (F(1), F(2), F(3)):
        for b in (F(0), F(1, 2), F(1), F(2)):
            alpha = F(1, 4)
            assert -b*energy-(energy-2*alpha*b)**2/(4*alpha) == -energy**2/(4*alpha)-alpha*b*b
            assert b*energy-(energy+2*alpha*b)**2/(4*alpha) == -energy**2/(4*alpha)-alpha*b*b
            filtered = scalar_filter(energy, b)
            target = refine.ex(-b*energy)
            error = target-filtered
            bound = exp_negative(energy*energy+b*b/4)
            assert error.lo > 0 and error.hi < bound.lo
            rows.append(dict(energy=energy,heat_time=b,filtered=filtered,target=target,
                             error_upper=error.hi,proved_error_upper=bound.hi))
    # Removing the normalization exp(-alpha*b^2) genuinely changes the filter.
    correct = scalar_filter(3, 1)
    wrong = correct*refine.ex(F(1, 4))
    target = refine.ex(-3)
    assert (wrong-target).lo > exp_negative(F(37, 4)).hi
    return dict(pi_enclosure=pi_interval(),scalar_spectral_filters=rows,
                omitted_gaussian_normalization_rejected=True,
                scope='Exact interval checks of the scalar spectral filter, not a YM spectrum')


def source_controls():
    h, v = [F(3, 5), 0, 0, F(4, 5)], [F(4, 5), 0, 0, F(-3, 5)]
    outer = lambda x: [[u*w for w in x] for u in x]
    P0 = outer(h)
    H = a.add(a.scale(outer(v), 2), a.diag((0, 3, 4, 0)))
    R = a.add(a.inverse(a.add(H, P0)), P0, -1)
    A, B = a.diag((1, 1, -1, -1)), a.diag((1, -1, 1, -1))
    assert a.multiply(A, B) == a.multiply(B, A)
    mean = F(-7, 25)
    x = [u-mean*w for u, w in zip(p.apply(A, h), h)]
    z = [u-mean*w for u, w in zip(p.apply(B, h), h)]
    assert x == z and p.dot(h, x) == 0
    moments = [p.dot(x, p.apply(a.power(R, n), z)) for n in range(7)]
    assert moments == [F(576, 625)*F(1, 2)**n for n in range(7)]
    # Exact two-dimensional lower eigenvalue: 1-sqrt(1-14j/25+j^2).
    direct_negative_second_derivative = 1-F(7, 25)**2
    assert 2*moments[1] == direct_negative_second_derivative
    return dict(fixture='Two local binary readouts and an exact entangled vacuum; not a YM truncation',
                mean=mean,inverse_moments=moments,
                independent_energy_hessian=direct_negative_second_derivative,
                source_reciprocity=True,vacuum_subtraction_required=True)


def tiling_controls():
    faces = lambda L: 3*L*(L-1)**2
    rows = []
    for block in (2, 3, 4):
        for count in (2, 3, 5):
            side = block*count
            crossing = faces(side)-count**3*faces(block)
            assert 0 <= crossing <= F(6*side**3, block)
            rows.append(dict(block_side=block,blocks_per_axis=count,crossing_faces=crossing,
                             crossing_per_volume=F(crossing, side**3)))
    lattice = model.lattice(6, False)
    internal = 0
    for word in lattice['faces']:
        vertices = {vertex for edge, sign in word for vertex in lattice['ends'][edge]}
        if len({tuple(c//3 for c in vertex) for vertex in vertices}) == 1:
            internal += 1
    assert internal == 8*faces(3)
    wraps = [dict(side=L,wrap_faces=3*L**3-faces(L)) for L in (4, 5, 8, 16)]
    assert all(row['wrap_faces'] <= 6*row['side']**2 for row in wraps)
    for q in (F(0), F(1, 4), F(1, 2)):
        moments = geometric_moments(q)
        for k in range(5):
            finite = sum((F(d**k)*q**d for d in range(1, 60)), F(0))
            assert finite <= moments[k]
    return dict(tilings=rows,independent_vertex_partition_internal_faces=internal,
                periodic_wrap_faces=wraps,periodic_energy_density_error='4*abs(theta)/L',
                bulk_energy_regular='C1,1 on compact subintervals of the certified window',
                no_second_derivative_convergence_claim=True)


def conditioning_controls():
    rows = []
    for distance in (4, 8, 16):
        rare = F(1, 2**distance)
        joint = {(x, y, z): (rare if z else 1-rare)*(1+F(1, 2)*x*y*z)/4
                 for x, y, z in itertools.product((-1, 1), (-1, 1), (0, 1))}
        assert sum(joint.values()) == 1 and min(joint.values()) > 0
        covariance = sum(x*y*w for (x, y, z), w in joint.items())
        conditioned = sum(x*y*w for (x, y, z), w in joint.items() if z)/rare
        total_variation = sum(abs(w-(rare if z else 1-rare)/4) for (x, y, z), w in joint.items())/2
        assert covariance == rare/2 and conditioned == F(1, 2)
        assert 2*total_variation == covariance
        rows.append(dict(distance=distance,rare_probability=rare,
                         all_bounded_left_right_covariance_bound=2*total_variation,
                         conditioned_covariance=conditioned))
    return dict(examples=rows,scope='Full-support probability counterexample to an inference; not the YM Hamiltonian')


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


def encode(value):
    """Serialize exact scalars and directed intervals without float conversion."""
    if isinstance(value, Iv):
        return {"lo": str(value.lo), "hi": str(value.hi)}
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def run():
    native.previous.runtime_check()
    theta = F(1, 16000)
    return encode(dict(
        certificate_type='YM65_SPATIAL_CLUSTERING_INVERSE_RESPONSE_LOCALITY_AND_BULK_CONTROL',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        incidence=incidence_controls(),filter=filter_controls(),sources=source_controls(),
        tiling=tiling_controls(),conditioning=conditioning_controls(),parameters=parameters(theta),
        budgets=[response_budget(theta, d, 4, n) for d in (16, 32, 64, 128) for n in (0, 1, 2, 4)],
        susceptibility_row_sum_upper=susceptibility_sum_upper(theta),
        claim_status='FULL_COUPLING_SPATIAL_BOUNDS_IN_EXISTING_LATTICE_WINDOW__CONTINUUM_OPEN',
        model=dict(YM62_Hamiltonian_retained=True,unbounded_on_site_terms_retained=True,
                   volume_uniform_propagation=True,volume_uniform_connected_clustering=True,
                   every_fixed_inverse_response_order_local=True,blocked_resolvent_readout_bounds=True,
                   bulk_plaquette_mean_unique=True,bulk_energy_C1_1=True,
                   uniform_conditional_density_locality=False,Schur_memory_kernel_locality=False,
                   bulk_second_derivative_limit=False,continuum_R4_QFT=False,
                   physical_light_speed_derived=False,physical_mass_prediction=False),
        proof_scope='Written all-content lattice proofs with exact finite and interval controls',
        open='Uniform conditional coarse-vacuum and Schur-kernel locality; nontrivial continuum trajectory'))


def check(cert, result=RESULT, pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM65 fresh certificate/pin mismatch')
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
    print('YM65 PASS', sha)
    rho = parameters(F(1, 16000))['rho']
    print(json.dumps(dict(runtime=sys.version.split()[0],spatial_rate_lower=str(rho.lo),
                          bulk_energy='C1,1',continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
