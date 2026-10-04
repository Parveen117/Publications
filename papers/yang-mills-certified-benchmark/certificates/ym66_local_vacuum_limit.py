"""YM66: local vacuum uniqueness, coarse marginals and inverse spectral moments.

Python 3.12 only. Exact finite/interval controls accompany full-space proofs.
No fixture below computes the interacting YM vacuum or its physical mass.
Default/--check are read-only; --write changes this result and hash only.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import itertools
import json
from math import comb
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import ym65_spatial_response_locality as previous

native, p, a, y = previous.native, previous.p, previous.a, previous.y
model, exact, Iv = previous.model, previous.exact, previous.Iv
RESULT = HERE / 'YM66_RESULT.json'
PIN = HERE / 'EXPECTED_YM66.sha256'
SOURCES = HERE / 'YM66_SOURCE_PINS.json'


def positive(value):
    value = exact(value)
    if value <= 0:
        raise ValueError('A positive exact scalar is required')
    return value


def parameters(eta):
    eta = previous.nonnegative(eta)
    return previous.parameters(eta)


def change_budget(eta, delta):
    parameters(eta)
    delta = eta if delta is None else previous.nonnegative(delta)
    if delta > 2*eta:
        raise ValueError('The profile difference exceeds the coupling cube diameter')
    return delta


def tail_moments(q, start, rounded=False):
    """T_k = sum_{d>=start} d^k q^d, k=0,...,4; include the endpoint."""
    q = exact(q)
    native.width(start)
    if not 0 <= q < 1 or not isinstance(rounded, bool):
        raise ValueError('Require an exact ratio in [0,1) and Boolean rounding flag')
    zerobased = (1/(1-q),)+previous.geometric_moments(q)[1:]
    power = previous.iv_power(Iv(q), start).hi if rounded else q**start
    values = tuple(power*sum(F(comb(k, j))*start**(k-j)*zerobased[j]
                             for j in range(k+1)) for k in range(5))
    return tuple(p._r(Iv(value)).hi for value in values) if rounded else values


def boundary_tail(eta, distance, links, delta=None):
    eta = previous.nonnegative(eta)
    cell = parameters(eta)
    native.width(distance)
    native.width(links)
    delta = change_budget(eta, delta)
    if delta == 0:
        return dict(bound=F(0),uncapped_upper=F(0),distance=distance,links=links,delta=delta)
    rho = cell['rho']
    q = previous.exp_negative(rho).hi
    if q >= 1:
        raise ArithmeticError('The interval precision does not resolve the decay')
    tails = tail_moments(q, distance, rounded=True)
    u = rho.hi
    coefficients = (625,2750*u+750,3300*u+300,1320*u+40,176*u)
    value = 2*links*delta/cell['gamma'].lo*sum(c*t for c, t in zip(coefficients, tails))
    upper = p._r(Iv(value)).hi
    return dict(bound=min(F(2),upper),uncapped_upper=upper,
                distance=distance,links=links,delta=delta)


def dynamic_tail(eta, distance, links, time, delta=None):
    eta = previous.nonnegative(eta)
    cell = parameters(eta)
    native.width(distance)
    native.width(links)
    time = abs(exact(time))
    delta = change_budget(eta, delta)
    v = cell['velocity']
    if not delta or not time:
        return F(0)
    if v*time > 8:
        return F(2)  # A valid universal cap outside the finite exponential routine.
    J = (previous.refine.ex(v*time).hi-1)/v-time
    tail = tail_moments(F(1, 2), distance)
    value = 24*links*delta*J*sum(c*t for c, t in zip((125,150,60,8),tail))
    return min(F(2),p._r(Iv(value)).hi)


def conditioning_bound(epsilon, minimum_probability):
    epsilon, minimum_probability = previous.nonnegative(epsilon), positive(minimum_probability)
    if epsilon > 1 or minimum_probability > 1:
        raise ValueError('Probability arguments must not exceed one')
    return min(F(1),2*epsilon/minimum_probability)


def poisson_tail_upper(heat_time, cutoff):
    heat_time, cutoff = previous.nonnegative(heat_time), positive(cutoff)
    return min(F(1),p._r(Iv(2*heat_time/cutoff)/previous.pi_interval()).hi)


def profile_controls():
    eta = F(1,16000)
    rows = []
    for side, periodic in ((4,True),(5,False)):
        lattice, _, faces = previous.graph(side, periodic)
        j = [F(i%5-2,2)*eta for i in range(len(faces))]
        k = [F(0) if i%3 else -value for i, value in enumerate(j)]
        for s in (F(0),F(1,3),F(1)):
            profile = [(1-s)*u+s*v for u, v in zip(j,k)]
            weighted = [Counter() for _ in lattice['edges']]
            for face, value in zip(faces,profile):
                for e in face:
                    for f in face:
                        weighted[e][f] += abs(value)
            assert max(map(abs,profile)) <= eta
            assert max(sum(row.values()) for row in weighted) <= 16*eta
            assert max(row[e] for e,row in enumerate(weighted)) <= 4*eta
            assert max(sum(value for f,value in row.items() if f != e)
                       for e,row in enumerate(weighted)) <= 12*eta
            rows.append(dict(side=side,periodic=periodic,path_parameter=s,
                             potential_norm_budget=sum(map(abs,profile)),
                             weighted_row_max=max(sum(row.values()) for row in weighted)))
    gate = model.block_gate(eta)
    assert gate['alpha'] == F(3,4) and gate['margin'] == F(1,8)
    assert gate['gamma'].lo > 0
    return dict(profiles=rows,alpha=gate['alpha'],margin=gate['margin'],
                independent_coupling_signs_allowed=True,removed_faces_keep_free_links=True)


def geometry_controls():
    rows = []
    for side, periodic in ((8,True),(9,False)):
        lattice, graph, faces = previous.graph(side,periodic)
        anchor = side//2
        index = {edge:i for i,edge in enumerate(lattice['edges'])}
        X = {index[((anchor,anchor,anchor),0)],index[((anchor,anchor,anchor),1)]}
        distances = previous.distances(graph,X)
        shells = Counter(min(distances[e] for e in face) for face in faces)
        assert all(count <= 3*len(X)*(2*d+5)**3 for d,count in shells.items())
        rows.append(dict(side=side,periodic=periodic,support_links=len(X),
                         shell_counts=sorted(shells.items())))
    lattice, graph, faces = previous.graph(9,False)
    inside = {i for i,ends in enumerate(lattice['ends'])
              if all(2 <= coordinate < 7 for vertex in ends for coordinate in vertex)}
    kept = [face for face in faces if face <= inside]
    removed = [face for face in faces if not face <= inside]
    index = {edge:i for i,edge in enumerate(lattice['edges'])}
    X = {index[(x,axis)] for x,axis in (((4,4,4),0),((5,4,4),1),
                                      ((4,5,4),0),((4,4,4),1))}
    distances = previous.distances(graph,X)
    separation = min(min(distances[e] for e in face) for face in removed)
    assert len(kept) == len(model.lattice(5,False)['faces']) == 240
    assert len(inside) == len(model.lattice(5,False)['edges']) == 300
    assert separation >= 1 and len(kept)+len(removed) == len(faces)
    return dict(shells=rows,interior_links=len(inside),interior_faces=len(kept),
                removed_faces=len(removed),distance_to_removed_faces=separation,
                discarded_links_remain_free=True)


def tail_controls():
    cases = []
    for q,R in itertools.product((F(1,4),F(1,2)),(1,4,16)):
        values = tail_moments(q,R)
        next_values = tail_moments(q,R+1)
        distant = tail_moments(q,R+40)
        for k,value in enumerate(values):
            assert value-next_values[k] == R**k*q**R
            direct = sum(F(d**k)*q**d for d in range(R,R+40))
            assert value == direct+distant[k]
        cases.append(dict(q=q,start=R,values=values))
    rho = F(2,7)
    coefficients = (625,2750*rho+750,3300*rho+300,1320*rho+40,176*rho)
    for d in (1,3,16,128):
        assert sum(c*d**k for k,c in enumerate(coefficients)) == (2*d+5)**3*(22*rho*d+5)
    assert tail_moments(F(1,2),1) == (1,2,6,26,150)
    return dict(exact_tail_cases=cases,polynomial_expansion_independent=True,
                endpoint_included=True)


def source_controls():
    # An independently solvable bounded fixture, not a Yang--Mills truncation.
    j,root,gap = F(24,7),F(25,7),F(50,7)
    h,v = [F(4,5),F(3,5)], [F(-3,5),F(4,5)]
    outer = lambda x:[[u*w for w in x] for u in x]
    P,excited = outer(h),outer(v)
    shifted = a.add([[0,-j],[-j,2]],a.diag((1,1)),root-1)
    assert shifted == a.scale(excited,gap)
    Z,B = a.diag((1,-1)),[[0,1],[1,0]]
    centered = lambda C:[u-p.dot(h,p.apply(C,h))*w for u,w in zip(p.apply(C,h),h)]
    response = 2*p.dot(centered(Z),p.apply(a.scale(excited,1/gap),centered(B)))
    assert response == -j/root**3 == F(-1176,15625)
    change = p.dot(h,p.apply(Z,h))-1
    assert change == 1/root-1 == F(-18,25)
    difference = a.add(P,a.diag((1,0)),-1)
    assert difference[0][0]+difference[1][1] == 0
    assert difference[0][0]*difference[1][1]-difference[0][1]*difference[1][0] == F(-9,25)
    return dict(scope='Exact two-level source fixture, not an interacting YM spectrum',
                derivative=response,integrated_readout_change=change,
                quantum_trace_distance=F(6,5),classical_L1_distance=F(18,25),
                off_diagonal_source_response_can_be_negative=True)


def time_controls():
    # H=diag(0,2,5), times pi,pi. All phases in this control are exactly +/-1.
    h = [1,0,0]
    U = a.diag((1,1,-1))
    B0 = [[0,1,0],[1,0,0],[0,0,0]]
    B1 = [[0,0,0],[0,0,1],[0,1,0]]
    B2 = [[0,0,1],[0,0,0],[1,0,0]]
    product = lambda matrices: _matrix_product(matrices)
    tau = lambda B,V: product((V,B,a.transpose(V)))
    word = p.dot(h,p.apply(product((B0,U,B1,U,B2)),h))
    cumulative = p.dot(h,p.apply(product((B0,tau(B1,U),tau(B2,a.multiply(U,U)))),h))
    reset = p.dot(h,p.apply(product((B0,tau(B1,U),tau(B2,U))),h))
    assert word == cumulative == -1 and reset == 1
    # The heat word equals exp(-2u-5v), with u=1/4 and v=1/10.
    first = previous.refine.ex(F(-1,2))
    second = previous.refine.ex(F(-1,2))
    target = previous.refine.ex(-1)
    assert (first*second).lo <= target.hi and target.lo <= (first*second).hi
    # Poisson extension of a negative energy is exp(-u*abs(E)), not exp(-u*E).
    negative_poisson = previous.refine.ex(F(-1,2))
    wrong_heat = previous.refine.ex(F(1,2))
    assert negative_poisson.hi < wrong_heat.lo
    tails = []
    for ratio in (F(1,4),F(1,2),F(3,4)):
        actual = Iv(2)*previous.atan_small(ratio)/previous.pi_interval()
        upper = poisson_tail_upper(ratio,1)
        assert 0 < actual.lo <= actual.hi < upper <= 1
        tails.append(dict(heat_to_cutoff_ratio=ratio,tail=actual,upper=upper))
    return dict(scope='Finite spectral and scalar Poisson controls, not a YM heat simulation',
                real_word=word,incorrect_reset_word=reset,heat_word=target,
                negative_energy_substitution_rejected=True,poisson_tails=tails)


def _matrix_product(matrices):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = a.multiply(result,matrix)
    return result


def total_variation(left,right):
    return sum(abs(left.get(key,0)-right.get(key,0)) for key in set(left)|set(right))/2


def pushforward(weights, transform):
    out = {}
    for key,value in weights.items():
        new = transform(key)
        out[new] = out.get(new,F(0))+value
    return out


def marginal_controls():
    signs = tuple(itertools.product((-1,1),repeat=3))
    epsilon = F(1,3)
    mu = {key:F(1,8) for key in signs}
    nu = {(x,y0,z):(1+epsilon*x*y0*z)/8 for x,y0,z in signs}
    for x,y0,z in signs:
        left = tuple(x*c for c in y.ONE)
        right = tuple(y0*c for c in y.ONE)
        assert y.qmul(left,right) == tuple(x*y0*c for c in y.ONE)
    coarse = lambda key:(key[0]*key[1],key[2])
    coarse_mu,coarse_nu = pushforward(mu,coarse),pushforward(nu,coarse)
    original = total_variation(mu,nu)
    after = total_variation(coarse_mu,coarse_nu)
    forgotten = total_variation(pushforward(coarse_mu,lambda k:k[0]),
                               pushforward(coarse_nu,lambda k:k[0]))
    assert original == after == epsilon/2 and forgotten == 0
    rare_cases = []
    for distance in (4,8,16):
        probability = F(1,2**distance)
        base = {(x,z):(probability if z else 1-probability)/2
                for x,z in itertools.product((-1,1),(0,1))}
        changed = {(x,z):weight*(1+F(1,2)*x*z) for (x,z),weight in base.items()}
        tv = total_variation(base,changed)
        conditional = total_variation({x:base[(x,1)]/probability for x in (-1,1)},
                                      {x:changed[(x,1)]/probability for x in (-1,1)})
        assert tv == probability/4 and conditional == F(1,4)
        assert conditional <= conditioning_bound(tv,probability)
        rare_cases.append(dict(probability=probability,unconditional_tv=tv,conditional_tv=conditional))
    return dict(scope='Native center-sign pushforwards and finite probability comparisons',
                original_tv=original,product_link_tv=after,forgotten_tv=forgotten,
                rare_conditioning=rare_cases)


def spectral_moments(energies,weights,count):
    energies,weights = tuple(map(positive,energies)),tuple(map(previous.nonnegative,weights))
    native.width(count,0)
    if not energies or len(energies) != len(weights) or sum(weights) <= 0:
        raise ValueError('A nonzero positive finite measure is required')
    return tuple(sum(w/E**n for E,w in zip(energies,weights)) for n in range(count+1))


def moment_controls():
    energies,weights = (F(2),F(5)),(F(1,3),F(2,3))
    moments = spectral_moments(energies,weights,12)
    ratios = [moments[n]/moments[n+1] for n in range(12)]
    assert all(x > y > 2 for x,y in zip(ratios,ratios[1:]))
    assert ratios[-1] < F(10001,5000)
    hankel = [[moments[i+j] for j in range(2)] for i in range(2)]
    recurrence = p.apply(a.inverse(hankel),[-moments[2],-moments[3]])
    assert recurrence == [F(1,10),F(-7,10)]
    assert all(moments[n+2]+recurrence[1]*moments[n+1]+recurrence[0]*moments[n] == 0
               for n in range(11))
    recovered_weights = p.apply(a.inverse([[F(1),F(1)],[F(1,2),F(1,5)]]),[moments[0],moments[1]])
    assert tuple(recovered_weights) == weights
    for c in ((1,-1,2),(3,0,-4),(0,1,-1)):
        gram = sum(c[i]*c[j]*moments[i+j] for i in range(3) for j in range(3))
        localizer = sum(c[i]*c[j]*(moments[i+j]/2-moments[i+j+1])
                        for i in range(3) for j in range(3))
        assert gram >= 0 and localizer >= 0
    # A level of energy one with zero probe weight changes the full gap, not this tower.
    invisible = spectral_moments((1,2,5),(0,F(1,3),F(2,3)),12)
    assert invisible == moments
    free = spectral_moments((3,),(F(1,4),),6)
    assert free == tuple(F(1,4)*F(1,3)**n for n in range(7))
    return dict(scope='Two-atom moment reconstruction and native free plaquette; no interacting mass',
                inverse_moments=moments,energy_ratio_upper_bounds=ratios,
                inverse_energy_polynomial=recurrence+[F(1)],
                recovered_energies=energies,recovered_weights=recovered_weights,
                invisible_full_gap=F(1),visible_threshold=F(2),free_plaquette_moments=free)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    checked = {}
    for path,expected in manifest['upstream_sha256'].items():
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
    eta = F(1,16000)
    return previous.encode(dict(
        certificate_type='YM66_BOUNDARY_RESPONSE_LOCAL_VACUUM_LIMIT_AND_SPECTRAL_TOWER',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        profiles=profile_controls(),geometry=geometry_controls(),tails=tail_controls(),
        source=source_controls(),time=time_controls(),marginals=marginal_controls(),
        moments=moment_controls(),parameters=parameters(eta),
        boundary_budgets=[boundary_tail(eta,R,4) for R in (64,128,192,256)],
        dynamics_budgets=[dict(distance=R,time=t,bound=dynamic_tail(eta,R,4,t))
                          for R in (16,32,64) for t in (1,10,100)],
        claim_status='UNIQUE_ADMISSIBLE_LATTICE_VACUUM_AND_HISTORIES__CONTINUUM_OPEN',
        model=dict(YM62_Hamiltonian_retained=True,inhomogeneous_gap_cube=True,
                   local_quantum_state_trace_norm_limit=True,admissible_vacuum_limit_unique=True,
                   heat_history_limit_unique=True,fixed_coarse_marginal_L1_limit=True,
                   fixed_inverse_moments_unique=True,probe_spectral_measure_determined_by_full_tower=True,
                   all_algebraic_ground_states_classified=False,pointwise_conditional_density_control=False,
                   Schur_memory_kernel_locality=False,time_zero_Markov_identified=False,
                   interacting_spectral_measure_computed=False,continuum_R4_QFT=False,
                   physical_mass_prediction=False),
        proof_scope='Written full-space lattice proofs with exact finite and interval controls',
        open='Pointwise conditional vacuum, full spatial Schur kernel, Markov identification, continuum trajectory'))


def check(cert,result=RESULT,pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM66 fresh certificate/pin mismatch')
    return sha


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args = parser.parse_args()
    cert = run()
    sha = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM66 PASS',sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],local_vacuum='unique admissible limit',
                          heat_histories='unique',coarse_marginals='L1 limit',continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
