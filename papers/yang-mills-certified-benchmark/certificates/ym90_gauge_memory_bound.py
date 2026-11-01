"""YM90: gauge-selected dynamic memory for the actual small-theta vacuum.

Written analytic proofs have exact finite controls, not formal verification.
The engineered native parent and matrix fixtures are not Wilson vacua.
Python 3.12 without optimization. Default/--check are read-only.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys

import ym89_gauge_block_contraction as previous
import ym64_vacuum_adapted_tower as adapted

p, algebra = previous.p, previous.algebra
local, poly, native = previous.local, previous.poly, previous.native
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM90_RESULT.json'
PIN = HERE/'EXPECTED_YM90.sha256'
SOURCES = HERE/'YM90_SOURCE_MANIFEST.json'
BASE_COMMIT = '39597244572339fb7ede1e552e2de0789467eab9'
BASE_TREE = 'a0324b4d6f42de5e59a10237b004c85b638afcc6'
PARENT_HASH = '5ce336c592fe53da4bb75e5a1ba2caf11089eaaf4950013a89b28fa5488a511d'
MODEL = dict(
    actual_Wilson_dynamic_relative_memory_bound=True,
    memory_source_is_internal_gauge_invariant=True,
    hidden_conditional_gap_uses_pointwise_geometry_not_assumed_full_gap=True,
    coupling_is_bounded_from_coarse_energy_domain_to_hidden_L2=True,
    volume_and_block_factor_uniform=True,
    coupling_window='0 <= theta <= 1/400',
    block_factor_at_least=2,
    coarse_periodic_side_at_least=4,
    frequency_dependent_resolvent_form_error_proved=True,
    time_kernel_energy_bound_includes_time_zero=True,
    actual_small_theta_dynamic_memory_gate_discharged=True,
    actual_weak_coupling_dynamic_memory_gate_discharged=False,
    arbitrary_fine_pins_preserve_internal_gauge_argument=False,
    coupling_is_bounded_on_coarse_L2=False,
    exact_reduced_dynamics_is_Markovian=False,
    resolvent_order_implies_semigroup_order=False,
    coarse_static_Haar_law_forces_zero_memory=False,
    relative_memory_improves_YM88_direct_gap_constant=False,
    engineered_parent_is_original_Wilson_model=False,
    matrix_fixtures_are_Wilson_vacuum_truncations=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def source_budget(cross_Hessian_norm_squared, conditional_floor):
    a2,d = p.nonnegative(cross_Hessian_norm_squared),p.positive(conditional_floor)
    sigma2 = 4*a2/d
    return dict(cross_Hessian_norm_squared=a2,conditional_floor=d,
                source_energy_factor=sigma2,relative_memory=sigma2/d)


def gauge_charge_floor(max_degree=6,max_internal_endpoints=2):
    p.size(max_degree)
    p.size(max_internal_endpoints)
    return F(3,4)/(max_degree*max_internal_endpoints)


def actual_budget(theta,factor,damping=0):
    fine = previous.previous.parameters(theta)
    coarse = previous.block_budget(theta,factor)
    theta,s = fine['theta'],p.nonnegative(damping)
    # Restrict the pointwise weighted Ricci tensor to unretained coordinates.
    d = F(1,2)-2*(fine['diagonal_log_Hessian_upper']+5*theta)
    out = source_budget((5*theta)**2,d)
    sigma2,eta = out['source_energy_factor'],out['relative_memory']
    c = coarse['compressed_diagonal_gap_lower']
    alpha = sigma2/(d+s)
    assert eta<=F(1,144) and alpha<=eta
    assembled = previous.two_block_floor(c,d,eta)
    return dict(**out,theta=theta,factor=factor,damping=s,
                coarse_diagonal_gap_lower=c,
                memory_frequency_factor=alpha,
                resolvent_relative_form_error=alpha/(1-alpha),
                physical_two_block_gap_lower=assembled,
                charged_sector_gap_lower=gauge_charge_floor(),
                unrestricted_two_block_gap_lower=min(gauge_charge_floor(),assembled),
                existing_YM88_direct_gap_lower=fine['global_bare_gap_lower'])


def remainder_budget(relative_memory,hidden_floor,damping,order,analytic=False):
    eta,d,s = p.nonnegative(relative_memory),p.positive(hidden_floor),p.exact(damping)
    p.size(order,0)
    r = abs(s)/d
    if analytic:
        if r>=1:
            raise ValueError('The analytic tower requires abs(damping) < hidden_floor')
        return eta*r**(order+1)/(1-r)
    if s<0:
        raise ValueError('The signed Stieltjes remainder requires nonnegative damping')
    return eta*r**(order+1)/(1+r)


def matrix_case(C,B,D,hidden_floor,source_factor,damping):
    """B maps coarse to hidden; this is the transpose of YM64's convention."""
    C,D = algebra.square(C),algebra.square(D)
    B = algebra.matrix(B,len(D),len(C))
    d,sigma2,s = p.positive(hidden_floor),p.nonnegative(source_factor),p.nonnegative(damping)
    eta = sigma2/d
    if eta>=1:
        raise ValueError('The resolvent comparison requires relative memory < 1')
    if C!=p.transpose(C) or D!=p.transpose(D):
        raise ValueError('Symmetric diagonal blocks required')
    if not p.psd(C) or not p.psd(p.add(D,algebra.scale(p.identity(len(D)),d),-1)):
        raise ValueError('Positive coarse block and certified hidden floor required')
    Bt = p.transpose(B)
    if not p.psd(p.add(algebra.scale(C,sigma2),p.multiply(Bt,B),-1)):
        raise ValueError('Source energy inequality fails')
    Sigma = adapted.self_energy(Bt,D,s)
    alpha = sigma2/(d+s)
    assert p.psd(Sigma) and p.psd(p.add(algebra.scale(C,alpha),Sigma,-1))
    shifted = p.add(C,algebra.scale(p.identity(len(C)),s))
    R0 = algebra.inverse(shifted)
    R = algebra.inverse(p.add(shifted,Sigma,-1))
    Ru = algebra.inverse(p.add(algebra.scale(C,1-alpha),algebra.scale(p.identity(len(C)),s)))
    assert p.psd(p.add(R,R0,-1)) and p.psd(p.add(Ru,R,-1))
    assert p.psd(p.add(algebra.scale(R0,1/(1-alpha)),Ru,-1))
    return dict(damping=s,memory=Sigma,frequency_factor=alpha,
                compressed_resolvent=R,Markov_resolvent=R0,resolvent_upper=Ru,
                relative_form_error=alpha/(1-alpha))


@lru_cache(maxsize=None)
def matrix_control():
    _,full_C,full_B,D = adapted.fixture()
    C = adapted.sub(full_C,[1,2])
    B = p.transpose(adapted.sub(full_B,[1,2],[0,1,2]))
    d,sigma2 = F(2),F(1,6)
    cases = [matrix_case(C,B,D,d,sigma2,s) for s in (F(0),F(1,3),F(2),F(20))]
    moments = [adapted.moment(p.transpose(B),D,n) for n in range(7)]
    for n,moment in enumerate(moments):
        assert p.psd(p.add(algebra.scale(C,sigma2/d**(n+1)),moment,-1))
    commutator = p.add(p.multiply(moments[0],moments[1]),p.multiply(moments[1],moments[0]),-1)
    assert any(v for row in commutator for v in row)
    return dict(kind='rational_operator_fixture_not_Wilson_vacuum',C=C,B_hidden_by_coarse=B,D=D,
                hidden_floor=d,source_energy_factor=sigma2,relative_memory=sigma2/d,
                frequencies=cases,moments=moments,first_moment_commutator=commutator)


def insufficient_zero_frequency_control():
    C,D,B2,d,eta,s = F(1),F(100),F(25),F(1),F(1,4),F(1)
    zero,shifted,unsupported = B2/D,B2/(D+s),eta*d/(d+s)*C
    assert zero==eta*C and shifted>unsupported and B2>eta*d*C
    return dict(kind='scalar_operator_counterexample',C=C,D=D,B_squared=B2,
                hidden_floor=d,relative_zero_frequency_memory=eta,
                zero_frequency_memory=zero,shifted_memory=shifted,
                unsupported_frequency_upper=unsupported,
                source_energy_hypothesis_satisfied=False)


def sector_control():
    C,D = [[F(1)]],[[F(1),F(0)],[F(0),F(1,16)]]
    B = [[F(1,8)],[F(0)]]
    good = matrix_case(C,[[F(1,8)]],[[F(1)]],1,F(1,64),0)
    good_memory = adapted.self_energy(p.transpose(B),D,0)[0][0]
    bad_B = [[F(0)],[F(1,8)]]
    bad_memory = adapted.self_energy(p.transpose(bad_B),D,0)[0][0]
    assert good_memory==good['memory'][0][0]==F(1,64)
    assert bad_memory==F(1,4)>good_memory
    return dict(kind='gauge_sector_matrix_comparison',D=D,B_invariant=B,B_charged=bad_B,
                invariant_floor=F(1),whole_hidden_floor=F(1,16),
                invariant_source_memory=good_memory,charged_source_memory=bad_memory,
                stronger_floor_requires_source_selection=True)


@lru_cache(maxsize=None)
def network_control(factor):
    geometry = previous.blocking.blocked_geometry(4,factor)
    ends,paths = geometry['fine']['ends'],geometry['paths']
    internal = [ends[e][1] for path in paths for e in path[:-1]]
    assert len(set(internal))==len(internal)
    vertices = set(internal)
    degrees = Counter(v for endpoints in ends for v in endpoints)
    maximum = max(degrees[v] for v in vertices)
    multiplicity = max(sum(v in vertices for v in endpoints) for endpoints in ends)
    assert maximum==6 and multiplicity<=2
    return dict(factor=factor,internal_vertices=len(internal),max_degree=maximum,
                max_internal_endpoints=multiplicity,universal_energy_factor=12,
                universal_charged_gap_lower=gauge_charge_floor())


PARENT_ENDS = ((0,1),(1,2),(0,3),(3,2),(0,2))
PARENT_W_WORD = ((0,1),(1,1),(4,-1))
PARENT_F_WORD = ((0,1),(1,1),(3,-1),(2,-1))


@lru_cache(maxsize=None)
def parent_polynomials():
    return (local.model.loop_polynomial(PARENT_W_WORD,5),
            local.model.loop_polynomial(PARENT_F_WORD,5))


def parent_point(quaternions,tau):
    tau = p.exact(tau)
    qs = tuple(tuple(p.exact(v) for v in q) for q in quaternions)
    if len(qs)!=5 or any(len(q)!=4 or sum(v*v for v in q)!=1 for q in qs):
        raise ValueError('Five unit rational quaternions required')
    W,f = parent_polynomials()
    evaluate = lambda value:local.model.evaluate(value,qs)
    dw = [[evaluate(native.derivative(W,e,k)) for k in (1,2,3)] for e in range(5)]
    df = [[evaluate(native.derivative(f,e,k)) for k in (1,2,3)] for e in range(5)]
    gamma = sum(x*y for row,col in zip(dw,df) for x,y in zip(row,col))
    Q,R = local.qproduct(qs[:2]),local.qproduct(qs[2:4])
    u = local.y.qmul(Q,local.y.qdagger(qs[4]))
    v = local.y.qmul(Q,local.y.qdagger(R))
    leakage = -2*tau*gamma
    assert leakage==-tau*sum(x*y for x,y in zip(u[1:],v[1:]))
    assert sum(x*x for row in dw for x in row)==F(3,4)*(1-u[0]**2)
    S = []
    for path in ((0,1),(2,3)):
        prefix,score = local.y.ONE,[F(0)]*3
        for e in path:
            rotated = previous.adjoint(prefix,dw[e])
            score = [x+tau*y for x,y in zip(score,rotated)]
            prefix = local.y.qmul(prefix,qs[e])
        S.append(score)
    return dict(W=evaluate(W),probe=evaluate(f),source=leakage,horizontal_score=S,
                parent_potential=F(9,4)*tau*u[0]-F(3,4)*tau*tau*(1-u[0]**2))


def positive_parent_series(tau,order,weighted=False):
    """All tails enclosed, with t_n=tau^(2n)/(n!(n+1)!)."""
    tau = p.nonnegative(tau)
    p.size(order,0)
    if tau>F(1,4):
        raise ValueError('The parent control uses 0 <= tau <= 1/4')
    term = F(3,4) if weighted else F(1)
    total = F(0)
    offset = 3 if weighted else 2
    for n in range(order+1):
        total += term
        term *= tau*tau/((n+1)*(n+offset))
    ratio = tau*tau/((order+2)*(order+1+offset))
    return dict(lower=total,upper=total+term/(1-ratio))


def parent_integral_control(tau,order=8):
    tau = p.nonnegative(tau)
    Z = positive_parent_series(tau,order)
    Fw = positive_parent_series(tau,order,weighted=True)
    d = (1-tau)/2
    budget = source_budget(tau*tau/8,d)
    ratio = dict(lower=Fw['lower']/Z['upper'],upper=Fw['upper']/Z['lower'])
    norm = {key:tau*tau*value/4 for key,value in ratio.items()}
    coarse_energy = F(3,4)
    assert norm['upper']<=budget['source_energy_factor']*coarse_energy
    assert budget['relative_memory']<1
    return dict(kind='engineered_gauge_parent_not_original_Wilson_vacuum',tau=tau,
                partial_order=order,normalizer=Z,weighted_normalizer=Fw,
                transverse_moment=ratio,source_norm_squared=norm,
                probe_coarse_energy=coarse_energy,**budget,
                coarse_marginal_is_exact_product_Haar=True,
                source_nonzero=(tau>0))


@lru_cache(maxsize=None)
def parent_native_control():
    W,f = parent_polynomials()
    assert local.model.electric(W,5)==poly.scale(W,F(9,4))
    assert local.model.electric(f,5)==poly.scale(f,F(3))
    for vertex in range(4):
        for axis in (1,2,3):
            assert not local.model.gauss(W,vertex,axis,ends=PARENT_ENDS)
            assert not local.model.gauss(f,vertex,axis,ends=PARENT_ENDS)
    one = local.y.ONE
    qs = (one,(F(3,5),F(4,5),F(0),F(0)),one,one,one)
    reference = parent_point(qs,F(1,10))
    assert reference['source']==-F(8,125)
    gauge = (F(1,2),)*4
    invariance = []
    for vertex in range(4):
        gauges = [one]*4
        gauges[vertex] = gauge
        changed = local.model.gauge_transform(qs,gauges,ends=PARENT_ENDS)
        out = parent_point(changed,F(1,10))
        for name in ('W','probe','source','parent_potential'):
            assert out[name]==reference[name]
        if vertex in (1,3):
            assert out['horizontal_score']==reference['horizontal_score']
        invariance.append(vertex)
    return dict(kind='native_five_link_engineered_parent',reference=reference,
                gauge_vertices_checked=invariance,Gauss_polynomial_identities=24,
                W_electric_eigenvalue=F(9,4),probe_electric_eigenvalue=F(3),
                arbitrary_fine_pins_used=False)


def controls():
    endpoint = previous.previous.THETA_MAX
    return dict(actual_endpoint=actual_budget(endpoint,2),
                actual_blocks=[actual_budget(endpoint,b) for b in (2,3,4,8,16,32)],
                actual_frequencies=[actual_budget(endpoint,3,s) for s in (0,F(1,10),1,10)],
                free_limit=actual_budget(0,2),
                gauge_networks=[network_control(b) for b in (2,3,4)],
                matrix_memory=matrix_control(),source_sector_selection=sector_control(),
                insufficient_zero_frequency=insufficient_zero_frequency_control(),
                parent_native=parent_native_control(),
                parent_integrals=[parent_integral_control(t) for t in (0,F(1,10),F(1,4))],
                actual_Wilson_leakage=adapted.native_controls())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage']!='YM90' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM90 source base mismatch')
    if len(manifest['upstream_sha256'])!=403 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM90 source manifest shape mismatch')
    if previous.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM89 result/pin mismatch')
    checked = {}
    for path,expected in manifest['upstream_sha256'].items():
        actual = sha256((ROOT/path).read_bytes()).hexdigest()
        if actual!=expected:
            raise ValueError('Changed upstream input: '+path)
        checked[path] = actual
    for path in manifest['local_inputs']:
        checked[path] = sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))] = sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2]!=(3,12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return p.encode(dict(certificate_type='YM90_GAUGE_MEMORY_BOUND',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='ACTUAL_SMALL_THETA_DYNAMIC_MEMORY_BOUND__WEAK_COUPLING_AND_CONTINUUM_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM90_GAUGE_MEMORY_BOUND':
        raise ValueError('YM90 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM90 fresh result/pin mismatch')
    return digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args = parser.parse_args()
    cert = run()
    digest = p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(digest+'\n')
    else:
        check(cert)
    print('YM90 PASS',digest)
    print('Actual small-theta dynamic memory bounded; weak coupling and 4D continuum remain open')


if __name__=='__main__':
    main()
