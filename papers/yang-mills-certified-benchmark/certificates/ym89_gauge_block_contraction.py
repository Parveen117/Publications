"""YM89: marginal curvature Schur bounds and actual gauge-product blocking.

Universal written proofs have exact finite controls, not formal verification.
Polynomial, Gaussian and Markov controls are not Wilson vacuum simulations.
Use Python 3.12 without optimization. Default/--check are read-only.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys

import ym88_spatial_score_contraction as previous
import ym63_spatial_blocking_scale as blocking

p, poly, native = previous.p, previous.poly, previous.native
algebra = previous.previous
local = previous.local
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM89_RESULT.json'
PIN = HERE/'EXPECTED_YM89.sha256'
SOURCES = HERE/'YM89_SOURCE_MANIFEST.json'
BASE_COMMIT = 'fcff60dc91f9dacd045f181c37e9e9716a0ff407'
BASE_TREE = '3c8feafa24b66b230c2fc82dbc3bc389e1c423e6'
PARENT_HASH = '688bc714744ed5a1a3dda26e12cda6dff8685b38ce6cd04c8b20c7ceb2cabca7'
MODEL = dict(
    actual_Wilson_marginal_curvature_Schur_bound=True,
    exact_native_coordinate_pins_uniform=True,
    comparison_matrix_elimination_is_associative=True,
    actual_gauge_product_block_score_decays_exponentially=True,
    coarse_conditionals_uniform_over_full_coarse_exterior=True,
    actual_coarse_diagonal_generator_gap_proved=True,
    block_geometry='disjoint straight paths with internal retained degree two',
    coupling_window='0 <= theta <= 1/400',
    coarse_periodic_side_at_least=4,
    block_factor_at_least=2,
    fine_arbitrary_pins_preserve_the_gauge_fiber_argument=False,
    static_Schur_matrix_is_dynamic_Schur_memory=False,
    coarse_diagonal_generator_is_exact_reduced_dynamics=False,
    coarse_static_Haar_limit_is_a_UV_continuum_construction=False,
    actual_weak_coupling_block_gate_discharged=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    finite_comparisons_are_actual_Wilson_vacua=False,
    analytic_proofs_formally_machine_verified=False)


def marginal_budget(theta):
    cell = previous.parameters(theta)
    theta,rho,gamma = cell['theta'],cell['conditional_gap_lower'],cell['global_bare_gap_lower']
    r,rw = 10*theta,20*theta
    return dict(theta=theta,rho=rho,fine_gap=gamma,
                unweighted_J_row=r,weighted_J_row=rw,
                unweighted_hidden_inverse_norm_upper=1/(rho-r),
                weighted_hidden_inverse_norm_upper=1/(rho-rw),
                marginal_mixed_Hessian_row_upper=rho*r/(2*(rho-r)),
                marginal_weighted_mixed_Hessian_row_upper=rho*rw/(2*(rho-rw)),
                marginal_weighted_fill_row_upper=rw*rw/(rho-rw))


def hidden_walk_tail(theta,degree):
    p.size(degree,0)
    cell = marginal_budget(theta)
    rho = cell['rho']
    q = cell['weighted_J_row']/rho
    return q**(degree+1)/(rho*(1-q))


def path_weight(factor):
    p.size(factor,2)
    candidates = [(F((factor-2*radius+1)*2**radius),radius)
                  for radius in range(1,(factor+1)//2+1)]
    weight,radius = max(candidates)
    return dict(weight=weight,radius=radius,guaranteed_far_links=factor-2*radius+1)


def block_budget(theta,factor,targets=1):
    p.size(targets)
    cell = marginal_budget(theta)
    path = path_weight(factor)
    t,tw = cell['marginal_mixed_Hessian_row_upper'],cell['marginal_weighted_mixed_Hessian_row_upper']
    decay = F(1,2**(factor-1))
    row = min(t/(factor-1),tw*(decay+1/path['weight']))
    score = 4*tw*decay
    floor = F(1,2)-2*row
    assert floor >= F(227,480)
    return dict(theta=cell['theta'],factor=factor,targets=targets,
                path_weight=path['weight'],path_radius=path['radius'],
                coarse_log_wave_gradient_upper=score/2,
                coarse_log_density_score_upper=score,
                coarse_diagonal_log_Hessian_upper=tw*decay,
                coarse_off_diagonal_Hessian_row_upper=min(t/factor,tw/path['weight']),
                coarse_full_Hessian_row_upper=row,
                unit_clock_coarse_gap_lower=floor,
                compressed_diagonal_gap_lower=factor*floor,
                conditional_log_oscillation_upper_pi_multiple=2*targets*score,
                conditional_TV_upper_pi_multiple=F(targets,2)*score,
                conditional_TV_rational_upper=min(F(1),F(22,7)*F(targets,2)*score),
                exact_reduced_dynamics_gap_equals_compressed_gap=False)


def submatrix(a,rows,cols):
    return [[a[i][j] for j in cols] for i in rows]


def schur_comparison(values,retained):
    a = algebra.square(values)
    n = len(a)
    if a!=p.transpose(a) or any(a[i][j]>0 for i in range(n) for j in range(n) if i!=j):
        raise ValueError('A symmetric matrix with nonpositive off-diagonal entries is required')
    margin = min(sum(row) for row in a)
    if margin<=0:
        raise ValueError('A strictly positive row margin is required')
    retained = list(retained)
    if not retained or any(type(i) is not int or not 0<=i<n for i in retained) or len(set(retained))!=len(retained):
        raise ValueError('Retained indices must be distinct valid integers')
    hidden = [i for i in range(n) if i not in retained]
    rr = submatrix(a,retained,retained)
    if hidden:
        inverse = algebra.inverse(submatrix(a,hidden,hidden))
        assert all(v>=0 for row in inverse for v in row)
        rz = submatrix(a,retained,hidden)
        fill = p.multiply(rz,p.multiply(inverse,p.transpose(rz)))
    else:
        inverse,fill = [],algebra.zeros(len(retained))
    effective = p.add(rr,fill,-1)
    assert all(v>=0 for row in fill for v in row)
    assert min(sum(row) for row in effective)>=margin
    assert p.psd(p.add(effective,algebra.scale(p.identity(len(retained)),margin),-1))
    return dict(retained=retained,hidden=hidden,hidden_inverse=inverse,
                fill=fill,effective=effective,original_row_margin=margin,
                effective_row_margin=min(sum(row) for row in effective))


def word_polynomial(indices,links):
    value = (local.model.source.constant(1,links),{},{},{})
    for e in indices:
        factor = tuple(native.lift(poly.COORD[k],e,links) for k in range(4))
        value = local.model.quaternion_product(value,factor)
    return value


def fixture(q,r):
    return poly.add(poly.scale(poly.add(poly.mul(q[0],r[0]),poly.mul(q[1],r[2])),F(1,7)),
                    poly.scale(poly.mul(q[0],q[0]),F(1,11)),poly.scale(r[1],F(1,13)))


def adjoint(q,direction):
    return local.y.qmul(local.y.qmul(q,(F(0),)+tuple(direction)),local.y.qdagger(q))[1:]


@lru_cache(maxsize=None)
def native_fiber_control(factor):
    if type(factor) is not int or factor not in (2,3):
        raise ValueError('The exact polynomial fixture uses block factor two or three')
    links = 2*factor
    Q,R = word_polynomial(range(factor),links),word_polynomial(range(factor,links),links)
    pulled = fixture(Q,R)
    coarse = fixture(word_polynomial([0],2),word_polynomial([1],2))
    carriers = [(F(3,5),F(4,5),F(0),F(0)),(F(0),F(5,13),F(12,13),F(0)),
                (F(1,2),F(1,2),F(1,2),F(1,2)),(F(0),F(0),F(0),F(1))]
    qs = [carriers[i%len(carriers)] for i in range(links)]
    coarse_qs = [local.qproduct(qs[:factor]),local.qproduct(qs[factor:])]
    units = [tuple(F(i==j) for i in range(3)) for j in range(3)]
    prefixes = {}
    for path in (range(factor),range(factor,links)):
        prefix = local.y.ONE
        for e in path:
            prefixes[e] = prefix
            prefix = local.y.qmul(prefix,qs[e])
    dirs = {e:[adjoint(local.y.qdagger(prefixes[e]),u) for u in units] for e in prefixes}
    evaluate = lambda value:local.model.evaluate(value,coarse_qs)
    grad = [evaluate(previous.D(coarse,i)) for i in range(3)]
    H = [[evaluate(v) for v in row[:3]] for row in previous.hessian(coarse,2)[:3]]
    raw_same = [[local.model.evaluate(previous.D(previous.D(pulled,3*(factor-1)+j),i),qs)
                 for j in range(3)] for i in range(3)]
    raw_distinct = [[local.model.evaluate(previous.D(previous.D(pulled,3*(links-1)+j),i),qs)
                     for j in range(3)] for i in range(3)]
    same = p.multiply(dirs[0],p.multiply(raw_same,p.transpose(dirs[factor-1])))
    distinct = p.multiply(dirs[0],p.multiply(raw_distinct,p.transpose(dirs[links-1])))
    for i in range(3):
        for j in range(3):
            assert same[i][j]==evaluate(previous.D(previous.D(coarse,i),j))
            assert distinct[i][j]==evaluate(previous.D(previous.D(coarse,3+j),i))
    symmetric = algebra.scale(p.add(same,p.transpose(same)),F(1,2))
    skew = algebra.scale(p.add(same,p.transpose(same),-1),F(1,2))
    assert symmetric==H
    reconstructed = [same[2][1]-same[1][2],same[0][2]-same[2][0],same[1][0]-same[0][1]]
    assert reconstructed==grad and any(grad)
    assert sum(v*v for row in skew for v in row)==sum(v*v for v in grad)/2
    return dict(kind='polynomial_product_pullback_not_Wilson_vacuum',factor=factor,
                same_path_cross_Hessian=same,coarse_covariant_Hessian=H,
                antisymmetric_part=skew,coarse_gradient=grad,
                gradient_reconstructed_from_cross_Hessian=reconstructed,
                distinct_path_cross_Hessian=distinct,
                exact_derivative_components_checked=18)


def gauge_flatten_control(factor):
    p.size(factor,2)
    options = [(F(3,5),F(4,5),F(0),F(0)),(F(1,2),)*4,(F(0),F(0),F(0),F(1))]
    qs = [options[i%3] for i in range(factor)]
    gauges = [local.y.ONE]
    for q in qs[:-1]:
        gauges.append(local.y.qmul(gauges[-1],q))
    gauges.append(local.y.ONE)
    changed = [local.qproduct((gauges[i],qs[i],local.y.qdagger(gauges[i+1]))) for i in range(factor)]
    assert all(q==local.y.ONE for q in changed[:-1])
    assert changed[-1]==local.qproduct(qs)==local.qproduct(changed)
    return dict(factor=factor,holonomy=changed[-1],internal_links_flattened=factor-1,
                endpoint_gauges_fixed=True)


@lru_cache(maxsize=None)
def network_control(factor):
    cell = blocking.blocked_geometry(4,factor)
    paths,edges,ends = cell['paths'],cell['fine']['edges'],cell['fine']['ends']
    retained_degree = Counter(v for path in paths for e in path for v in ends[e])
    interior = [ends[e][1] for path in paths for e in path[:-1]]
    assert all(retained_degree[v]==2 for v in interior)
    rows = [{} for _ in edges]
    _,pairs = local.model.incidence(cell['fine']['faces'],len(edges))
    for (i,j),count in pairs.items():
        rows[i][j] = count
    checks = 0
    for seed in paths[0]:
        distances = [None]*len(edges)
        distances[seed] = 0
        queue = deque([seed])
        while queue:
            i = queue.popleft()
            for j in rows[i]:
                if distances[j] is None:
                    distances[j] = distances[i]+1
                    queue.append(j)
        for path in paths:
            for radius in range(1,(factor+1)//2+1):
                near = sum(distances[e]<radius for e in path)
                assert near<=2*radius-1
                assert sum(2**distances[e] for e in path)>=(factor-2*radius+1)*2**radius
                checks += 1
        if seed==paths[0][0]:
            assert distances[paths[0][-1]]>=factor-1
    return dict(factor=factor,fine_links=len(edges),retained_links=factor*len(paths),
                coarse_links=len(paths),internal_vertices=len(interior),
                every_internal_retained_degree_two=True,ball_path_incidence_checks=checks)


def memory_control(fast,slow):
    fast,slow = p.positive(fast),p.positive(slow)
    # Three states with uniform law and edges 0--1 (fast), 1--2 (slow).
    generator = [[fast,-fast,0],[-fast,fast+slow,-slow],[0,-slow,slow]]
    C,D = F(3,2)*(fast+slow),F(1,2)*(fast+slow)
    B2 = F(3,4)*(fast-slow)**2
    schur = C-B2/D
    eta = B2/(C*D)
    trace,product = 2*(fast+slow),3*fast*slow
    assert schur==6*fast*slow/(fast+slow)
    assert schur*D==product
    # Nonzero eigenvalues solve z^2-trace*z+product=0.
    return dict(kind='reversible_three_state_comparison_not_Wilson_model',
                fast=fast,slow=slow,positive_generator=generator,
                retained_partition=[[0,2],[1]],coarse_stationary_law=[F(2,3),F(1,3)],
                compressed_centered_diagonal=C,hidden_floor=D,coupling_squared=B2,
                relative_memory=eta,zero_frequency_Schur=schur,
                nonzero_eigenvalue_trace=trace,nonzero_eigenvalue_product=product,
                full_gap_lower=product/trace,full_gap_upper=2*product/trace,
                coarse_and_hidden_gaps_alone_are_sufficient=False)


def two_block_floor(coarse,hidden,relative_memory):
    coarse,hidden = p.positive(coarse),p.positive(hidden)
    eta = p.nonnegative(relative_memory)
    if eta>=1:
        raise ValueError('Strictly subunit relative memory is required')
    return coarse*hidden*(1-eta)/(coarse+hidden)


def controls():
    a = [[F(4),F(-1),F(-1,2),F(0)],
         [F(-1),F(5),F(-1),F(-1,3)],
         [F(-1,2),F(-1),F(4),F(-1)],
         [F(0),F(-1,3),F(-1),F(3)]]
    direct = schur_comparison(a,[0,3])
    first = schur_comparison(a,[0,2,3])
    second = schur_comparison(first['effective'],[0,2])
    assert direct['effective']==second['effective']
    return dict(marginal_endpoint=marginal_budget(previous.THETA_MAX),
                hidden_inverse_tails={m:hidden_walk_tail(previous.THETA_MAX,m) for m in (0,1,4,8,16)},
                block_endpoints=[block_budget(previous.THETA_MAX,b) for b in (2,3,4,8,16,32)],
                free_blocks=[block_budget(0,b) for b in (2,4,16)],
                Schur_associativity=dict(direct=direct,sequential=second),
                native_fibers=[native_fiber_control(b) for b in (2,3)],
                gauge_flattening=[gauge_flatten_control(b) for b in (2,3,5)],
                native_networks=[network_control(b) for b in (2,3,4)],
                memory_obstruction=[memory_control(b,F(1,b*b)) for b in (2,4,16,64)])


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage']!='YM89' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM89 source base mismatch')
    if len(manifest['upstream_sha256'])!=396 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM89 source manifest shape mismatch')
    if previous.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM88 result/pin mismatch')
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
    return p.encode(dict(certificate_type='YM89_GAUGE_BLOCK_CONTRACTION',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='ACTUAL_SMALL_THETA_STATIC_BLOCK_CONTRACTION__DYNAMIC_AND_CONTINUUM_GATES_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM89_GAUGE_BLOCK_CONTRACTION':
        raise ValueError('YM89 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM89 fresh result/pin mismatch')
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
    print('YM89 PASS',digest)
    print('Actual small-theta gauge-block contraction; exact dynamic memory and 4D continuum remain open')


if __name__=='__main__':
    main()
