"""YM88: covariant mixed-Hessian decay for the actual small-coupling vacuum.

Exact finite controls accompany the universal proofs in the paper. Polynomial
parent-density fixtures and finite matrices are not Wilson vacuum computations.
Python 3.12 only, without optimization. Default/--check are read-only.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import permutations
import json
from math import comb
from pathlib import Path
import sys

import ym87_conditional_energy_transport as previous
import ym67_pointwise_coarse_vacuum as local
import ym81_local_third_moment_vacuum_response as geometry

p, poly, native = previous.p, local.p, local.native
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM88_RESULT.json'
PIN = HERE/'EXPECTED_YM88.sha256'
SOURCES = HERE/'YM88_SOURCE_MANIFEST.json'
BASE_COMMIT = 'b7e7cf5e9d5e92b2bd6a6071a0e23f427475a800'
BASE_TREE = '213e7b9e30887235ef56c25fa83ade6f6b5b27eb'
PARENT_HASH = 'b272f7faebfb8f3e887777be5b53b25840e7085f6bb8cd0f459343a377eb772f'
THETA_MAX = F(1,400)
MODEL = dict(
    actual_Wilson_vacuum_spatial_Hessian_bound=True,
    actual_full_exterior_conditional_score_summability=True,
    actual_small_coupling_sweep_gate_discharged=True,
    actual_untruncated_finite_lattice_gap_uniform_in_volume=True,
    coupling_window='0 <= theta <= 1/400',
    periodic_spatial_dimension=3,
    minimum_periodic_side=4,
    native_fundamental_Casimir='3/4',
    native_Ricci='1/2',
    global_gap_lower='1/2 - 80 theta',
    full_exterior_conditioning_replaced_by_averaging=False,
    small_quadratic_root_selected_without_continuation=False,
    finite_walk_truncation_alone_is_an_upper_bound=False,
    polynomial_parent_fixtures_are_Wilson_vacua=False,
    analytic_proofs_formally_machine_verified=False,
    actual_weak_coupling_sweep_gate_discharged=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False)


def parameters(theta):
    theta = p.nonnegative(theta)
    if theta > THETA_MAX:
        raise ValueError('Outside the proved window 0 <= theta <= 1/400')
    M = 4*theta
    exp9 = 1/(1-9*M)
    x = 45*M
    exp45 = 1+x+x*x/(2*(1-x/3))
    contraction = exp9*(F(1,20)+9*M)
    B = F(9,2)*M*F(21,19)**2*exp45/(1-contraction)
    assert contraction < 1 and B <= 42*theta
    assert B/2+3*B*B <= 35*theta
    rho = F(1,2)-70*theta
    gap = F(1,2)-80*theta
    c = 1/(1-140*theta)
    radius = 10*theta
    assert c*(6*theta+2*radius*radius) <= radius
    assert gap >= F(3,10)
    return dict(theta=theta,one_link_potential_budget=M,
                exp9M_upper=exp9,exp45M_upper=exp45,
                one_link_derivative_contraction_upper=contraction,
                wave_log_gradient_upper=B,wave_log_gradient_linear_upper=42*theta,
                diagonal_log_Hessian_upper=35*theta,
                weighted_mixed_Hessian_row_upper=radius,
                weighted_source_row=6*theta,conditional_gap_lower=rho,
                weighted_transfer_row_upper=20*theta/rho,
                unweighted_transfer_row_upper=10*theta/rho,
                global_bare_gap_lower=gap,walk_resolvent_factor=c,
                Picard_contraction_upper=4*c*radius,
                Catalan_geometric_ratio=48*c*c*theta)


def endpoint_control():
    exp9,exp45 = F(100,91),F(1067,680)
    assert exp9 < F(11,10) and exp45 < F(8,5)
    c0 = F(11,10)*(F(1,20)+F(9,100))
    B = F(9,2)*F(1,100)*F(21,19)**2*F(8,5)/(1-c0)
    assert c0 == F(77,500)
    assert B == F(1764,16967) and B < F(21,200)
    assert F(21,200)/2+3*F(21,200)**2 == F(3423,40000) < F(7,80)
    out = parameters(THETA_MAX)
    assert out['conditional_gap_lower'] == F(13,40)
    assert out['weighted_transfer_row_upper'] == F(2,13)
    assert out['unweighted_transfer_row_upper'] == F(1,13)
    assert out['Picard_contraction_upper'] == F(2,13)
    assert out['Catalan_geometric_ratio'] == F(48,169)
    return dict(exp9M_rational_upper=exp9,exp45M_rational_upper=exp45,
                simplified_c0_upper=c0,simplified_gradient_upper=B,
                diagonal_Hessian_endpoint_upper=F(3423,40000),
                quadratic_roots=[F(1,40),F(3,10)],
                continuity_from_free_vacuum_required=True,
                old_window_endpoint=F(1,8000),window_endpoint_ratio=20,
                endpoint=out)


def barrier_polynomial(u):
    u = p.nonnegative(u)
    value = F(3,200)-F(13,20)*u+2*u*u
    assert value == 2*(u-F(1,40))*(u-F(3,10))
    return value


def score_tail(theta,distance):
    cell = parameters(theta)
    p.size(distance,1)
    return cell['weighted_transfer_row_upper']/2**distance


def walk_tail(theta,degree):
    cell = parameters(theta)
    p.size(degree,0)
    c,q = cell['walk_resolvent_factor'],cell['Catalan_geometric_ratio']
    return 6*c*cell['theta']*q**degree/(1-q)


def catalan_coefficients(theta,degree):
    cell = parameters(theta)
    p.size(degree,0)
    c = cell['walk_resolvent_factor']
    return [F(comb(2*n,n),n+1)*2**n*c**(2*n+1) for n in range(degree)]


def connection(i,j,k):
    return geometry.epsilon(i%3,j%3,k%3)/2 if i//3==j//3==k//3 else F(0)


def D(value,index):
    return native.derivative(value,index//3,index%3+1)


def variable(index,links):
    p.size(links)
    p.size(index,0)
    if index >= 4*links:
        raise ValueError('Coordinate outside carrier')
    key = [0]*(4*links)
    key[index] = 1
    return {tuple(key):F(1)}


def hessian(value,links):
    gradient = [D(value,i) for i in range(3*links)]
    return [[poly.add(D(gradient[j],i),*(poly.scale(gradient[k],-connection(i,j,k))
                    for k in range(3*links))) for j in range(3*links)] for i in range(3*links)]


def covariant_derivative(tensor,k):
    n = len(tensor)
    return [[poly.add(D(tensor[i][j],k),
                     *(poly.scale(tensor[v][j],-connection(k,i,v)) for v in range(n)),
                     *(poly.scale(tensor[i][v],-connection(k,j,v)) for v in range(n)))
             for j in range(n)] for i in range(n)]


@lru_cache(maxsize=None)
def tensor_control(links=3):
    if type(links) is not int or links not in (2,3):
        raise ValueError('This exact comparison uses two or three links')
    x = [variable(i,links) for i in range(4*links)]
    chi = poly.add(poly.mul(x[0],x[4]),poly.scale(poly.mul(x[1],x[6]),F(2,3)),
                   poly.scale(poly.mul(x[0],x[0]),F(1,7)),
                   poly.scale(poly.mul(x[5],x[5]),F(-1,5)))
    if links==3:
        chi = poly.add(chi,poly.scale(poly.mul(x[5],x[10]),F(3,7)),
                       poly.scale(poly.mul(x[8],x[1]),F(-2,9)),
                       poly.scale(poly.mul(poly.mul(x[0],x[5]),x[10]),F(1,11)))
    n = 3*links
    gradient = [D(chi,i) for i in range(n)]
    H = hessian(chi,links)
    assert all(H[i][j]==H[j][i] for i in range(n) for j in range(n))
    # h=exp(chi) is a positive comparison ground state for V=-Delta chi-|grad chi|^2.
    V = poly.scale(poly.add(*(D(gradient[k],k) for k in range(n)),
                            *(poly.mul(g,g) for g in gradient)),-1)
    rough = [[{} for _ in range(n)] for _ in range(n)]
    drift = [[{} for _ in range(n)] for _ in range(n)]
    for k in range(n):
        dH = covariant_derivative(H,k)
        ddH = covariant_derivative(dH,k)
        for i in range(n):
            for j in range(n):
                rough[i][j] = poly.add(rough[i][j],ddH[i][j])
                drift[i][j] = poly.add(drift[i][j],poly.mul(gradient[k],dH[i][j]))
    checked,missing_Ricci,ordinary_laplacian = 0,0,0
    for i in range(n):
        for j in range(i+1,n):
            if i//3==j//3:
                continue
            rhs = poly.add(D(D(V,j),i),
                           *(poly.scale(poly.mul(H[i][k],H[k][j]),2) for k in range(n)))
            lhs = poly.add(poly.scale(rough[i][j],-1),poly.scale(drift[i][j],-2),H[i][j])
            assert lhs==rhs
            checked += 1
            missing_Ricci += bool(H[i][j])
            raw = poly.add(*(D(D(H[i][j],k),k) for k in range(n)))
            ordinary_laplacian += raw!=rough[i][j]
    assert missing_Ricci==checked and ordinary_laplacian==checked
    return dict(kind='polynomial_parent_density_not_actual_Wilson_vacuum',links=links,
                chi_sha256=p.canonical_sha(p.encode(chi)),
                potential_sha256=p.canonical_sha(p.encode(V)),
                exact_mixed_components_checked=checked,
                residuals_nonzero_if_Ricci_is_dropped=missing_Ricci,
                components_changed_by_ordinary_instead_of_rough_Laplacian=ordinary_laplacian,
                native_Ricci_each_factor=F(1,2),mixed_tensor_damping=F(1))


def lattice_graph(side):
    cell = local.model.lattice(side,True)
    n = len(cell['edges'])
    degree,pairs = local.model.incidence(cell['faces'],n)
    rows = [{} for _ in range(n)]
    for (i,j),multiplicity in pairs.items():
        rows[i][j] = multiplicity
    assert all(d==4 for d in degree)
    assert all(sum(row.values())==12 for row in rows)
    distances = [None]*n
    distances[0] = 0
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in rows[i]:
            if distances[j] is None:
                distances[j] = distances[i]+1
                queue.append(j)
    assert all(d is not None for d in distances)
    return rows,distances


def graph_control(side,degree=5):
    p.size(degree,1)
    rows,distances = lattice_graph(side)
    values = {0:1}
    moments = []
    for step in range(1,degree+1):
        out = {}
        for i,v in values.items():
            for j,multiplicity in rows[i].items():
                out[j] = out.get(j,0)+v*multiplicity
        values = out
        weighted = sum(2**distances[i]*v for i,v in values.items())
        assert sum(values.values())==12**step
        assert weighted <= 24**step
        assert all(distances[i]<=step for i in values)
        moments.append(dict(degree=step,total_walk_multiplicity=sum(values.values()),
                            weighted_walk_sum=weighted,weighted_upper=24**step,
                            diagonal_return_multiplicity=values.get(0,0)))
    return dict(periodic_side=side,links=len(rows),off_diagonal_incidence_row_sum=12,
                maximum_graph_distance=max(distances),walk_controls=moments)


def sweep_metric(coefficients,rho,order=None):
    k = previous.square(coefficients)
    rho = p.positive(rho)
    if k!=p.transpose(k) or any(x<0 for row in k for x in row) or any(k[i][i] for i in range(len(k))):
        raise ValueError('A symmetric nonnegative matrix with zero diagonal is required')
    row = max(sum(v) for v in k)
    if row>=1:
        raise ValueError('The row sum must be strictly below one')
    metric = previous.scale(previous.inverse(p.add(p.identity(len(k)),k,-1)),1/rho)
    out = previous.gap_certificate(k,[rho]*len(k),metric,1/(rho*(1-row)),order)
    assert out['Lyapunov_remainder']==previous.zeros(len(k))
    return out


def controls():
    endpoint = parameters(THETA_MAX)
    k = [[0,F(1,26),F(1,26)],[F(1,26),0,F(1,52)],[F(1,26),F(1,52),0]]
    metrics = [sweep_metric(k,endpoint['conditional_gap_lower'],order)
               for order in permutations(range(3))]
    assert all(v['conditional_bare_gap_lower']==F(3,10) for v in metrics)
    return dict(endpoint=endpoint_control(),native_geometry=geometry.lie_ricci_control(),
                tensor_identities=[tensor_control(2),tensor_control(3)],
                native_lattice=[graph_control(side) for side in (4,5)],
                budgets=[parameters(theta) for theta in (F(0),F(1,8000),F(1,1000),THETA_MAX)],
                endpoint_score_tails={r:score_tail(THETA_MAX,r) for r in (1,2,4,8,16)},
                endpoint_walk_tails={m:walk_tail(THETA_MAX,m) for m in (0,1,4,8,16)},
                finite_matrix_identity_controls=metrics,
                matrix_controls_are_vacuum_computations=False,
                quadratic_large_branch_is_also_admissible=barrier_polynomial(F(3,10))==0,
                quadratic_middle_branch_is_forbidden=barrier_polynomial(F(1,10))<0)


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage']!='YM88' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM88 source base mismatch')
    if len(manifest['upstream_sha256'])!=389 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM88 source manifest shape mismatch')
    if previous.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM87 result/pin mismatch')
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
    return p.encode(dict(certificate_type='YM88_SPATIAL_SCORE_CONTRACTION',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='ACTUAL_SMALL_THETA_SPATIAL_SUMMABILITY_AND_GAP__CONTINUUM_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM88_SPATIAL_SCORE_CONTRACTION':
        raise ValueError('YM88 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM88 fresh result/pin mismatch')
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
    print('YM88 PASS',digest)
    print('Actual spatial score summability and gap for 0 <= theta <= 1/400; 4D continuum open')


if __name__=='__main__':
    main()
