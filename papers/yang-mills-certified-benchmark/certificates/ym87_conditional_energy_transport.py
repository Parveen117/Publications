"""YM87: actual score budgets, conditional energy and sweep certificates.

Python 3.12 exact finite controls accompany the written universal proofs.
Gaussian, Poisson-matrix and Chebyshev controls are comparison models, not
Wilson vacuum computations. Default/--check are read-only; --write updates
only this stage's result and pin.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import isqrt
from pathlib import Path
import sys

import ym86_finite_energy_response as previous

p, poly = previous.p, previous.poly
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM87_RESULT.json'
PIN = HERE/'EXPECTED_YM87.sha256'
SOURCES = HERE/'YM87_SOURCE_MANIFEST.json'
BASE_COMMIT = 'dc4d7e7946e8015117c280ac7688af6fdf9599f0'
BASE_TREE = '1248dfb7eeb222c086af55f9ea1102b42153c67d'
PARENT_HASH = '021dc5a1d78b099d8a7dd0de89c9638f4e3eb82eb24d04ce0120c6af9c8d071d'
MODEL = dict(
    actual_all_coupling_one_link_score_upper_bound=True,
    actual_fixed_block_transfer_coefficients_finite_and_total_volume_uniform=True,
    actual_score_bound_has_small_coupling_restriction=False,
    score_bound_uniform_in_coupling=False,
    actual_spatial_summability_under_full_conditioning_proved=False,
    conditional_energy_dual_score_identity=True,
    finite_poisson_trial_alone_is_an_upper_bound=False,
    full_residual_and_cross_terms_required=True,
    directed_sweep_Lyapunov_gap_criterion_proved=True,
    actual_weak_coupling_sweep_gate_discharged=False,
    projection_preserves_original_local_support=False,
    comparison_families_are_original_Wilson_vacua=False,
    fixed_L2_response_implies_uniform_projection_energy=False,
    Gaussian_example_beats_every_weighted_row_criterion=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def matrix(values,rows=None,cols=None):
    if not values or not values[0]:
        raise ValueError('A nonempty rectangular matrix is required')
    out = [[p.exact(x) for x in row] for row in values]
    width = len(out[0])
    if any(len(row)!=width for row in out) or (rows is not None and len(out)!=rows) or (cols is not None and width!=cols):
        raise ValueError('Matrix dimensions disagree')
    return out


def square(values):
    out = matrix(values)
    if len(out)!=len(out[0]):
        raise ValueError('A square matrix is required')
    return out


def zeros(n):
    return [[F(0) for _ in range(n)] for _ in range(n)]


def scale(a,c):
    return [[c*x for x in row] for row in a]


def outer(a):
    return [[x*y for y in a] for x in a]


def inverse(a):
    a = square(a)
    return p.transpose([previous.solve(a,column) for column in p.identity(len(a))])


def actual_score_budget(theta,source_links=1,source_faces=4,target_links=1):
    theta = p.nonnegative(theta)
    for n in (source_links,source_faces,target_links):
        p.size(n)
    if not theta:
        return dict(theta=theta,wave_log_gradient_squared_upper=F(0),
                    conditional_energy_transfer_squared_upper=F(0),
                    exact_free_conditional_gap=F(3,4),actual_Wilson_model=True)
    density = previous.density_budget(1,4,theta)
    conditional = p.conditional_budget(source_links,2*source_faces*theta)
    d_options = [(name,density[name]['prefactor'],density[name][key])
                 for name,key in (('fixed_time','exp_exponent'),('optimized','exp_exponent_upper'))]
    r_options = [(name,conditional[name]['prefactor'],conditional[name][key])
                 for name,key in (('fixed_time','exponent'),('optimized','exponent_upper'))]
    wave = [dict(method=name,prefactor=32*theta*pref*pref,exp_exponent_upper=2*exponent)
            for name,pref,exponent in d_options]
    transfer = [dict(density_method=dn,gap_method=rn,
                     prefactor=128*target_links*theta*dp*dp/rp,
                     exp_exponent_upper=2*de+re)
                for dn,dp,de in d_options for rn,rp,re in r_options]
    return dict(theta=theta,source_links=source_links,source_incident_faces=source_faces,
                target_links=target_links,one_link_potential_bound=4*theta,
                heat_time=1/(16*theta),wave_log_gradient_squared_upper_options=wave,
                conditional_energy_transfer_squared_upper_options=transfer,
                interpretation='minimum of prefactor times exp(exponent_upper)',
                actual_Wilson_model=True,total_volume_uniform_for_fixed_blocks=True,
                spatial_decay_or_sweep_contraction_certified=False)


def residual_enclosure(energy_gram,cross,residual_gram,rho):
    g = square(energy_gram)
    n = len(g)
    c,r = matrix(cross,n,n),matrix(residual_gram,n,n)
    rho = p.positive(rho)
    if not p.psd(g) or not p.psd(r):
        raise ValueError('Energy and full residual Grams must be positive semidefinite')
    lower = p.add(p.add(g,c),p.transpose(c))
    upper = p.add(lower,r,1/rho)
    return dict(lower=lower,upper=upper,full_residual_required=True,
                Galerkin_orthogonality_assumed=False,
                actual_conditional_integrals_and_floor_must_be_certified=True)


def poisson_control(operator,scores,trials,rho):
    """A finite comparison; L is positive on its declared centered space."""
    operator = square(operator)
    n = len(operator)
    scores = matrix(scores,n)
    trials = matrix(trials,n,len(scores[0]))
    rho = p.positive(rho)
    if not previous.positive_definite(operator) or not p.psd(p.add(operator,p.identity(n),-rho)):
        raise ValueError('The supplied conditional gap floor fails')
    lt = p.multiply(operator,trials)
    residual = p.add(scores,lt,-1)
    g = p.multiply(p.transpose(trials),lt)
    c = p.multiply(p.transpose(trials),residual)
    r = p.multiply(p.transpose(residual),residual)
    enclosure = residual_enclosure(g,c,r,rho)
    exact = p.multiply(p.transpose(scores),p.multiply(inverse(operator),scores))
    assert p.psd(p.add(exact,enclosure['lower'],-1))
    assert p.psd(p.add(enclosure['upper'],exact,-1))
    return dict(kind='finite_Poisson_comparison_not_actual_Wilson_conditional',
                exact_energy_dual_score=exact,energy_gram=g,cross=c,residual_gram=r,
                residual=residual,enclosure=enclosure)


def sweep_data(coefficients,gaps,order=None):
    k = square(coefficients)
    n = len(k)
    rho = [p.positive(x) for x in gaps]
    if len(rho)!=n or any(x<0 for row in k for x in row) or any(k[i][i] for i in range(n)):
        raise ValueError('Nonnegative directed coefficients, zero diagonal and matching positive floors required')
    order = list(range(n)) if order is None else list(order)
    if any(type(i) is not int for i in order) or sorted(order)!=list(range(n)):
        raise ValueError('The sweep order must be a permutation of all blocks')
    prefix = p.identity(n)
    loss,updates,rows = zeros(n),[],[]
    for i in order:
        row = prefix[i][:]
        rows.append(row)
        loss = p.add(loss,outer(row),1/rho[i])
        update = p.identity(n)
        update[i] = [F(0)]*n
        for j in range(n):
            if j!=i:
                update[j][i] = k[j][i]
        updates.append(update)
        prefix = p.multiply(update,prefix)
    assert previous.positive_definite(loss)
    return dict(transfer=prefix,loss=loss,updates=updates,loss_rows=rows,order=order)


def lyapunov_solution(transfer,loss):
    m,d = square(transfer),square(loss)
    n = len(m)
    if len(d)!=n:
        raise ValueError('Lyapunov dimensions disagree')
    pairs = list(product(range(n),repeat=2))
    system = [[F(i==a and j==b)-m[a][i]*m[b][j] for a,b in pairs] for i,j in pairs]
    values = previous.solve(system,[d[i][j] for i,j in pairs])
    out = [values[i*n:(i+1)*n] for i in range(n)]
    assert out==p.transpose(out)
    return out


def gap_certificate(coefficients,gaps,metric,beta,order=None):
    data = sweep_data(coefficients,gaps,order)
    m,d = data['transfer'],data['loss']
    g = matrix(metric,len(m),len(m))
    beta = p.positive(beta)
    if not previous.positive_definite(g):
        raise ValueError('A positive definite Lyapunov metric is required')
    remainder = p.add(p.add(g,p.multiply(p.transpose(m),p.multiply(g,m)),-1),d,-1)
    if not p.psd(remainder):
        raise ValueError('The full sweep Lyapunov inequality fails')
    if not p.psd(p.add(scale(p.identity(len(m)),beta),g,-1)):
        raise ValueError('The claimed metric norm upper bound fails')
    return dict(transfer=m,loss=d,metric=g,Lyapunov_remainder=remainder,
                beta=beta,conditional_bare_gap_lower=1/beta,order=data['order'],
                finite_matrix_gate_passed=True,
                actual_conditional_coefficient_hypotheses_discharged=False,
                continuum_gap_proved=False)


def gaussian_star(leaves=4,coupling=F(3,10)):
    p.size(leaves)
    a = p.positive(coupling)
    root = isqrt(leaves)
    if root*root!=leaves or leaves*a*a>=1:
        raise ValueError('Use a positive definite star with a square number of leaves')
    n = leaves+1
    precision = p.identity(n)
    for j in range(1,n):
        precision[0][j] = precision[j][0] = -a
    covariance = inverse(precision)
    coefficients = [[abs(x) if i!=j else F(0) for j,x in enumerate(row)] for i,row in enumerate(precision)]
    gap = 1-a*root
    out = gap_certificate(coefficients,[1]*n,covariance,1/gap)
    witness = [F(root)]+[F(1)]*leaves
    assert [row[0] for row in p.multiply(precision,[[v] for v in witness])]==[gap*v for v in witness]
    assert out['Lyapunov_remainder']==zeros(n)
    return dict(kind='Gaussian_comparison_not_Wilson_vacuum',leaves=leaves,coupling=a,
                precision=precision,covariance=covariance,
                unweighted_conditional_correlation_row=leaves*a,
                exact_Gaussian_gap=gap,sharp_gradient=witness,certificate=out,
                weighted_row_and_Hessian_criteria_can_also_succeed=True)


def chebyshev(order):
    p.size(order,0)
    if not order:
        return [F(1)]
    old,current = [F(1)],[F(0),F(1)]
    for _ in range(1,order):
        old,current = current,poly.poly_add(poly.poly_multiply([0,2],current),old,-1)
    return current


def gauge_energy_counterexample(order,epsilon=F(1,2)):
    p.size(order,3)
    epsilon = p.positive(epsilon)
    if epsilon>=1:
        raise ValueError('The comparison density must be strictly positive')
    t = chebyshev(order)
    derivative = previous.derivative(t)
    norm = poly.haar_integral(poly.poly_multiply(t,t))
    energy = poly.haar_integral(poly.poly_multiply([1,0,-1],poly.poly_multiply(derivative,derivative)))
    assert t==poly.poly_add(scale_polynomial(poly.character(order),F(1,2)),poly.character(order-2),F(-1,2))
    assert not poly.haar_integral(t) and norm==F(1,2) and energy==F(order*order,2)
    ratio = (1+epsilon)/(1-epsilon)
    return dict(kind='smooth_gauge_invariant_SU2_comparison_not_Wilson_vacuum',
                order=order,epsilon=epsilon,density_ratio=ratio,
                conditional_L2_correlation_squared=epsilon*epsilon/8,
                input_variance=F(1,4),input_full_energy=F(3,4),
                output_variance=epsilon*epsilon/32,
                output_full_energy=epsilon*epsilon*order*order/32,
                energy_amplification=epsilon*epsilon*order*order/24,
                energy_dual_score_norm_at_declared_exterior=epsilon*epsilon*order*order/12,
                full_and_conditional_Poincare_floor=F(3,4)/ratio,
                native_polynomial=t)


def scale_polynomial(values,factor):
    return [factor*x for x in values]


def controls():
    directed = [[0,F(1,5),F(1,10)],[F(1,3),0,F(1,6)],[F(1,8),F(1,7),0]]
    gaps = [F(3,4),F(1),F(2)]
    data = sweep_data(directed,gaps)
    metric = lyapunov_solution(data['transfer'],data['loss'])
    beta = max(sum(abs(x) for x in row) for row in metric)
    return dict(
        actual_score_budgets=[actual_score_budget(t,n,f,j) for t in (F(0),F(1,16000),F(1),F(25))
                             for n,f,j in ((1,4,1),(4,13,4))],
        heat_constant=dict(chosen_M_times_t=F(1,4),squared_prefactor_before_exp=F(9,2),
                           e_upper=F(3),sqrt_e_upper=F(16,9),final_squared_coefficient=F(8)),
        gauge_energy_examples=[gauge_energy_counterexample(n) for n in (3,5,9,17)],
        poisson_residual_examples=[
            poisson_control([[3,0,0],[0,8,0],[0,0,15]],[[1,2],[3,1],[2,-1]],
                            [[F(1,3),F(2,3)],[F(3,8),F(1,8)],[0,0]],3),
            poisson_control([[3]],[[3]],[[F(1,2)]],3),
            poisson_control([[1]],[[1]],[[3]],1)],
        directed_matrix_comparison=gap_certificate(directed,gaps,metric,beta),
        Gaussian_star=gaussian_star(),
        free_full_link_comparison=gap_certificate(zeros(3),[F(3,4)]*3,scale(p.identity(3),F(4,3)),F(4,3)))


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage']!='YM87' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM87 source base mismatch')
    if len(manifest['upstream_sha256'])!=382 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM87 source manifest shape mismatch')
    if previous.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM86 result/pin mismatch')
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
    return p.encode(dict(certificate_type='YM87_CONDITIONAL_ENERGY_TRANSPORT',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='ACTUAL_SCORE_UPPER__CONDITIONAL_ENERGY__GAP_CRITERION__CONTINUUM_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM87_CONDITIONAL_ENERGY_TRANSPORT':
        raise ValueError('YM87 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM87 fresh result/pin mismatch')
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
    print('YM87 PASS',digest)
    print('Actual score upper bounds and conditional energy proved; sweep gap conditional; 4D continuum open')


if __name__=='__main__':
    main()
