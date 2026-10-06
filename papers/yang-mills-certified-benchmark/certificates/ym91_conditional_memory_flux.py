"""YM91: conditional Poisson and exact-flux upper bounds on dynamic memory.

Native rapidly varying parents are comparison Hamiltonians, not Wilson vacua.
Exact finite controls support the written proofs, not formal verification.
Use Python 3.12 without optimization. Default/--check are read-only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import sys

import ym90_gauge_memory_bound as previous
import ym84_gauge_resolved_scaling as representation

p,a = previous.p,previous.algebra
poly,native,local = previous.poly,previous.native,previous.local
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM91_RESULT.json'
PIN = HERE/'EXPECTED_YM91.sha256'
SOURCES = HERE/'YM91_SOURCE_MANIFEST.json'
BASE_COMMIT = '5e3b853251e5b3d7acee0e8926c5b18f2b69d83d'
BASE_TREE = '8639fa40802cd96d7ac71f5ff03e4bf915e15ef6'
PARENT_HASH = '8b5b659762f13eeb5e93652cf4aadcdc4aecfd4ee0433e48df0680e3667bafbc'
MODEL = dict(
    actual_finite_Wilson_conditional_memory_domination=True,
    conditional_Poisson_tensor_is_an_upper_bound_not_exact_dynamic_memory=True,
    exact_flux_upper_certificate_requires_no_numerical_conditional_gap=True,
    full_physical_gap_still_requires_a_hidden_floor=True,
    divergence_free_corrections_monotonically_improve_flux_Gram=True,
    finite_Poisson_trial_alone_bounds_dynamic_memory_from_above=False,
    lower_conditional_tensor_bounds_dynamic_memory_from_below=False,
    exterior_averaged_tensor_certifies_unrestricted_coarse_energy=False,
    relative_zero_frequency_memory_implies_frequency_decay=False,
    native_comparison_memory_decreases_with_representation_order=True,
    native_comparison_has_arbitrarily_negative_pointwise_weighted_Ricci=True,
    native_comparison_is_original_Wilson_model=False,
    actual_Wilson_small_theta_window_extended=False,
    actual_Wilson_numerical_window='0 <= theta <= 1/400, inherited from YM90',
    actual_weak_coupling_memory_gate_discharged=False,
    continuum_state_survival_proved=False,
    positive_four_dimensional_continuum_mass_gap=False,
    analytic_proofs_formally_machine_verified=False)


def gram(x):
    return p.multiply(p.transpose(x),x)


def memory_comparison(conditional,hidden,scores,damping=0):
    A,D = a.square(conditional),a.square(hidden)
    R = a.matrix(scores,len(A))
    s = p.nonnegative(damping)
    if len(D)!=len(A) or A!=p.transpose(A) or D!=p.transpose(D):
        raise ValueError('Compatible symmetric conditional and hidden operators required')
    if not a.previous.positive_definite(A) or not p.psd(p.add(D,A,-1)):
        raise ValueError('A positive conditional operator and D >= A are required')
    inverse = lambda L:a.inverse(p.add(L,a.scale(p.identity(len(L)),s)))
    tensor = p.multiply(p.transpose(R),p.multiply(inverse(A),R))
    memory = p.multiply(p.transpose(R),p.multiply(inverse(D),R))
    assert p.psd(p.add(tensor,memory,-1))
    return dict(kind='finite_centered_operator_comparison_not_Wilson_vacuum',
                damping=s,conditional_tensor=tensor,dynamic_memory=memory)


def flux_certificate(gradient,scores,flux,cycles=None):
    N = a.matrix(gradient)
    R = a.matrix(scores,len(N[0]))
    J = a.matrix(flux,len(N),len(R[0]))
    if p.multiply(p.transpose(N),J)!=R:
        raise ValueError('The full divergence residual must vanish')
    initial = gram(J)
    if cycles is not None:
        Z = a.matrix(cycles,len(N))
        if any(v for row in p.multiply(p.transpose(N),Z) for v in row):
            raise ValueError('Cycle corrections must be divergence free')
        G = gram(Z)
        if not a.previous.positive_definite(G):
            raise ValueError('Cycle columns must be independent')
        coefficients = p.multiply(a.inverse(G),p.multiply(p.transpose(Z),J))
        J = p.add(J,p.multiply(Z,coefficients),-1)
    upper = gram(J)
    A = gram(N)
    exact = p.multiply(p.transpose(R),p.multiply(a.inverse(A),R))
    assert p.multiply(p.transpose(N),J)==R
    assert p.psd(p.add(upper,exact,-1)) and p.psd(p.add(initial,upper,-1))
    return dict(kind='finite_exact_flux_control',initial_Gram=initial,
                corrected_flux=J,upper_Gram=upper,exact_conditional_tensor=exact)


@lru_cache(maxsize=None)
def matrix_control():
    N = a.matrix([[1,0,0],[0,2,0],[0,0,3],[1,1,0],[0,1,1]])
    A = gram(N)
    extra = a.matrix([[1,1,0],[1,1,0],[0,0,2]])
    D = p.add(A,extra)
    R = a.matrix([[F(1,4),F(1,5)],[F(1,3),0],[0,F(1,6)]])
    Z = a.matrix([[-1,0],[F(-1,2),F(-1,2)],[0,F(-1,3)],[1,0],[0,1]])
    optimum = p.multiply(N,p.multiply(a.inverse(A),R))
    J = p.add(optimum,p.multiply(Z,a.matrix([[F(1,5),F(1,7)],[F(-1,4),F(1,6)]])))
    hierarchy = [flux_certificate(N,R,J),flux_certificate(N,R,J,[[row[0]] for row in Z]),
                 flux_certificate(N,R,J,Z)]
    assert hierarchy[-1]['upper_Gram']==hierarchy[-1]['exact_conditional_tensor']
    assert p.multiply(A,D)!=p.multiply(D,A)
    trials = a.scale(p.multiply(a.inverse(A),R),F(2,3))
    residuals = [a.poisson_control(p.add(A,a.scale(p.identity(3),s)),R,trials,1+s)
                 for s in (F(0),F(1,2),F(4))]
    return dict(gradient=N,conditional_operator=A,hidden_operator=D,scores=R,
                cycles=Z,initial_flux=J,flux_hierarchy=hierarchy,
                frequencies=[memory_comparison(A,D,R,s) for s in (0,F(1,2),4)],
                conditional_residual_enclosures=residuals)


def conditional_lower_is_not_dynamic_lower():
    out = memory_comparison([[1]],[[100]],[[1]])
    trial = a.poisson_control([[1]],[[1]],[[1]],1)
    assert trial['enclosure']['lower'][0][0]>out['dynamic_memory'][0][0]
    return dict(**out,exact_conditional_trial_lower=trial['enclosure']['lower'])


def invisible_hidden_mode(epsilon):
    epsilon = p.positive(epsilon)
    if epsilon>=F(1,2):
        raise ValueError('Use 0 < epsilon < 1/2')
    H = a.matrix([[1,F(1,2),0],[F(1,2),1,0],[0,0,epsilon]])
    assert p.psd(p.add(H,a.scale(p.identity(3),epsilon),-1))
    return dict(kind='unobserved_slow_mode_counterexample',operator=H,
                coarse_floor=F(1),relative_memory=F(1,4),exact_full_gap=epsilon,
                fixed_projected_zero_resolvent=F(4,3))


def averaged_exterior_counterexample(width,height):
    width,height = p.positive(width),p.positive(height)
    if width>=1:
        raise ValueError('Use an exceptional arc of length below one')
    # A periodic piecewise-linear function: slopes +1,-1,0 on three arcs.
    lengths,slopes,tensor = [width/2,width/2,1-width],[1,-1,0],[height,height,0]
    assert sum(w*u for w,u in zip(lengths,slopes))==0
    energy = sum(w*u*u for w,u in zip(lengths,slopes))
    response = sum(w*t*u*u for w,t,u in zip(lengths,tensor,slopes))
    return dict(kind='periodic_H1_gradient_counterexample',arc_lengths=lengths,
                slopes=slopes,tensor_values=tensor,exterior_average=width*height,
                probe_energy=energy,probe_response=response,energy_relative_response=response/energy)


def validate_native(order,epsilon):
    p.size(order)
    if order%2!=1:
        raise ValueError('The symmetric reciprocal series uses an odd representation order')
    epsilon = p.nonnegative(epsilon)
    if epsilon>=1:
        raise ValueError('A strictly positive density requires epsilon < 1')
    return epsilon


def character_coefficients(order):
    p.size(order)
    return [v/(order+1) for v in representation.character(order)]


def derivative(coefficients):
    return [i*v for i,v in enumerate(coefficients)][1:] or [F(0)]


def evaluate(coefficients,x):
    return sum((v*x**k for k,v in enumerate(coefficients)),F(0))


def native_budget(order,epsilon=F(1,4),factor=2):
    epsilon = validate_native(order,epsilon)
    p.size(factor,2)
    casimir = F(order*(order+2),4)
    d = F(3,4)*(1-epsilon)/(1+epsilon)
    q = epsilon**2/((1-epsilon**2)*(order+1)**2)
    eta = factor*q/(1+q)
    k = epsilon*casimir/(6*(1-epsilon))
    cross_lower = factor*k*k
    ricci_at_minus_one = F(1,2)-2*k
    impossible_ratio = 4*cross_lower/F(3,4)**2
    c = F(3*factor,2)
    floor = previous.previous.two_block_floor(c,d,eta) if eta<1 else F(0)
    return dict(kind='engineered_native_parent_not_Wilson_vacuum',order=order,
                epsilon=epsilon,factor=factor,Casimir=casimir,
                normalized_character_L2_squared=F(1,(order+1)**2),
                conditional_floor=d,conditional_gap_upper=F(3,4),
                reciprocal_excess_upper=q,relative_memory_upper=eta,
                physical_coarse_floor=c,physical_assembled_gap_lower=floor,
                pointwise_weighted_Ricci_at_minus_one=ricci_at_minus_one,
                cross_Hessian_norm_squared_lower=cross_lower,
                best_possible_YM90_Hessian_gate_ratio_lower=impossible_ratio,
                YM90_Hessian_gate_impossible=(impossible_ratio>=1),
                strict_memory_gate=(eta<1))


@lru_cache(maxsize=None)
def reciprocal_integral(order,epsilon=F(1,4),degree=8):
    epsilon = validate_native(order,epsilon)
    p.size(degree,0)
    Fpoly = character_coefficients(order)
    square = representation.poly_multiply(Fpoly,Fpoly)
    power,total = [F(1)],F(1)
    moments = [F(1)]
    for k in range(1,degree+2):
        power = representation.poly_multiply(power,square)
        moment = representation.haar_integral(power)
        assert 0<=moment<=F(1,(order+1)**2)
        moments.append(moment)
        if k<=degree:
            total += epsilon**(2*k)*moment
    tail = epsilon**(2*(degree+1))*moments[-1]/(1-epsilon**2)
    upper = total+tail
    delta = dict(lower=1-1/total,upper=1-1/upper)
    return dict(order=order,epsilon=epsilon,partial_degree=degree,
                even_Haar_moments=moments,reciprocal_integral=dict(lower=total,upper=upper),
                flux_defect=delta,certified_flux_constant=1/upper,
                full_tail_upper=tail)


def native_flux_point(order,epsilon,quaternions,constant=1):
    epsilon = validate_native(order,epsilon)
    c = p.exact(constant)
    W,f = previous.parent_polynomials()
    qs = tuple(tuple(p.exact(v) for v in q) for q in quaternions)
    if len(qs)!=5 or any(len(q)!=4 or sum(x*x for x in q)!=1 for q in qs):
        raise ValueError('Five native unit quaternions required')
    value = lambda g:local.model.evaluate(g,qs)
    w = value(W)
    Fp = character_coefficients(order)
    fw,dfw,ddfw = evaluate(Fp,w),evaluate(derivative(Fp),w),evaluate(derivative(derivative(Fp)),w)
    rho = 1+epsilon*fw
    chi_prime = epsilon*dfw/(2*rho)
    chi_second = epsilon*ddfw/(2*rho)-epsilon**2*dfw**2/(2*rho**2)
    Q = local.qproduct(qs[:2])
    flux_checks = []
    for axis in (1,2,3):
        direction = tuple(F(i==axis-1) for i in range(3))
        rotated = previous.previous.adjoint(local.y.qdagger(Q),direction)
        v_w = -sum(rotated[k-1]*value(local.model.right_derivative(W,4,k)) for k in (1,2,3))
        self_w = value(native.derivative(W,0,axis))
        assert v_w==self_w
        v_rho = epsilon*dfw*v_w
        j = -2*(1-c/rho)
        v_j = -2*c*v_rho/rho**2
        negative_divergence = -v_j-j*v_rho/rho
        score = 2*v_rho/rho
        assert negative_divergence==score
        flux_checks.append(dict(score=score,negative_weighted_divergence=negative_divergence))
    gamma = sum(value(native.derivative(W,e,k))*value(native.derivative(f,e,k))
                for e in range(5) for k in (1,2,3))
    source = -2*chi_prime*gamma
    gamma_F = (1-w*w)*dfw*dfw/4
    potential = 3*(epsilon*F(order*(order+2),4)*fw/(2*rho)+epsilon**2*gamma_F/(4*rho*rho))
    J = [[chi_second*value(native.derivative(W,4,k))*value(native.derivative(W,e,l))
          +chi_prime*value(native.derivative(native.derivative(W,e,l),4,k))
          for e in range(4) for l in (1,2,3)] for k in (1,2,3)]
    return dict(W=w,normalized_character=fw,density=rho,probe=value(f),
                source=source,parent_potential=potential,flux_components=flux_checks,
                cross_Hessian=J)


@lru_cache(maxsize=None)
def native_control(order):
    budget = native_budget(order)
    coeffs = character_coefficients(order)
    polynomial = {(k,0,0,0):v for k,v in enumerate(coeffs) if v}
    assert not poly.sphere_reduce(poly.add(poly.lap(polynomial),poly.scale(polynomial,-budget['Casimir'])))
    assert local.y.phi(poly.mul(polynomial,polynomial))==F(1,(order+1)**2)
    # Every imaginary coordinate has mean zero and Rayleigh quotient 3/4.
    density = poly.add(poly.ONE,poly.scale(polynomial,F(1,4)))
    coordinate = poly.COORD[1]
    variance = local.y.phi(poly.mul(density,poly.mul(coordinate,coordinate)))
    energy_poly = poly.scale(poly.add(poly.ONE,poly.scale(poly.mul(coordinate,coordinate),-1)),F(1,4))
    energy = local.y.phi(poly.mul(density,energy_poly))
    assert not local.y.phi(poly.mul(density,coordinate))
    assert variance==F(1,4) and energy/variance==F(3,4)
    one = local.y.ONE
    center = (one,tuple(-v for v in one),one,one,one)
    central = native_flux_point(order,F(1,4),center)
    cross_gram = p.multiply(central['cross_Hessian'],p.transpose(central['cross_Hessian']))
    assert cross_gram==a.scale(p.identity(3),budget['cross_Hessian_norm_squared_lower'])
    qs = (one,(F(3,5),F(4,5),0,0),one,one,one)
    reference = native_flux_point(order,F(1,4),qs,F(7,8))
    for vertex in range(4):
        gauges = [one]*4
        gauges[vertex] = (F(1,2),)*4
        changed = local.model.gauge_transform(qs,gauges,ends=previous.PARENT_ENDS)
        out = native_flux_point(order,F(1,4),changed,F(7,8))
        for name in ('density','probe','source','parent_potential'):
            assert out[name]==reference[name]
    return dict(budget=budget,native_Casimir_identity=True,Rayleigh_variance=variance,
                Rayleigh_energy=energy,central_cross_Hessian_Gram=cross_gram,
                reference=reference,gauge_vertices_checked=4)


def controls():
    endpoint = previous.actual_budget(F(1,400),2)
    return dict(matrix=matrix_control(),
                conditional_lower_counterexample=conditional_lower_is_not_dynamic_lower(),
                invisible_modes=[invisible_hidden_mode(F(1,k)) for k in (4,16,256)],
                averaged_exteriors=[averaged_exterior_counterexample(F(1,k*k),k) for k in (4,16,64)],
                native_budgets=[native_budget(n) for n in (1,3,5,9,31,101)],
                native_controls=[native_control(n) for n in (1,3,5,9)],
                reciprocal_integrals=[reciprocal_integral(n) for n in (1,3,5,9)],
                inherited_actual_Wilson_endpoint=endpoint,
                inherited_frequency_counterexample=previous.insufficient_zero_frequency_control())


def source_checks():
    manifest = json.loads(SOURCES.read_text())
    if manifest['stage']!='YM91' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM91 source base mismatch')
    if len(manifest['upstream_sha256'])!=410 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM91 source manifest shape mismatch')
    if previous.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(previous.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM90 result/pin mismatch')
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
    return p.encode(dict(certificate_type='YM91_CONDITIONAL_MEMORY_FLUX',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='ACTUAL_CONDITIONAL_MEMORY_CRITERION__NATIVE_FLUX_COMPARISON__WILSON_WEAK_COUPLING_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM91_CONDITIONAL_MEMORY_FLUX':
        raise ValueError('YM91 scope mismatch')
    digest = p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM91 fresh result/pin mismatch')
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
    print('YM91 PASS',digest)
    print('Conditional memory and native flux bounds; original Wilson weak coupling and 4D continuum remain open')


if __name__=='__main__':
    main()
