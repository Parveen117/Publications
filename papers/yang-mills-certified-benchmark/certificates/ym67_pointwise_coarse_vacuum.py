"""YM67: positive coarse densities, native derivatives and finite exact pins.

Python 3.12 only. Full-space proofs are accompanied by exact polynomial
and directed-interval controls. No interacting vacuum is discretized.
Default/--check are read-only; --write updates this result and hash only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
from math import comb, factorial
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0,str(HERE))
import ym66_local_vacuum_limit as previous

spatial = previous.previous
adapted = spatial.previous
blocking = adapted.previous
native,p,a,y,model,exact,Iv = (previous.native,previous.p,previous.a,previous.y,
                              previous.model,previous.exact,previous.Iv)
RESULT = HERE/'YM67_RESULT.json'
PIN = HERE/'EXPECTED_YM67.sha256'
SOURCES = HERE/'YM67_SOURCE_PINS.json'
DECAY_DIGITS = 120
GRID = 10**DECAY_DIGITS


def up(value):
    value = exact(value)
    return F(-((-value.numerator*GRID)//value.denominator),GRID)


def power_upper(value,order):
    value = spatial.nonnegative(value)
    native.width(order,0)
    result = F(1)
    while order:
        if order&1:
            result = up(result*value)
        order //= 2
        if order:
            value = up(value*value)
    return result


def decay_upper(exponent):
    """Scale to [0,1] before squaring with 120-place directed rounding."""
    exponent = spatial.nonnegative(exponent)
    if not exponent:
        return F(1)
    halvings = 0
    while exponent > 1:
        exponent /= 2
        halvings += 1
    value = min(F(1),spatial.refine.ex(-exponent).hi)
    for _ in range(halvings):
        value = up(value*value)
    return value


def exp_interval(value):
    if not isinstance(value,Iv):
        value = Iv(exact(value))
    if value.lo < 0 or value.hi > 8:
        raise ValueError('Require a positive exponential interval inside [0,8]')
    return p._r(Iv(spatial.refine.ex(value.lo).lo,spatial.refine.ex(value.hi).hi))


def local_constants(M):
    M = spatial.nonnegative(M)
    if M > F(1,2000):
        raise ValueError('Outside the declared one-link interaction budget')
    contraction = p._r(exp_interval(9*M)*Iv(F(1,20)+9*M))
    denominator = Iv(1)-contraction
    if denominator.lo <= 0:
        raise ArithmeticError('The derivative contraction is unresolved')
    B = p._r(Iv(F(9,2)*M*F(21,19)**2)*exp_interval(45*M)/denominator)
    if not M:
        B = Iv(0)
    B2 = p._r(B/Iv(2)+Iv(2)*B*B)
    score = p._r(Iv(2)*B)
    density_second = p._r(B+Iv(6)*B*B)
    log_second = p._r(B+Iv(10)*B*B)
    kappa = exp_interval(p._r(Iv(4)*spatial.pi_interval()*B))
    return dict(interaction_budget=M,contraction=contraction,
                wave_log_gradient=B,wave_second_ratio=B2,score=score,
                density_second_ratio=density_second,log_second=log_second,kappa=kappa)


@lru_cache(maxsize=None,typed=True)
def parameters(eta):
    eta = spatial.nonnegative(eta)
    if eta >= F(1,8000):
        raise ValueError('Outside the strict inherited coupling window')
    return dict(eta=eta,**local_constants(4*eta))


def tail_moments_upper(q,start):
    q = exact(q)
    native.width(start)
    if not 0 <= q < 1:
        raise ValueError('A geometric ratio in [0,1) is required')
    moments = (1/(1-q),)+spatial.geometric_moments(q)[1:]
    power = power_upper(q,start)
    return tuple(up(power*sum(F(comb(k,j))*start**(k-j)*moments[j]
                              for j in range(k+1))) for k in range(5))


def boundary_l1(eta,distance,fine_links):
    parameters(eta)
    native.width(distance)
    native.width(fine_links)
    if not eta:
        return F(0)
    cell = previous.parameters(eta)
    rho = cell['rho']
    q = decay_upper(rho.lo)
    tails = tail_moments_upper(q,distance)
    u = rho.hi
    coefficients = (625,2750*u+750,3300*u+300,1320*u+40,176*u)
    value = 2*fine_links*eta/cell['gamma'].lo*sum(c*t for c,t in zip(coefficients,tails))
    return min(F(2),up(value))


def covariance_upper(eta,distance,smaller_fine_support):
    parameters(eta)
    native.width(distance)
    native.width(smaller_fine_support)
    if not eta:
        return F(0)
    r = previous.parameters(eta)['rho'].lo*distance
    return min(F(1),up((1+F(2*smaller_fine_support,3)+1/(3*r))*decay_upper(r)))


def density_upgrade(eta,coarse_links,error):
    """Explicit dyadic radii realize valid instances of the interpolation bound."""
    cell = parameters(eta)
    native.width(coarse_links)
    error = spatial.nonnegative(error)
    if not error or not eta:
        return dict(bound=F(0),radius=F(1))
    k = power_upper(cell['kappa'].hi,coarse_links)
    lipschitz = 4*coarse_links*cell['wave_log_gradient'].hi*k
    options = [(F(2)*k,F(1))]
    for n in range(81):
        radius = F(1,2**n)
        options.append((up(40**coarse_links*error/radius**(3*coarse_links)+lipschitz*radius),radius))
    value,radius = min(options)
    return dict(bound=value,radius=radius)


def score_upgrade(eta,coarse_links,density_error):
    cell = parameters(eta)
    native.width(coarse_links)
    density_error = spatial.nonnegative(density_error)
    if not density_error or not eta:
        return dict(bound=F(0),step=F(1))
    log_error = power_upper(cell['kappa'].hi,coarse_links)*density_error
    options = [(2*cell['score'].hi,F(1))]
    for n in range(81):
        step = F(1,2**n)
        options.append((up(log_error/step+cell['log_second'].hi*step),step))
    value,step = min(options)
    return dict(bound=value,step=step)


def boundary_upgrade(eta,distance,coarse_links,block_factor):
    native.width(block_factor)
    error = boundary_l1(eta,distance,coarse_links*block_factor)
    density = density_upgrade(eta,coarse_links,error)
    score = score_upgrade(eta,coarse_links,density['bound'])
    return dict(distance=distance,coarse_links=coarse_links,block_factor=block_factor,
                L1_error=error,density=density,score=score)


def pin_budget(eta,distance,target_links,pinned_links,block_factor):
    cell = parameters(eta)
    native.width(target_links)
    native.width(pinned_links)
    native.width(block_factor)
    m = target_links+pinned_links
    covariance = covariance_upper(eta,distance,block_factor*min(target_links,pinned_links))
    density = density_upgrade(eta,m,covariance)
    kY = power_upper(cell['kappa'].hi,target_links)
    kZ = power_upper(cell['kappa'].hi,pinned_links)
    relative_error = up(kY*kZ*density['bound'])
    log_cap = up(8*target_links*spatial.pi_interval().hi*cell['wave_log_gradient'].hi)
    return dict(distance=distance,target_links=target_links,pinned_links=pinned_links,
                block_factor=block_factor,covariance=covariance,joint_density=density,
                conditional_density=min(2*kY,up(kZ*density['bound'])),
                log_ratio=min(log_cap,relative_error),
                score=score_upgrade(eta,m,density['bound']),
                mutual_information=min(log_cap,up(relative_error**2)))


def derivative(poly,edge,direction):
    direction = tuple(exact(x) for x in direction)
    if len(direction) != 3 or sum(x*x for x in direction) != 1:
        raise ValueError('A rational unit direction is required')
    return p.add(*(p.scale(native.derivative(poly,edge,i+1),x)
                   for i,x in enumerate(direction)))


def integrate_link(poly,edge):
    out = {}
    for monomial,coefficient in poly.items():
        local = monomial[4*edge:4*edge+4]
        weight = y.phi({local:F(1)})
        key = monomial[:4*edge]+(0,0,0,0)+monomial[4*edge+4:]
        out[key] = out.get(key,F(0))+coefficient*weight
    return p.clean(out)


def qproduct(factors):
    out = y.ONE
    for factor in factors:
        out = y.qmul(out,factor)
    return out


def word_jet(qs,i,j):
    """D_i D_j of the ordered face, independently from coefficient differentiation."""
    e,u = i
    f,v = j
    U,V = (F(0),)+u,(F(0),)+v
    factors = [q if sign == 1 else y.qdagger(q) for q,(_,sign) in zip(qs,model.ORIENTED_WORD)]
    def first(edge,direction):
        q,sign = qs[edge],model.ORIENTED_WORD[edge][1]
        return tuple(x*F(-1,2) for x in y.qmul(direction,q)) if sign == 1 else tuple(
            x/2 for x in y.qmul(y.qdagger(q),direction))
    if e != f:
        factors[e],factors[f] = first(e,U),first(f,V)
    elif model.ORIENTED_WORD[e][1] == 1:
        factors[e] = tuple(x/4 for x in qproduct((V,U,qs[e])))
    else:
        factors[e] = tuple(x/4 for x in qproduct((y.qdagger(qs[e]),U,V)))
    return qproduct(factors)[0]


def native_controls():
    W = model.loop_polynomial()
    axes = ((F(1),F(0),F(0)),(F(0),F(1),F(0)),(F(3,5),F(4,5),F(0)))
    turns = tuple(y.turns())
    configurations = [tuple(turns[(offset+k)%len(turns)] for k in range(4))
                      for offset in range(4)]
    cases = 0
    for e,f,u,v in itertools.product(range(4),range(4),axes,axes):
        ordered = derivative(derivative(W,f,v),e,u)
        assert not integrate_link(ordered,e)
        for qs in configurations:
            value = model.evaluate(ordered,qs)
            assert value == word_jet(qs,(e,u),(f,v))
            assert abs(value) <= F(1,4)
            cases += 1
    radius = model.source.constant(1,4)
    for e in range(4):
        radius = p.mul(radius,native.lift(y.RADIUS,e,4))
    for e in range(4):
        gradients = [native.derivative(W,e,axis) for axis in (1,2,3)]
        squares = p.add(*(p.mul(d,d) for d in gradients))
        assert not p.add(p.scale(squares,4),p.mul(W,W),p.scale(radius,-1))
        assert all(not integrate_link(d,e) for d in gradients)
    q0 = p.COORD[0]
    commutator = p.add(y.native_generator(y.native_generator(q0,2),1),
                       p.scale(y.native_generator(y.native_generator(q0,1),2),-1))
    assert commutator
    return dict(ordered_word_jet_checks=cases,one_link_gradient_square_identities=4,
                ordered_derivatives_not_symmetric=True,zero_link_averages=True,
                first_derivative_norm_bound=F(1,2),second_derivative_bound=F(1,4))


def pivot_controls():
    turns = tuple(y.turns())
    unit = ((F(0),F(1),F(0),F(0)),(F(0),F(0),F(1),F(0)),(F(0),F(0),F(0),F(1)))
    rows = []
    for length in (1,2,3,5):
        qs = tuple(turns[(i+2)%len(turns)] for i in range(length))
        prefix,coarse = qproduct(qs[:-1]),qproduct(qs)
        recovered = y.qmul(y.qdagger(prefix),coarse)
        assert recovered == qs[-1]
        rotated = [qproduct((y.qdagger(prefix),U,prefix)) for U in unit]
        assert all(R[0] == 0 for R in rotated)
        gram = [[sum(u*v for u,v in zip(U,V)) for V in rotated] for U in rotated]
        assert gram == a.identity(3)
        for U,R in zip(unit,rotated):
            assert qproduct((prefix,R,qs[-1])) == y.qmul(U,coarse)
        for U,V,R,S in zip(unit,unit[1:]+unit[:1],rotated,rotated[1:]+rotated[:1]):
            assert qproduct((prefix,S,R,qs[-1])) == qproduct((V,U,coarse))
        rows.append(dict(length=length,coarse=coarse,pivot=recovered,rotation_gram=gram))
    return dict(cases=rows,pivot_count_per_coarse_link=1,block_length_factor_in_score=1)


def constant_controls():
    endpoint = local_constants(F(1,2000))
    assert endpoint['contraction'].hi < F(11,200)
    assert endpoint['wave_log_gradient'].hi < F(6,2000)
    cases = [parameters(eta) for eta in (F(0),F(1,32000),F(1,16000),F(1,10000))]
    assert cases[0]['wave_log_gradient'].hi == 0 and cases[0]['kappa'].lo == 1
    assert all(left['wave_log_gradient'].hi < right['wave_log_gradient'].lo
               for left,right in zip(cases,cases[1:]))
    for row in cases:
        assert row['contraction'].hi < 1
        assert row['density_second_ratio'].hi <= row['log_second'].hi
    heat = model.blocks.heat_controls()
    assert F(heat['all_content_tail_upper']) < F(1,20)
    return dict(cases=cases,closed_endpoint_contraction=endpoint['contraction'],
                endpoint_B_over_M_upper=endpoint['wave_log_gradient'].hi*2000,
                simple_score_bound='48*eta',inherited_heat_tail=heat['all_content_tail_upper'])


def ball_controls():
    rows = []
    pi = spatial.pi_interval()
    for r in (F(1),F(1,2),F(1,8),F(1,32)):
        upper_sine = r-r**3/6+r**5/120
        lower_sine = upper_sine-r**7/5040
        cap = Iv(r-upper_sine,r-lower_sine)/(Iv(2)*pi)
        assert cap.lo > r**3/40
        rows.append(dict(radius=r,haar_cap=cap,proved_floor=r**3/40))
    assert F(19,240)/pi.hi > F(1,40)
    # Positive density spikes: L1 convergence alone permits persistent pointwise error.
    spikes = [dict(n=n,L1=F(2,n)*(1-F(1,n))**2,sup=1-F(1,n),slope=n)
              for n in (4,16,64)]
    assert all(row['sup'] >= F(3,4) for row in spikes)
    assert spikes[-1]['L1'] < spikes[0]['L1']
    return dict(native_haar_caps=rows,positive_spike_comparison=spikes,
                shrinking_sine_comparison='sup=1/n, first derivative=1, second derivative=n',
                comparison_scope='Continuous density comparisons, not YM vacua')


def conditional_fixture():
    epsilon = F(1,10)
    X = native.lift(p.COORD[0],0,2)
    Y = native.lift(p.COORD[0],1,2)
    w = p.add(model.source.constant(1,2),p.scale(p.mul(X,Y),epsilon))
    dx = native.derivative(w,0,1)
    assert integrate_link(w,0) == model.source.constant(1,2)
    assert integrate_link(w,1) == model.source.constant(1,2)
    probe = (F(0),F(1),F(0),F(0))
    scores = []
    for sign in (1,-1):
        pin = tuple(sign*x for x in y.ONE)
        value = model.evaluate(w,(probe,pin))
        scores.append(model.evaluate(dx,(probe,pin))/value)
    assert scores == [epsilon/2,-epsilon/2]
    covariance = native.phi(p.mul(p.mul(X,Y),w),2)
    assert covariance == epsilon/16
    return dict(scope='Smooth positive native two-link density, not an interacting YM vacuum',
                coupling=epsilon,covariance=covariance,exact_point_scores=scores,
                density_floor=1-epsilon,conditional_mean_at_positive_pin=epsilon/4)


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
    return spatial.encode(dict(
        certificate_type='YM67_POINTWISE_COARSE_VACUUM_NATIVE_SCORES_AND_FINITE_PINS',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        constants=constant_controls(),native=native_controls(),pivots=pivot_controls(),
        balls=ball_controls(),conditional_fixture=conditional_fixture(),
        decay_arithmetic_digits=DECAY_DIGITS,
        marginal_budgets=[boundary_upgrade(eta,R,4,2) for R in (256,512,1024)],
        pin_budgets=[pin_budget(eta,d,4,4,2) for d in (256,512,1024)],
        model=dict(actual_YM62_vacuum_retained=True,unbounded_link_spaces_retained=True,
                   pointwise_positive_coarse_marginals=True,score_bound_independent_of_volume_and_block_length=True,
                   fixed_marginals_C1_limit=True,finite_exact_pin_conditionals=True,
                   finite_pin_spatial_score_decay=True,local_diagonal_drift_energy_bound=True,
                   arbitrary_growing_conditioning_set_uniform_locality=False,
                   full_Schur_memory_kernel_locality=False,time_zero_Markov_identified=False,
                   continuum_R4_QFT=False,physical_mass_prediction=False),
        claim_status='FINITE_POINTWISE_CONDITIONAL_CONTROL__FULL_EXTERIOR_AND_MEMORY_OPEN',
        proof_scope='Written full-space native heat proofs with finite exact and interval checks'))


def check(cert,result=RESULT,pin=PIN):
    sha = p.canonical_sha(cert)
    if json.loads(result.read_text()) != cert or pin.read_text().strip() != sha:
        raise ValueError('YM67 fresh certificate/pin mismatch')
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
    print('YM67 PASS',sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],fixed_coarse_marginals='C1 positive',
                          exact_finite_pins=True,full_memory_locality=False,continuum_R4_QFT=False)))


if __name__ == '__main__':
    main()
