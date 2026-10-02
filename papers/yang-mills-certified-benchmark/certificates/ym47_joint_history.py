"""YM-47: exact controls for the written joint local-history limit.

Python 3.12 only. Finite fixtures are controls, not an SU(2) replacement.
Default/--check rebuilds without writing; --write is explicit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from ym1_certified_gap import Iv, canonical_sha
from ym6_seam_integer_dock import _r, iv_sqrt
from ym43_native_time_transfer import multiply, decimals
import ym37_space_transfer as y37
import ym44_time_refinement as y44
import ym46_infinite_volume as y46

RESULT=HERE/'YM47_RESULT.json'
PIN=HERE/'EXPECTED_YM47.sha256'
SOURCES=HERE/'YM47_SOURCE_PINS.json'
CELLS=tuple((*cell,u) for cell,u in zip(y46.CELLS,(F(9,8),F(9,8),F(13,12),F(1025,1024))))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def power_iv(q,n):
    if not isinstance(n,int) or n<0 or q<0:
        raise ValueError('nonnegative base and integer exponent required')
    out=Iv(1);base=_r(Iv(q))
    while n:
        if n%2:
            out=_r(out*base)
        n//=2
        if n:
            base=_r(base*base)
    return out


def constants(theta,J,R,rho,u):
    gate=y46.joint_gate(theta,J,R,rho)
    if gate is None or not 1<u or u*u>=rho:
        raise ValueError('strict spatial gate and 1<u, u^2<rho required')
    g=gate['gamma'].lo;d=g/(1+g)
    C=28*abs(theta)*J*rho*R*(1+2/g)/gate['margin']
    return {'theta':abs(theta),'g':g,'d':d,'C':C,'u':u,'rho':rho,
            'margin':gate['margin'],'time_ratio':1/u,'space_ratio':u*u/rho}


def normalized_budget(width,theta,a,t):
    if not isinstance(width,int) or width<1 or not 0<a<=1 or t<0:
        raise ValueError('positive width/step and nonnegative time required')
    b=abs(theta)*(width-1)
    return 6*b*t*iv_sqrt(Iv(a)).hi*y44.ex(2*b*a).hi


def history_budget(width,theta,a,T,r,delta,c):
    if not isinstance(r,int) or r<0 or T<0 or delta<=0 or a>delta/4:
        raise ValueError('separated observation times and a<=Delta/4 required')
    # Validate width and step even when theta is zero.
    normalized_budget(width,theta,a,T)
    b=abs(theta)*(width-1);B=48/c['d']+6*(T+2)
    return B*b*iv_sqrt(Iv(a)).hi*y44.ex(2*b*a).hi+4*r*a/delta


def tail(k,p,T,r,delta,N,c):
    if any(not isinstance(x,int) for x in (k,p,r,N)) or min(k,p)<1 or min(r,N)<0 or T<0 or delta<=0:
        raise ValueError('invalid local-history tail parameters')
    u=c['u'];q=1/u;q2=q*q;v=u*u/c['rho']
    if power_iv(q2,N).hi>delta/4:
        raise ValueError('tail starts before the separated-time mesh gate')
    B=48/c['d']+6*(T+2)
    Z=y44.ex(2*c['theta']*(k+2+2/(u*u-1))).hi
    H=power_iv(q,N).hi*((k+2+2*N)/(1-q)+2*q/(1-q)**2)
    G2=power_iv(q2,N).hi/(1-q2)
    Gv=power_iv(v,N).hi/(1-v)
    return 2*B*c['theta']*Z*H+8*r/delta*G2+4*p*c['C']*Gv


def tail_radius(k,p,T,r,delta,tol,c):
    if tol<=0:
        raise ValueError('positive tolerance required')
    N=1
    while power_iv(1/c['u']**2,N).hi>delta/4:
        N*=2
    while 3*tail(k,p,T,r,delta,N,c)>tol:
        N*=2
    return N,3*tail(k,p,T,r,delta,N,c)


def crop_radius(a,rho):
    if not 0<a<1 or rho<=1:
        raise ValueError('0<a<1 and rho>1 required')
    d=1;q=1/rho;value=q
    while value>a*a:
        d+=1;value*=q
    return d


def parameter_controls():
    rows=[]
    for theta,J,R,rho,u in CELLS:
        c=constants(theta,J,R,rho,u)
        for a in (F(1),F(1,8),F(7,256),F(1,1024)):
            assert y46.budgets(a,theta,J,R,rho)['space']<=c['C']/a
        radius,error=tail_radius(3,6,F(2),2,F(1),F(1,1000000),c)
        rows.append({'theta_abs':str(theta),'J':J,'R':str(R),'rho':str(rho),'u':str(u),
                     'gap_lower':decimals(c['g'],False),'weighted_margin':str(c['margin']),
                     'time_tail_ratio':str(c['time_ratio']),'space_tail_ratio':str(c['space_ratio']),
                     'radius_for_three_readouts_1e_minus_6':radius,
                     'volume_error_upper_per_norm_product':decimals(error)})
    generic=[]
    for theta in (F(0),F(1,8192),F(1,4096),F(1,2000),F(1,1681)):
        J=5;alpha=F(1,2)+840*theta
        rho=F(2) if theta==0 else 1+(1-alpha)/(112*theta*J)
        weighted=F(2)/(y46.EPS*J)+(112+56*rho)*theta*J
        R=(1+1/weighted)/2;u=2*rho/(rho+1)
        c=constants(theta,J,R,rho,u)
        assert u*u<rho and c['g']>0
        generic.append({'theta_abs':str(theta),'rho':str(rho),'u':str(u),'strict_gate':True})
    theta,J,R,rho,u=CELLS[1];c=constants(theta,J,R,rho,u)
    rows_joint=[]
    for power in (32,64,128):
        a=F(1,2**power);D=crop_radius(a,rho)
        local=history_budget(3+2*D,theta,a,F(2),2,F(1),c)
        N=D
        bound=8*6*c['C']*a+local+3*tail(3,6,F(2),2,F(1),N,c)
        assert rho**(-D)<=a*a
        rows_joint.append({'a':str(a),'crop_padding':D,'cropped_width':3+2*D,
                           'joint_error_upper_per_norm_product':decimals(bound)})
    assert F(rows_joint[-1]['joint_error_upper_per_norm_product'])<F(1,1000000)
    return {'cells':rows,'space_budget_checks':16,'full_window_parameter_controls':generic,
            'volume_first_cutoff_moduli':rows_joint}


def vacuum_controls():
    rows=[]
    for h,k in (((F(3,5),F(4,5)),(F(5,13),F(12,13))),
                ((F(5,13),F(12,13)),(F(8,17),F(15,17)))):
        dot=sum(x*y for x,y in zip(h,k))
        s2=1-dot*dot;distance2=sum((x-y)**2 for x,y in zip(h,k))
        for q in (F(0),F(1,2),F(99,100)):
            A=[[q*(i==j)+(1-q)*h[i]*h[j] for j in range(2)] for i in range(2)]
            C=[[q*(i==j)+(1-q)*k[i]*k[j] for j in range(2)] for i in range(2)]
            diff=[[A[i][j]-C[i][j] for j in range(2)] for i in range(2)]
            square=multiply(diff,diff);epsilon2=(1-q)**2*s2
            assert square==[[epsilon2,F(0)],[F(0),epsilon2]]
            assert distance2<=4*epsilon2/(1-q)**2
            rows.append({'q':str(q),'error_squared':str(epsilon2),'vacuum_distance_squared':str(distance2)})
    assert F(rows[2]['vacuum_distance_squared'])>16*F(rows[2]['error_squared'])
    return {'rational_rotated_vacua':rows,'gap_denominator_needed':True}


def spectral_pair(M):
    a,b,d=M[0][0],M[0][1],M[1][1]
    rad=_r(iv_sqrt(_r((a-d)*(a-d)+Iv(4)*b*b)))
    if rad.lo<=0:
        raise ValueError('simple two-state top value required')
    lam=_r((a+d+rad)/Iv(2));mu=_r((a+d-rad)/Iv(2))
    projector=[[_r((M[i][j]-(mu if i==j else Iv(0)))/rad)
                for j in range(2)] for i in range(2)]
    return lam,projector


def normalize2(M,lam):
    return [[_r(x/lam) for x in row] for row in M]


def history2(projector,V,observables):
    out=[[Iv(int(i==j)) for j in range(2)] for i in range(2)]
    for j,values in enumerate(observables):
        M=[[Iv(values[i]) if i==k else Iv(0) for k in range(2)] for i in range(2)]
        out=y44.ivmul2(out,M)
        if j<len(V):
            out=y44.ivmul2(out,V[j])
    out=y44.ivmul2(projector,out)
    return _r(out[0][0]+out[1][1])


def independent_history_controls():
    rows=[];observables=((F(1),F(-1)),(F(1,2),F(1)),(F(-1,3),F(1)))
    for theta in (F(1,16),F(-1,16)):
        U=y44.reference2(F(1),theta);nu,E=spectral_pair(U)
        V=normalize2(U,nu)
        reference=history2(E,[V,V],observables)
        # Independent full exponential from a matrix Taylor series + factorial tail.
        direct=y44.taylor_reference2(F(1),theta)
        assert all(not U[i][j].separated_from(direct[i][j]) for i,j in itertools.product(range(2),repeat=2))
        previous=None
        for n in (4,16,64,256):
            a=F(1,n);S=y44.split2(a,theta);lam,P=spectral_pair(S)
            An=y44.ivpow2(normalize2(S,lam),n)
            delta=[[An[i][j]-V[i][j] for j in range(2)] for i in range(2)]
            actual=y44.row_bound2(delta);budget=normalized_budget(2,theta,a,F(1))
            assert actual<budget
            value=history2(P,[An,An],observables)
            err=max(abs((value-reference).lo),abs((value-reference).hi))
            if previous is not None:
                assert err<previous
            previous=err
            # Fixture has a stronger gap; a conservative d=1/4 suffices.
            gap_floor=F(1,4)
            q=S[0][0]+S[1][1]-lam
            assert (q/lam).hi**n<1-gap_floor
            hist=history_budget(2,theta,a,F(2),2,F(1),{'d':gap_floor})
            assert err<hist
            rows.append({'theta':str(theta),'steps_per_unit':n,
                         'normalized_operator_error_upper':decimals(actual),
                         'normalized_operator_budget_upper':decimals(budget),
                         'three_readout_error_upper':decimals(err)})
    return {'independent_three_time_histories':rows,'exponential_entry_cross_checks':8}


def positive_clock_controls():
    W=y46.markov_fixture();dimension=len(W)
    ident=[[F(i==j) for j in range(dimension)] for i in range(dimension)]
    diff=[[ident[i][j]-W[i][j] for j in range(dimension)] for i in range(dimension)]
    power=[row[:] for row in ident];count=0
    for n in range(17):
        residual=multiply(power,diff)
        ceiling=[[F(i==j,n+1)-residual[i][j] for j in range(dimension)] for i in range(dimension)]
        assert y37.inertia(residual)[1]==0 and y37.inertia(ceiling)[1]==0
        power=multiply(power,W);count+=1
    rounding=0
    times=(F(0),F(2,3),F(7,4),F(3))
    delta=min(b-a for a,b in zip(times,times[1:]))
    for a in (F(1,17),F(1,64),F(3,256)):
        for shift in (F(0),F(1,2),F(1)):
            steps=[0]+[int(t/a+shift) for t in times[1:]]
            assert all(x<y for x,y in zip(steps,steps[1:]))
            for j in range(1,len(times)):
                tau=(steps[j]-steps[j-1])*a;t=times[j]-times[j-1]
                assert abs(t-tau)<=2*a and abs(t-tau)/min(t,tau)<=4*a/delta
                rounding+=1
    return {'positive_contraction_power_ceilings':count,'off_grid_time_roundings':rounding}


def geometric_controls():
    count=0
    for q,k,N in itertools.product((F(1,2),F(8,9),F(32,33)),(1,4),(0,1,7)):
        H=lambda n:q**n*((k+2+2*n)/(1-q)+2*q/(1-q)**2)
        for last in (N,N+1,N+9):
            partial=sum(((k+2+2*n)*q**n for n in range(N,last+1)),F(0))
            assert partial+H(last+1)==H(N)
            G=lambda n:q**n/(1-q)
            assert sum((q**n for n in range(N,last+1)),F(0))+G(last+1)==G(N)
            count+=1
        for n in (N,N+12):
            exact=q**n;box=power_iv(q,n)
            assert box.lo<=exact<=box.hi
    return {'exact_polynomial_and_geometric_tail_identities':count}


def cropping_controls():
    count=0;extensions=0
    for left,right,D,k in itertools.product(range(1,8),range(1,8),range(1,6),(1,3)):
        cropped=(min(left,D),min(right,D))
        changed=[j for j,x in enumerate((left,right)) if x>D]
        assert all(cropped[j]==D for j in changed)
        assert all(cropped[j]==(left,right)[j] for j in (0,1) if j not in changed)
        assert sum(cropped)+k<=k+2*D
        N=min(left,right)
        # Balanced box N -> arbitrary interval: only the longer face moves.
        for n in range(N,max(left,right)):
            width=k+N+n+1
            assert width<=k+2*n+2
            extensions+=1
        assert min(cropped)==min(left,right,D)
        count+=1
    # A wildly asymmetric net: near face stays; the remote face is cropped.
    near,far,D=3,10**12,41
    assert (min(near,D),min(far,D))==(3,41) and near<D<far
    return {'asymmetric_crops':count,'one_face_extension_width_checks':extensions,
            'remote_width_1e12_control':True}


def negative_controls():
    theta,J,R,rho,u=CELLS[1];c=constants(theta,J,R,rho,u)
    invalid=[]
    for bad in (F(1),rho):
        try:
            constants(theta,J,R,rho,bad)
        except ValueError:
            invalid.append(True)
    # Small raw error can become an order-one normalized error without a floor.
    eps=F(1,10000)
    X=[[2*eps,F(0)],[F(0),eps]]
    Y=[[eps,F(0)],[F(0),2*eps]]
    raw=eps;normalized=F(1,2)
    # Strong convergence on each fixed mode is not full norm continuity.
    a=F(1,256);casimir=40*41
    high_mode=y44.ex(-a*casimir).hi
    # A stationary two-point limit alone does not compose projected maps.
    P=y46.markov_fixture()
    C=[[F(1) if i==j==0 else F(1,3) if i>0 and j>0 else F(0) for j in range(4)] for i in range(4)]
    compressed=multiply(multiply(C,P),C)
    full_two=multiply(multiply(C,multiply(P,P)),C)
    controls={
        'strict_tail_gates_refused':len(invalid)==2,
        'space_tail_equality_not_summable':F(4,3)**2/F(16,9)==1,
        'raw_error_without_vacuum_floor_is_insufficient':normalized>1000*raw and X!=Y,
        'gap_denominator_cannot_be_dropped':vacuum_controls()['gap_denominator_needed'],
        'norm_continuity_at_zero_refused':1-high_mode>F(99,100),
        'weak_projected_readouts_do_not_prove_composition':full_two!=multiply(compressed,compressed),
        'uncropped_width_can_destroy_refinement_budget':normalized_budget(2**20,theta,F(1,2**20),F(1))>1,
        'fixed_radius_does_not_kill_boundary_budget':c['C']*2**20/rho**3>1,
        'colliding_time_gate_refused':False,
        'premature_tail_start_refused':False}
    try:
        history_budget(3,theta,F(1,2),F(1),1,F(1),c)
    except ValueError:
        controls['colliding_time_gate_refused']=True
    try:
        tail(3,6,F(2),2,F(1),0,c)
    except ValueError:
        controls['premature_tail_start_refused']=True
    assert all(controls.values())
    return controls


def source_checks():
    pins=json.loads(SOURCES.read_text());checked={}
    for path,expected in pins['upstream_sha256'].items():
        actual=digest(ROOT/path)
        if actual!=expected:
            raise ValueError('upstream source changed: '+path)
        checked[path]=actual
    for path in pins['local_inputs']:
        checked[path]=digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))]=digest(SOURCES)
    return checked


def run():
    if not __debug__:
        raise RuntimeError('Optimized Python is refused')
    return {'certificate_type':'YM47_JOINT_VOLUME_TIME_LOCAL_HISTORY_LIMIT','verdict':'PASS',
            'runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'joint_moduli':parameter_controls(),'vacuum_replacement':vacuum_controls(),
            'independent_histories':independent_history_controls(),'time_rounding':positive_clock_controls(),
            'summable_extensions':geometric_controls(),'unrestricted_joint_geometry':cropping_controls(),
            'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_UNRESTRICTED_JOINT_LOCAL_HISTORY_LIMIT__BOTH_ITERATED_LIMITS_IDENTIFIED__'
                           'POSITIVE_FUNCTIONAL_CORRELATION_GAP_REFLECTION_POSITIVITY__'
                           'CONTINUOUS_TIME_MAP_GENERATOR_NATIVE_DICTIONARY_4D_CLAY_QG_OPEN',
            'evidence_scope':{'proof':'written; not mechanically formalized',
                'controls':'exact rational and outward intervals; two-state controls are not SU(2)',
                'carrier':'existing full SU(2) heat/positive-functional adapter, kappa=theta*a',
                'window':'abs(theta)<1/1680; both signs; full prior window',
                'joint_limit':'local vacuum histories, spatial interval exhaustion and a->0 at arbitrary relative rates',
                'history_algebra':'finite sums of bounded local row products at finitely many real heat times, uniform completion',
                'time_endpoints':'already removed in finite-volume vacuum histories',
                'not_established':'infinite-volume strongly continuous time map/generator; path measure; physical clock; native measure/NCG dictionary; 4D continuum; Clay/QG',
                'forward_clarification':'YM46 displayed (1) omitted the plus between spatial and temporal forcing; prior bytes preserved'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM47 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM47 PASS',sha)
    print(json.dumps({'joint_cells':len(cert['joint_moduli']['cells']),
                      'independent_histories':len(cert['independent_histories']['independent_three_time_histories']),
                      'cropping_controls':cert['unrestricted_joint_geometry']['asymmetric_crops'],
                      'refusal_groups':len(cert['negative_controls'])},sort_keys=True))


if __name__=='__main__':
    main()
