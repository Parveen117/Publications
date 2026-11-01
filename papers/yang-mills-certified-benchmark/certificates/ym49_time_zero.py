"""YM-49: exact controls for time-zero observables and row-memory closure.

Python 3.12 only. Finite controls do not evaluate infinite-chain row closure.
Default/--check is read-only; --write regenerates this chapter's evidence.
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
from ym6_seam_integer_dock import iv_sqrt
import ym37_space_transfer as y37
import ym46_infinite_volume as y46
import ym48_reflected_time as y48

RESULT=HERE/'YM49_RESULT.json'
PIN=HERE/'EXPECTED_YM49.sha256'
SOURCES=HERE/'YM49_SOURCE_PINS.json'
multiply=y48.multiply
apply=y48.apply
dot=y48.dot
I4=y48.I4
VAC=y48.VAC


def identity(n):
    return [[F(i==j) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(row) for row in zip(*A)]


def add(A,B,scale=F(1)):
    return [[a+scale*b for a,b in zip(row,other)] for row,other in zip(A,B)]


def scale(A,c):
    return [[c*a for a in row] for row in A]


def power(A,n):
    if not isinstance(n,int) or n<0:
        raise ValueError('nonnegative integer power required')
    out=identity(len(A))
    for _ in range(n):
        out=multiply(out,A)
    return out


def inverse(A):
    n=len(A)
    out=[list(row)+list(unit) for row,unit in zip(A,identity(n))]
    for j in range(n):
        pivot=next((i for i in range(j,n) if out[i][j]),None)
        if pivot is None:
            raise ValueError('singular matrix')
        out[j],out[pivot]=out[pivot],out[j]
        v=out[j][j];out[j]=[x/v for x in out[j]]
        for i in range(n):
            if i!=j:
                v=out[i][j];out[i]=[x-v*y for x,y in zip(out[i],out[j])]
    return [row[n:] for row in out]


def psd(A):
    return A==transpose(A) and y37.inertia(A)[1]==0


def diag(values):
    return [[F(a) if i==j else F(0) for j in range(len(values))]
            for i,a in enumerate(values)]


def complex_multiplier(real,imag):
    """Realification of an exact Gaussian-rational diagonal multiplier."""
    R,S=diag(real),diag(imag);n=len(R)
    return [R[i]+[-v for v in S[i]] for i in range(n)]+[
        S[i]+R[i] for i in range(n)]


def jump_gradient_ceiling(values):
    H=y48.generator()
    return max(sum((-H[i][j]*(values[j]-values[i])**2/F(2)
                    for j in range(4) if i!=j),F(0)) for i in range(4))


def observable_controls():
    obs=[(tuple(map(F,y48.SIGNS[1])),(F(0),)*4),
         ((F(1),F(2),F(-1),F(0)),tuple(map(F,y48.SIGNS[2]))),
         ((F(1,2),F(-1,3),F(2,3),F(1,4)),
          (F(2,3),F(0),F(-1,2),F(1)))]
    matrices=[complex_multiplier(*v) for v in obs]
    products=daggers=bounds=0
    for (a,b),M in zip(obs,matrices):
        assert transpose(M)==complex_multiplier(a,tuple(-v for v in b));daggers+=1
        ceiling=max(x*x+y*y for x,y in zip(a,b))
        assert psd(add(scale(identity(8),ceiling),multiply(transpose(M),M),-1))
        bounds+=1
        for (c,d),N in zip(obs,matrices):
            expected=complex_multiplier(tuple(x*z-y*w for x,y,z,w in zip(a,b,c,d)),
                                        tuple(x*w+y*z for x,y,z,w in zip(a,b,c,d)))
            assert multiply(M,N)==expected==multiply(N,M);products+=1
    H=y48.generator();form_bounds=0
    for values in (y48.SIGNS[1],(F(1),F(2),F(0),F(-1)),
                   (F(1,2),F(-1,3),F(2,3),F(1,4))):
        M=diag(values);f=max(abs(x) for x in values);C=jump_gradient_ceiling(values)
        difference=add(add(scale(H,2*f*f),scale(I4,2*C)),
                       multiply(multiply(M,H),M),-1)
        assert psd(difference);form_bounds+=1
    return {'complex_product_checks':products,'complex_dagger_checks':daggers,
            'complex_norm_bounds':bounds,'all_source_form_multiplier_bounds':form_bounds}


def path_controls():
    f,g=y48.SIGNS[1:3]
    histories=[[(F(1),[])],[(F(1),[(1,f)])],
               [(F(1),[(1,f),(2,g)])],
               [(F(2,3),[(2,g)]),(F(-1,5),[(1,f)])]]
    checks=0
    for left,right,values in itertools.product(histories,histories,(f,g,(1,2,-1,0))):
        at_zero=[(F(1),[(0,values)])]
        direct=y48.direct_pair(left,y48.product_history(at_zero,right))
        folded=dot(y48.fold(left),apply(diag(values),y48.fold(right)))
        assert direct==folded;checks+=1
    null=[(F(1),[(2,f)]),(F(-1,2),[(1,f)])]
    assert apply(diag(g),y48.fold(null))==[F(0)]*4
    illegal=y48.product_history(null,[(F(1),[(2,f)])])
    assert y48.direct_pair(illegal,illegal)==F(9,16)
    return {'independent_time_zero_path_checks':checks,
            'null_time_zero_action_is_zero':True,'future_multiplier_null_violation':'9/16'}


def insertion_budget(delta,epsilon,f,C,M=F(1)):
    if epsilon<=0 or not 0<delta<epsilon/2 or min(f,C,M)<0:
        raise ValueError('positive separated collision margin and nonnegative bounds required')
    return 2*M*iv_sqrt(Iv(delta*(f*f/epsilon+C))).hi+2*M*f*delta/epsilon


def collision_controls():
    rows=[];checks=0;previous=None
    values=(F(1),F(2),F(0),F(-1));M=diag(values)
    f=max(abs(x) for x in values);C=jump_gradient_ceiling(values)
    sources=[tuple(F(i==j) for i in range(4)) for j in range(4)]
    sources.append(tuple(F(x,6) for x in (1,-2,3,-4)))
    target=multiply(M,y48.transfer(1))
    for n in (8,32,128,512):
        # delta=log(1+1/n), epsilon=log(2). Exponentials are exact rationals.
        left=y48.spectral_matrix([(F(n,n+1))**e for e in y48.ENERGIES])
        right=y48.spectral_matrix([(F(n+1,2*n))**e for e in y48.ENERGIES])
        assert multiply(left,right)==y48.transfer(1)
        assert F(n+1,n)**2<2
        actual=multiply(multiply(left,M),right)
        error=add(actual,target,-1)
        # delta<=1/n and epsilon>=1/2 give a conservative rational budget.
        budget=insertion_budget(F(1,n),F(1,2),f,C)
        max_error=F(0)
        for w in sources:
            d=apply(error,w);squared=dot(d,d)
            assert squared<=budget*budget*dot(w,w);checks+=1
            max_error=max(max_error,squared)
        if previous is not None:
            old,old_budget=previous;difference=add(actual,old,-1)
            for w in sources:
                d=apply(difference,w)
                assert dot(d,d)<=(budget+old_budget)**2*dot(w,w);checks+=1
        previous=actual,budget
        rows.append({'n':n,'delta_upper':str(F(1,n)),'epsilon_lower':'1/2',
                     'source_error_squared_max':str(max_error),'collision_bound_upper':str(budget)})
    return {'exact_collision_and_cauchy_checks':checks,'samples':rows}


def weighted_form_controls():
    x=y48.COORD
    polynomials=[x[0],y48.add(x[0],y48.scale(x[1],F(1,3))),
                 y48.mul(x[0],x[2])]
    checks=0
    for tilt,f,g in itertools.product((F(-1,3),F(1,3)),polynomials,polynomials):
        h=y48.add(y48.ONE,y48.scale(x[0],tilt));h2=y48.mul(h,h)
        fg=y48.mul(f,g);hfg=y48.mul(h,fg)
        lhs=y48.integral(y48.grad_square(hfg))-y48.integral(
            y48.mul(y48.mul(fg,fg),y48.mul(h,y48.lap(h))))
        rhs=y48.integral(y48.mul(h2,y48.grad_square(fg)))
        assert lhs==rhs>=0
        fbound=sum(abs(c) for c in f.values())
        C=sum((sum(abs(c) for c in y48.deriv(f,a).values())**2
               for a in range(1,4)),F(0))
        qg=y48.integral(y48.mul(h2,y48.grad_square(g)))
        normg=y48.integral(y48.mul(h2,y48.mul(g,g)))
        assert rhs<=2*fbound*fbound*qg+2*C*normg;checks+=1
    return {'nonconstant_weight_product_energy_and_bound_checks':checks}


def row_cut(closed=False):
    if closed:
        return add(y48.PROJECTORS[0],y48.PROJECTORS[1])
    return [[F(1,3) if i<3 and j<3 else F(i==j==3) for j in range(4)]
            for i in range(4)]


def memory_blocks(closed=False):
    U=y48.transfer(1);P=row_cut(closed);Q=add(I4,P,-1)
    return U,P,Q,multiply(multiply(P,U),P),multiply(multiply(P,U),Q),\
        multiply(multiply(Q,U),P),multiply(multiply(Q,U),Q)


def memory_tail_bound(ell_squared,r,N):
    if ell_squared<0 or not 0<=r<1 or not isinstance(N,int) or N<0:
        raise ValueError('nonnegative coupling, strict contraction and integer cut required')
    return ell_squared*r**N/(1-r)


def schur_tail_bound(ell_squared,r,z,N):
    if z<=1:
        raise ValueError('z>1 required for full visible inverse')
    memory_tail_bound(ell_squared,r,N)
    return ell_squared/z*(r/z)**N/(1-r/z)


def closure_controls():
    U,P,Q,A,B,C,D=memory_blocks()
    delta=add(multiply(multiply(P,power(U,2)),P),multiply(A,A),-1)
    target=multiply(transpose(C),C)
    assert delta==target and psd(delta)
    ell_squared=sum(delta[i][i] for i in range(4))
    assert ell_squared==F(7,288)>0
    assert multiply(delta,delta)==scale(delta,ell_squared)
    assert apply(delta,VAC)==[F(0)]*4
    Pc=row_cut(True)
    assert multiply(multiply(add(I4,Pc,-1),U),Pc)==scale(I4,F(0))
    centered=add(P,y48.PROJECTORS[0],-1)
    assert psd(add(scale(centered,F(1,4)),delta,-1))
    # Time-zero observables respect this row cut, while time evolution does not.
    observable=diag((0,0,0,1))
    assert multiply(observable,P)==multiply(P,observable)
    source=apply(observable,VAC)
    cyclic=[list(VAC)]+[apply(power(U,n),source) for n in range(3)]
    gram=[[dot(a,b) for b in cyclic] for a in cyclic]
    assert y37.inertia(gram)==(4,0,0)
    mixed=0
    for s,t in itertools.product((1,2,3),repeat=2):
        Us,Ut=y48.transfer(s),y48.transfer(t)
        As=multiply(multiply(P,Us),P);At=multiply(multiply(P,Ut),P)
        lhs=add(multiply(multiply(P,y48.transfer(s+t)),P),multiply(As,At),-1)
        rhs=multiply(multiply(multiply(multiply(P,Us),Q),Ut),P)
        assert lhs==rhs;mixed+=1
    return {'mixed_time_defect_checks':mixed,'leakage_norm_squared':str(ell_squared),
            'nonclosing_defect_rank':y37.inertia(delta)[0],
            'vacuum_probe_is_blind':True,'closing_control_defect_is_zero':True,
            'observable_time_cyclic_gram_inertia':[4,0,0],
            'actual_chain_row_closure':'OPEN'}


def memory_controls():
    U,P,Q,A,B,C,D=memory_blocks();r=F(1,2);ell_squared=F(7,288)
    centered=add(P,y48.PROJECTORS[0],-1)
    kernels=[];positive=0
    for n in range(12):
        K=multiply(multiply(B,power(D,n)),C);kernels.append(K)
        assert psd(K) and psd(add(scale(centered,ell_squared*r**n),K,-1))
        positive+=1
    seeds=[tuple(map(F,x)) for x in ((1,0,0,0),(1,-2,3,-4),(0,1,-1,0),(1,1,1,1))]
    recursions=tail_checks=0
    for seed in seeds:
        observed=[apply(P,apply(power(U,n),seed)) for n in range(10)]
        y0=apply(Q,seed)
        for n in range(9):
            rhs=apply(A,observed[n])
            initial=apply(B,apply(power(D,n),y0))
            rhs=[a+b for a,b in zip(rhs,initial)]
            assert dot(initial,initial)<=ell_squared*r**(2*n)*dot(y0,y0)
            for j in range(n):
                term=apply(kernels[n-1-j],observed[j])
                rhs=[a+b for a,b in zip(rhs,term)]
            assert rhs==observed[n+1];recursions+=1
            for N in (0,1,3):
                omitted=[F(0)]*4;max_norm=F(0)
                for j in range(n):
                    x=observed[j];mean=dot(VAC,x)
                    xc=[a-mean*b for a,b in zip(x,VAC)]
                    max_norm=max(max_norm,dot(xc,xc))
                    age=n-1-j
                    if age>=N:
                        term=apply(kernels[age],x)
                        omitted=[a+b for a,b in zip(omitted,term)]
                bound=memory_tail_bound(ell_squared,r,N)
                assert dot(omitted,omitted)<=bound*bound*max_norm;tail_checks+=1
    schur_checks=schur_tails=0
    for z in (F(3,2),F(2),F(7,3)):
        full=inverse(add(scale(I4,z),U,-1))
        hidden=inverse(add(scale(I4,z),D,-1))
        sigma=multiply(multiply(B,hidden),C)
        assert psd(sigma)
        assert psd(add(scale(centered,ell_squared/(z-r)),sigma,-1))
        augmented=add(add(add(scale(P,z),A,-1),sigma,-1),Q)
        reduced=multiply(multiply(P,inverse(augmented)),P)
        assert reduced==multiply(multiply(P,full),P);schur_checks+=1
        for N in (0,1,3,8):
            partial=scale(I4,F(0))
            for n in range(N):
                partial=add(partial,scale(kernels[n],z**(-n-1)))
            tail=add(sigma,partial,-1);budget=schur_tail_bound(ell_squared,r,z,N)
            assert psd(tail) and psd(add(scale(centered,budget),tail,-1))
            schur_tails+=1
    return {'positive_memory_coefficient_bounds':positive,'independent_recursion_checks':recursions,
            'memory_tail_checks':tail_checks,'full_inverse_schur_checks':schur_checks,
            'schur_tail_checks':schur_tails,'sampled_gap_ratio':'1/2',
            'leakage_norm_squared':str(ell_squared)}


def noncommuting_controls():
    values=y48.SIGNS[1];M=diag(values);rows=[]
    for n in (1,2,4,8):
        U=y48.transfer(n);comm=add(multiply(U,M),multiply(M,U),-1)
        v=apply(comm,VAC)
        assert dot(v,v)==(1-F(1,2)**n)**2
    for theta,J,R,rho in y46.CELLS:
        gate=y46.joint_gate(theta,J,R,rho)
        assert gate is not None
        gamma=gate['gamma'].lo
        lower=(1-y48.exp_negative(gamma).hi)/2
        assert 0<lower and lower*lower<F(3,8)
        rows.append({'theta_abs':str(theta),'time':'1',
                     'chain_commutator_norm_lower':str(lower),
                     'chain_commutator_norm_squared_upper':'3/8'})
    return {'exact_fixture_noncommutator_checks':4,'inherited_chain_parameter_bounds':rows}


def negative_controls():
    U,P,Q,A,B,C,D=memory_blocks();delta=multiply(B,C)
    hidden=apply(Q,(F(1),F(0),F(0),F(0)))
    visible=apply(P,(F(0),F(0),F(0),F(1)))
    x0=y48.COORD[0]
    true_cross=y48.integral(y48.grad_square(y48.mul(x0,x0)))
    dropped_cross=2*y48.integral(y48.mul(y48.mul(x0,x0),y48.grad_square(x0)))
    M=complex_multiplier((1,2,3,4),(1,0,-1,2))
    z=F(2);full=multiply(multiply(P,inverse(add(scale(I4,z),U,-1))),P)
    naive=multiply(multiply(P,inverse(add(add(scale(P,z),A,-1),Q))),P)
    controls={
        'future_null_multiplier_is_illegal':path_controls()['future_multiplier_null_violation']=='9/16',
        'complex_dagger_is_not_plain_identity':transpose(M)!=M,
        'product_gradient_cross_term_required':true_cross!=dropped_cross,
        'compressed_time_is_not_automatically_a_semigroup':
            multiply(multiply(P,power(U,2)),P)!=multiply(A,A),
        'vacuum_probe_cannot_certify_closure':apply(delta,VAC)==[F(0)]*4 and psd(delta)
            and sum(delta[i][i] for i in range(4))>0,
        'returned_memory_cannot_be_dropped':dot(apply(delta,visible),apply(delta,visible))>0,
        'initial_hidden_source_cannot_be_dropped':dot(apply(B,hidden),apply(B,hidden))>0,
        'finite_memory_cut_is_not_exact_zero':sum(
            multiply(multiply(B,power(D,3)),C)[i][i] for i in range(4))>0,
        'schur_self_energy_cannot_be_dropped':full!=naive,
        'zero_collision_margin_rejected':False,
        'nonstrict_memory_contraction_rejected':False,
        'vacuum_resolvent_pole_rejected':False}
    for key,call in (
        ('zero_collision_margin_rejected',lambda:insertion_budget(F(1,8),F(0),F(1),F(1))),
        ('nonstrict_memory_contraction_rejected',lambda:memory_tail_bound(F(1),F(1),3)),
        ('vacuum_resolvent_pole_rejected',lambda:schur_tail_bound(F(1),F(1,2),F(1),3))):
        try:
            call()
        except ValueError:
            controls[key]=True
    assert all(controls.values())
    return controls


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    return {'certificate_type':'YM49_TIME_ZERO_OBSERVABLE_ACTION_AND_ROW_MEMORY_DIAGNOSTIC',
            'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'observable_algebra':observable_controls(),'time_zero_paths':path_controls(),
            'collision_limits':collision_controls(),'weighted_form':weighted_form_controls(),
            'row_closure_diagnostic':closure_controls(),'bounded_memory':memory_controls(),
            'noncommuting_time_action':noncommuting_controls(),'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_BOUNDED_TIME_ZERO_REPRESENTATION__ORDERED_READOUTS__'
                'POSITIVE_ROW_CLOSURE_CRITERION_AND_GAP_CONTROLLED_MEMORY__'
                'ACTUAL_CHAIN_ROW_CLOSURE_DICT_PHYSICAL_OBSERVABLES_4D_OPEN',
            'evidence_scope':{
                'proof':'written; not mechanically formalized',
                'carrier':'YM48 reflected completion of YM47 histories, admitted full SU(2) heat/reference-functional adapter',
                'window':'abs(theta)<1/1680; kappa(a)=theta*a',
                'observable_algebra':'local coefficient polynomials and their uniform completion; time zero only',
                'coefficient_sector':'R subset K subset H_ref; R=K iff all time defects vanish; K=H_ref not asserted',
                'memory':'exact bounded T(s) blocks; no unproved blocks of the unbounded generator H',
                'fixtures':'independent rational four-state observer cuts and nonconstant S3 test weights, not an infinite-chain leakage calculation',
                'noncommutation':'same-chain coefficient witness; not an interaction, native curvature or physical gauge-invariance identification',
                'open':'actual chain row closure; full bounded-history density; native measure/NCG; physical clock/gauge-observable selection; relativistic/spatial continuum; AF/Clay/QG'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM49 fresh certificate/pin mismatch')
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
    print('YM49 PASS',sha)
    print(json.dumps({'time_zero_path_checks':cert['time_zero_paths']['independent_time_zero_path_checks'],
        'collision_checks':cert['collision_limits']['exact_collision_and_cauchy_checks'],
        'memory_recursions':cert['bounded_memory']['independent_recursion_checks'],
        'refusal_groups':len(cert['negative_controls']),
        'actual_chain_row_closure':'OPEN'},sort_keys=True))


if __name__=='__main__':
    main()
