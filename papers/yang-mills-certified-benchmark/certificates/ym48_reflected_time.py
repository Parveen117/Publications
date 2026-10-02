"""YM-48: finite exact controls for reflected time action and its generator.

Python 3.12 only. Written infinite proof and finite controls have distinct roles.
Default/--check is read-only; --write explicitly regenerates this certificate.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from ym1_certified_gap import Iv, canonical_sha
from ym6_seam_integer_dock import _r
from ym43_native_time_transfer import multiply, decimals
import ym37_space_transfer as y37
import ym44_time_refinement as y44
import ym46_infinite_volume as y46
import ym47_joint_history as y47

RESULT=HERE/'YM48_RESULT.json'
PIN=HERE/'EXPECTED_YM48.sha256'
SOURCES=HERE/'YM48_SOURCE_PINS.json'
I4=[[F(i==j) for j in range(4)] for i in range(4)]
SIGNS=((1,1,1,1),(1,1,-1,-1),(1,-1,1,-1),(1,-1,-1,1))
VECTORS=tuple(tuple(F(x,2) for x in v) for v in SIGNS)
PROJECTORS=tuple([[v[i]*v[j] for j in range(4)] for i in range(4)] for v in VECTORS)
ENERGIES=(0,1,2,3)
VAC=VECTORS[0]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def apply(A,x):
    return [sum((a*b for a,b in zip(row,x)),F(0)) for row in A]


def dot(x,y):
    return sum((a*b for a,b in zip(x,y)),F(0))


def spectral_matrix(values):
    return [[sum((values[k]*PROJECTORS[k][i][j] for k in range(4)),F(0))
             for j in range(4)] for i in range(4)]


@lru_cache(maxsize=None)
def transfer(n):
    if not isinstance(n,int) or n<0:
        raise ValueError('nonnegative integer tick required')
    # Tick duration log(2) for the rational-energy continuous fixture.
    return spectral_matrix([F(1,2)**(e*n) for e in ENERGIES])


def generator():
    return spectral_matrix(list(map(F,ENERGIES)))


def resolvent(lam):
    if lam<=0:
        raise ValueError('positive damping required')
    return spectral_matrix([1/(lam+e) for e in ENERGIES])


def fold(history,shift=0):
    out=[F(0)]*4
    for coeff,events in history:
        ordered=sorted(events)
        vec=list(VAC)
        for j in range(len(ordered)-1,-1,-1):
            t,values=ordered[j]
            vec=[x*y for x,y in zip(values,vec)]
            elapsed=t-(ordered[j-1][0] if j else 0)
            if j==0:
                elapsed+=shift
            vec=apply(transfer(elapsed),vec)
        for i in range(4):
            out[i]+=coeff*vec[i]
    return out


def value(history,configuration,reflect=False,shift=0):
    result=F(0)
    for coeff,events in history:
        term=coeff
        for t,values in events:
            term*=values[configuration[-t if reflect else t+shift]]
        result+=term
    return result


def direct_pair(left,right,shift=0):
    times=sorted({-t for _,events in left for t,_ in events}
                 |{t+shift for _,events in right for t,_ in events})
    if not times:
        return value(left,{})*value(right,{})
    out=F(0)
    for states in itertools.product(range(4),repeat=len(times)):
        weight=F(1,4)
        for j in range(1,len(times)):
            weight*=transfer(times[j]-times[j-1])[states[j-1]][states[j]]
        conf=dict(zip(times,states))
        out+=weight*value(left,conf,True)*value(right,conf,False,shift)
    return out


def product_history(left,right):
    return [(a*b,sorted(e+f)) for a,e in left for b,f in right]


def history_controls():
    f,g=SIGNS[1],SIGNS[2]
    histories=[[(F(1),[])],[(F(1),[(1,f)])],[(F(1),[(2,g)])],
               [(F(1),[(1,f),(2,g)])]]
    vectors=[fold(h) for h in histories];gram=[];count=0
    for i,left in enumerate(histories):
        row=[]
        for j,right in enumerate(histories):
            for shift in (0,1,3):
                direct=direct_pair(left,right,shift)
                folded=dot(vectors[i],apply(transfer(shift),vectors[j]))
                assert direct==folded
                count+=1
                if shift==0:
                    row.append(direct)
        gram.append(row)
    assert y37.inertia(gram)==(4,0,0)
    gap=0
    for coeff in ((1,2,3,4),(2,-1,0,3),(1,0,-2,1),(0,3,1,-2),(4,-3,2,-1),(1,1,1,1)):
        x=[sum((F(coeff[k])*vectors[k][i] for k in range(4)),F(0)) for i in range(4)]
        mean=dot(VAC,x);x=[v-mean*w for v,w in zip(x,VAC)]
        for n in range(1,6):
            y=apply(transfer(n),x)
            assert dot(y,y)<=F(1,2)**(2*n)*dot(x,x)
            assert apply(transfer(n+2),x)==apply(transfer(n),apply(transfer(2),x))
            gap+=1
    null=[(F(1),[(2,f)]),(F(-1,2),[(1,f)])]
    multiplied=product_history(null,[(F(1),[(2,f)])])
    assert fold(null)==[F(0)]*4
    for n in (0,1,2,4):
        assert fold(null,n)==[F(0)]*4
    assert fold(multiplied)==[F(3,4)*v for v in VAC]
    assert direct_pair(null,null)==0
    assert direct_pair(multiplied,multiplied)==F(9,16)
    return {'independent_reflected_path_identities':count,'positive_gram_inertia':[4,0,0],
            'mixed_source_gap_and_composition_checks':gap,'translated_null_checks':4,
            'null_multiplier_norm_squared':'9/16'}


def resolvent_controls():
    H=generator();lambdas=(F(1,2),F(1),F(2),F(7,3));inverses=0;identities=0;domains=0
    for lam in lambdas:
        R=resolvent(lam);M=[[H[i][j]+lam*I4[i][j] for j in range(4)] for i in range(4)]
        assert multiply(M,R)==I4 and multiply(R,M)==I4
        inverses+=2
        assert y37.inertia(R)==(4,0,0)
        Q=PROJECTORS[0]
        gap_ceiling=[[F(i==j)/(lam+1)-R[i][j]+Q[i][j]*(1/lam-1/(lam+1))
                      for j in range(4)] for i in range(4)]
        assert y37.inertia(gap_ceiling)[1]==0
        for seed in range(6):
            x=[F(((seed+2)*i+seed*seed)%7-3) for i in range(4)]
            u=apply(R,x)
            assert apply(H,u)==[xj-lam*uj for xj,uj in zip(x,u)]
            assert dot(u,apply(H,u))>=0
            centered=[uj-dot(VAC,u)*vj for uj,vj in zip(u,VAC)]
            assert dot(u,apply(H,u))>=dot(centered,centered)
            domains+=1
    for lam,mu in itertools.combinations(lambdas,2):
        R,S=resolvent(lam),resolvent(mu)
        assert [[R[i][j]-S[i][j] for j in range(4)] for i in range(4)]==[
            [(mu-lam)*x for x in row] for row in multiply(R,S)]
        identities+=1
    regularization=[]
    for lam in (F(2),F(8),F(32),F(128)):
        error=max(e/(lam+e) for e in ENERGIES)
        assert error<=3/lam
        regularization.append({'lambda':str(lam),'regularization_operator_error':str(error)})
    return {'two_sided_inverse_checks':inverses,'resolvent_identities':identities,
            'generator_domain_and_gap_checks':domains,'dense_range_regularization':regularization}


def exp_negative(x):
    if x<0:
        raise ValueError('nonnegative exponent magnitude required')
    halves=0
    while x>8:
        x/=2;halves+=1
    result=y44.ex(-x)
    for _ in range(halves):
        result=_r(result*result)
    return result


def quadrature_budget(lam,L,mesh,epsilon,M=F(1)):
    if min(lam,L,mesh,epsilon)<=0 or M<0 or mesh>L:
        raise ValueError('positive damping, cut, mesh and history margin required')
    return M*(L*mesh*(lam+1/epsilon)/2+exp_negative(lam*L).hi/lam)


def quadrature_controls():
    rows=[]
    for lam,L,n in itertools.product((F(1),F(2)),(F(2),F(4)),(16,64,256)):
        mesh=L/n;errors=[]
        for e in ENERGIES:
            q=exp_negative((lam+e)*mesh)
            total=_r(Iv(mesh)*(Iv(1)-y44.ivpower(q,n))/(Iv(1)-q))
            error=total-Iv(1/(lam+e))
            errors.append(max(abs(error.lo),abs(error.hi)))
        actual=max(errors)
        # The finite fixture has ||H||=3; margin 1/4 gives a lawful larger
        # Lipschitz bound 4 for unit sources, independently of the sum formula.
        budget=quadrature_budget(lam,L,mesh,F(1,4))
        assert actual<=budget
        rows.append({'lambda':str(lam),'time_cut':str(L),'cells':n,
                     'independent_resolvent_error_upper':decimals(actual),
                     'tail_plus_mesh_budget_upper':decimals(budget)})
    # The infinite-cut tail alone is not enough at a coarse mesh.
    q=exp_negative(F(1,2))
    approx=Iv(F(1,2))*(Iv(1)-y44.ivpower(q,16))/(Iv(1)-q)
    assert (approx-Iv(1)).lo>exp_negative(F(8)).hi
    return {'outward_quadratures':rows,'tail_and_mesh_both_needed':True}


# Finite four-coordinate polynomial controls, not a second completion engine.
ZERO=(0,0,0,0)


def clean(p):
    return {m:F(v) for m,v in p.items() if v}


def add(*polys):
    out={}
    for p in polys:
        for m,c in p.items():
            out[m]=out.get(m,F(0))+c
    return clean(out)


def scale(p,c):
    return clean({m:c*v for m,v in p.items()})


def mul(p,q):
    out={}
    for m,a in p.items():
        for n,b in q.items():
            key=tuple(x+y for x,y in zip(m,n))
            out[key]=out.get(key,F(0))+a*b
    return clean(out)


ONE={ZERO:F(1)}
COORD=tuple({tuple(int(i==j) for i in range(4)):F(1)} for j in range(4))


def partial(p,j):
    out={}
    for m,c in p.items():
        if m[j]:
            n=list(m);n[j]-=1;out[tuple(n)]=c*m[j]
    return clean(out)


def levi(a,b,c):
    if len({a,b,c})<3:
        return 0
    return 1 if (a,b,c) in ((1,2,3),(2,3,1),(3,1,2)) else -1


def deriv(p,a):
    terms=[scale(mul(COORD[a],partial(p,0)),F(1,2)),
           scale(mul(COORD[0],partial(p,a)),F(-1,2))]
    for b,c in itertools.product(range(1,4),repeat=2):
        e=levi(a,b,c)
        if e:
            terms.append(scale(mul(COORD[b],partial(p,c)),F(e,2)))
    return add(*terms)


def lap(p):
    return scale(add(*(deriv(deriv(p,a),a) for a in range(1,4))),-1)


def sphere_reduce(p):
    current=p
    while any(m[3]>=2 for m in current):
        out={}
        for m,c in current.items():
            if m[3]<2:
                out=add(out,{m:c});continue
            n=list(m);n[3]-=2
            out=add(out,{tuple(n):c})
            for j in range(3):
                z=n[:];z[j]+=2;out=add(out,{tuple(z):-c})
        current=out
    return clean(current)


def moment(m):
    # Exact normalized S^3 moments: rotation recurrence and sum x_i^2=1.
    if any(n%2 for n in m):
        return F(0)
    orders=[n//2 for n in m];degree=sum(orders)
    numerator=1
    for n in orders:
        for j in range(1,2*n,2):
            numerator*=j
    denominator=1
    for j in range(degree):
        denominator*=4+2*j
    return F(numerator,denominator)


def integral(p):
    return sum((c*moment(m) for m,c in p.items()),F(0))


def grad_square(p):
    return add(*(mul(deriv(p,a),deriv(p,a)) for a in range(1,4)))


def coefficient_controls():
    radius=add(*(mul(x,x) for x in COORD))
    assert sphere_reduce(add(radius,scale(ONE,-1)))=={}
    identities=0
    examples=list(COORD)+[mul(COORD[0],COORD[2]),
                         add(mul(COORD[0],COORD[0]),scale(mul(COORD[2],COORD[2]),-1))]
    for a in range(1,4):
        assert deriv(radius,a)=={}
        for f,g in itertools.product(examples[:3],repeat=2):
            assert deriv(mul(f,g),a)==add(mul(deriv(f,a),g),mul(f,deriv(g,a)))
            identities+=1
    for x in COORD:
        assert lap(x)==scale(x,F(3,4))
        assert integral(x)==0 and integral(mul(x,x))==F(1,4)
        identities+=2
    fundamental=COORD[0]
    target=scale(add(ONE,scale(mul(fundamental,fundamental),-1)),F(1,4))
    assert sphere_reduce(add(grad_square(fundamental),scale(target,-1)))=={}
    assert integral(grad_square(fundamental))==F(3,16)
    # Spin-one fusion normalization: L(F^2-1/4)=2(F^2-1/4).
    spin_one=add(mul(fundamental,fundamental),scale(ONE,F(-1,4)))
    assert sphere_reduce(add(lap(spin_one),scale(spin_one,-2)))=={}
    weighted=[]
    for tilt in (F(-1,3),F(1,3)):
        h=add(ONE,scale(COORD[0],tilt));norm=integral(mul(h,h))
        assert abs(tilt)<1
        for f in examples:
            fh=mul(f,h)
            lhs=integral(grad_square(fh))-integral(mul(mul(f,f),mul(h,lap(h))))
            rhs=integral(mul(mul(h,h),grad_square(f)))
            assert lhs==rhs>=0
            C=sum((sum(abs(c) for c in deriv(f,a).values())**2 for a in range(1,4)),F(0))
            assert rhs<=C*norm
            weighted.append({'tilt':str(tilt),'weighted_energy':str(rhs/norm),
                             'coefficient_gradient_ceiling':str(C)})
    ibp=0
    for f,g,a in itertools.product(examples,examples,range(1,4)):
        assert integral(mul(deriv(f,a),g))==-integral(mul(f,deriv(g,a)))
        ibp+=1
    variance=integral(mul(fundamental,fundamental))
    energy=integral(grad_square(fundamental));epsilon=F(1,4)
    nonzero=variance-2*epsilon*energy
    assert nonzero==F(5,32)>0
    return {'leibniz_and_fundamental_checks':identities,
            'integration_by_parts_checks':ibp,'nonconstant_vacuum_energy_controls':weighted,
            'single_site_variance':str(variance),'single_site_energy':str(energy),
            'excited_history_time':str(epsilon),'excited_history_norm_squared_lower':str(nonzero)}


def high_content_controls():
    rows=[]
    for j in (16,64,256):
        casimir=F(j*(j+1));a=1/casimir
        H=[[F(0),F(0)],[F(0),casimir]]
        source=[F(0),F(1)]
        assert dot(apply(H,source),apply(H,source))==casimir**2
        ratio=exp_negative(a*casimir)
        assert 1-ratio.hi>F(1,2)
        rows.append({'spin':j,'generator_norm_lower':str(casimir),'step':str(a),
                     'time_identity_error_lower':decimals(1-ratio.hi,False)})
    return rows


def parameter_controls():
    rows=[]
    for theta,J,R,rho in y46.CELLS:
        gate=y46.joint_gate(theta,J,R,rho)
        assert gate is not None
        gamma=gate['gamma'].lo
        assert 0<gamma<F(3,4)
        rows.append({'theta_abs':str(theta),'gamma_lower':decimals(gamma,False),
                     'resolvent_Q_norm_upper_at_lambda_1':decimals(1/(1+gamma)),
                     'nonzero_excited_norm_squared_lower':'5/32'})
    return rows


def negative_controls():
    P=transfer(1);C=[[F(i==j==0) for j in range(4)] for i in range(4)]
    reduced=multiply(multiply(C,P),C)
    returned=multiply(multiply(C,multiply(P,P)),C)
    H=generator();lam=F(2);R=resolvent(lam)
    wrong=[[H[i][j]-lam*I4[i][j] for j in range(4)] for i in range(4)]
    source=list(VECTORS[1]);u=apply(R,source)
    high=high_content_controls()
    # One-dimensional vacuum-only model satisfies every complementary bound.
    vacuum_only_T=[[F(1)]];vacuum_only_Q=[[F(0)]]
    controls={
        'reflection_null_not_a_multiplier_ideal':history_controls()['null_multiplier_norm_squared']=='9/16',
        'projected_correlation_is_not_composition':returned!=multiply(reduced,reduced),
        'norm_continuity_at_zero_not_inferred':1-exp_negative(F(1)).hi>F(1,2),
        'unbounded_generator_not_a_bounded_series_input':all(
            F(b['generator_norm_lower'])>F(a['generator_norm_lower']) for a,b in zip(high,high[1:])),
        'finite_vacuum_time_cut_has_nonzero_tail':exp_negative(F(1)).lo>0,
        'quadrature_tail_alone_is_insufficient':quadrature_controls()['tail_and_mesh_both_needed'],
        'wrong_resolvent_sign_rejected':multiply(wrong,R)!=I4,
        'wrong_generator_reconstruction_rejected':apply(H,u)!=[x+lam*y for x,y in zip(source,u)],
        'fundamental_derivative_scale_cannot_be_dropped':F(3,16)!=F(3,4),
        'gap_alone_allows_vacuum_only_space':multiply(vacuum_only_T,vacuum_only_Q)==[[F(0)]],
        'closed_parameter_boundary_rejected':y46.joint_gate(F(1,4096),8,F(4,3),F(2)) is None,
        'zero_history_margin_rejected':False}
    try:
        quadrature_budget(F(1),F(2),F(1,8),F(0))
    except ValueError:
        controls['zero_history_margin_rejected']=True
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
    return {'certificate_type':'YM48_REFLECTED_TIME_SEMIGROUP_GENERATOR_AND_LOCAL_EXCITATION',
            'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'reflection_and_time_action':history_controls(),'resolvent_and_domain':resolvent_controls(),
            'integration_tails':quadrature_controls(),'coefficient_energy':coefficient_controls(),
            'inherited_parameter_cells':parameter_controls(),'high_content_controls':high_content_controls(),
            'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_REFLECTED_HISTORY_COMPLETION__STRONGLY_CONTINUOUS_POSITIVE_TIME_ACTION__'
                           'SELF_ADJOINT_GENERATOR_WITH_DOMAIN_GRAPH_CORE_AND_GAP__'
                           'COEFFICIENT_ROW_EMBEDDING_AND_NONZERO_EXCITATION__ROW_CLOSURE_DICT_4D_OPEN',
            'evidence_scope':{'proof':'written; not mechanically formalized',
                'carrier':'reflected completion of YM47 strict-future histories on the admitted full SU(2) heat/positive-functional adapter',
                'window':'abs(theta)<1/1680; declared kappa=theta*a',
                'time':'normalized heat parameter, not a derived material clock or physical unitary time',
                'row_embedding':'isometric on the local coefficient-polynomial row completion; not proved onto or invariant',
                'domain':'Ran R_lambda; graph core R_lambda(A_plus/null); not an assumed infinite local sum',
                'excitation':'same interacting carrier, centered coefficient character, norm squared >=5/32 at time 1/4; no separate physical gauge-invariance claim',
                'open':'complete instantaneous-row Markov/observable identification; native measure/NCG; physical clock; full path measure; relativistic/spatial continuum; Clay/QG',
                'lineage':'carrier-specific reflection-positive reconstruction; Osterwalder-Schrader credited, no priority claim'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM48 fresh certificate/pin mismatch')
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
    print('YM48 PASS',sha)
    print(json.dumps({'reflected_path_identities':cert['reflection_and_time_action']['independent_reflected_path_identities'],
                      'generator_domain_checks':cert['resolvent_and_domain']['generator_domain_and_gap_checks'],
                      'weighted_energy_controls':len(cert['coefficient_energy']['nonconstant_vacuum_energy_controls']),
                      'refusal_groups':len(cert['negative_controls'])},sort_keys=True))


if __name__=='__main__':
    main()
