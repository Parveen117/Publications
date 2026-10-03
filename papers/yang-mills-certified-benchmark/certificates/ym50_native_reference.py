"""YM-50: native counted-reference and symmetric-turn heat bridge.

Python 3.12 only. Exact finite controls support the separately written proof.
Default/--check is read-only; --write regenerates this chapter only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
from ym1_certified_gap import Iv, canonical_sha
from ym6_seam_integer_dock import _r
import ym37_space_transfer as y37
import ym44_time_refinement as y44
import ym48_reflected_time as y48
import ym49_time_zero as y49
import ymf1_chain_fabric as fabric

RESULT=HERE/'YM50_RESULT.json'
PIN=HERE/'EXPECTED_YM50.sha256'
SOURCES=HERE/'YM50_SOURCE_PINS.json'
ORIGIN=HERE/'YM50_ORIGIN_LEDGER.json'
ONE=(F(1),F(0),F(0),F(0))
I4=y49.I4
mulmat=y49.multiply
apply=y49.apply
dot=y49.dot
RADIUS=y48.add(*(y48.mul(x,x) for x in y48.COORD))


def qtuple(q):
    return (q.a,q.b,q.c,q.d)


def qmul(a,b):
    return qtuple(fabric.Quat(*a)*fabric.Quat(*b))


def qdagger(q):
    return (q[0],-q[1],-q[2],-q[3])


def turns(axes=(1,2,3)):
    out=[]
    for a in axes:
        for sign in (1,-1):
            q=[F(3,5),F(0),F(0),F(0)];q[a]=F(sign*4,5);out.append(tuple(q))
    return tuple(out)


ALPHABET=(ONE,)*6+turns()


def qmatrix(q,right=False):
    basis=[tuple(F(i==j) for i in range(4)) for j in range(4)]
    cols=[qmul(b,q) if right else qmul(q,b) for b in basis]
    return y49.transpose(cols)


@lru_cache(maxsize=None)
def basis(d):
    if not isinstance(d,int) or d<0:
        raise ValueError('nonnegative integer degree required')
    return tuple((a,b,c,d-a-b-c) for a in range(d+1)
                 for b in range(d-a+1) for c in range(d-a-b+1))


def poly_power(p,n):
    out=y48.ONE
    for _ in range(n):
        out=y48.mul(out,p)
    return out


def substitute(p,A):
    linear=[{m:a for m,a in zip((next(iter(x)) for x in y48.COORD),row) if a}
            for row in A]
    powers=[[poly_power(linear[j],k) for k in range(max(
        (m[j] for m in p),default=0)+1)] for j in range(4)]
    result={}
    for m,c in p.items():
        term={y48.ZERO:c}
        for j,n in enumerate(m):
            term=y48.mul(term,powers[j][n])
        result=y48.add(result,term)
    return result


def evaluate(p,x):
    return sum((c*prod(xj**n for xj,n in zip(x,m)) for m,c in p.items()),F(0))


def prod(values):
    out=F(1)
    for v in values:
        out*=v
    return out


def vector(p,d):
    if any(sum(m)!=d for m in p):
        raise ValueError('homogeneous degree mismatch')
    return [p.get(m,F(0)) for m in basis(d)]


def polynomial(v,d):
    return {m:c for m,c in zip(basis(d),v) if c}


def action_matrix(d,A):
    columns=[vector(substitute({m:F(1)},A),d) for m in basis(d)]
    return y49.transpose(columns)


@lru_cache(maxsize=None)
def degree_data(d):
    mons=basis(d);size=len(mons)
    weights=tuple(F(prod(factorial(n) for n in m),factorial(d)) for m in mons)
    G=y49.diag(weights);I=y49.identity(size)
    actions=[action_matrix(d,qmatrix(q)) for q in turns()]
    P=y49.scale(I,F(1,2))
    for A in actions:
        P=y49.add(P,y49.scale(A,F(1,12)))
    radial=vector(poly_power(RADIUS,d//2),d) if d%2==0 else [F(0)]*size
    norm=dot(radial,apply(G,radial))
    pi=[[radial[i]*weights[j]*radial[j]/norm if norm else F(0)
         for j in range(size)] for i in range(size)]
    S=y49.add(I,P,-1)
    inverse=y49.inverse(y49.add(S,pi))
    return {'basis':mons,'G':G,'P':P,'Pi':pi,'S':S,'inverse':inverse,
            'actions':actions,'radial':radial}


def native_moment(m):
    """Formula derived in T4; matrix projections independently check it."""
    if any(n%2 for n in m):
        return F(0)
    numerator=1
    for n in m:
        for j in range(1,n,2):
            numerator*=j
    degree=sum(m)//2
    denominator=prod(range(4,2*degree+3,2))
    return F(numerator,denominator)


def phi(p):
    return sum((c*native_moment(m) for m,c in p.items()),F(0))


def split_degrees(p):
    return {d:{m:c for m,c in p.items() if sum(m)==d}
            for d in sorted({sum(m) for m in p})}


def poisson(p):
    answer={};bound=F(0)
    for d,part in split_degrees(p).items():
        data=degree_data(d);v=vector(part,d)
        centered=[a-b for a,b in zip(v,apply(data['Pi'],v))]
        u=apply(data['inverse'],centered)
        assert apply(data['S'],u)==centered
        assert apply(data['Pi'],u)==[F(0)]*len(u)
        answer=y48.add(answer,polynomial(u,d));bound+=sum(abs(x) for x in u)
    return answer,bound


def cesaro(p,N):
    if not isinstance(N,int) or N<1:
        raise ValueError('positive integer average depth required')
    out={}
    for d,part in split_degrees(p).items():
        data=degree_data(d);v=vector(part,d);total=[F(0)]*len(v)
        for _ in range(N):
            total=[a+b for a,b in zip(total,v)]
            v=apply(data['P'],v)
        out=y48.add(out,polynomial([v/N for v in total],d))
    return out


@lru_cache(maxsize=None)
def endpoints(k):
    if not isinstance(k,int) or k<0:
        raise ValueError('nonnegative integer word depth required')
    if k==0:
        return (ONE,)
    return tuple(qmul(g,x) for x in endpoints(k-1) for g in ALPHABET)


def direct_count(p,k):
    return sum((evaluate(p,x) for x in endpoints(k)),F(0))/12**k


def native_generator(p,a):
    q=[F(0)]*4;q[a]=F(-1,2);A=qmatrix(tuple(q))
    out={}
    for j in range(4):
        field={next(iter(y48.COORD[k])):A[j][k] for k in range(4) if A[j][k]}
        out=y48.add(out,y48.mul(field,y48.partial(p,j)))
    return out


def lap(p):
    return y48.scale(y48.add(*(native_generator(native_generator(p,a),a)
                              for a in range(1,4))),-1)


def quaternion_controls():
    points=[ONE,*turns(),qtuple(fabric.rational_unit(F(1,3),F(-1,4),F(2,5)))]
    checks=0
    for q in points:
        A=qmatrix(q);B=qmatrix(q,True)
        assert mulmat(y49.transpose(A),A)==I4
        assert mulmat(y49.transpose(B),B)==I4
        assert qmul(q,qdagger(q))==ONE
        for r in points:
            matrix_product=fabric.mmul(fabric.to_matrix(fabric.Quat(*q)),
                                      fabric.to_matrix(fabric.Quat(*r)))
            assert fabric.meq(matrix_product,fabric.to_matrix(fabric.Quat(*qmul(q,r))))
            checks+=1
    # Infinite-order argument's exact monic-recurrence congruence controls.
    prev,cur=F(2),F(6,5);congruences=0
    for n in range(1,17):
        value=cur*5**n
        assert value.denominator==1 and value.numerator%5==1
        prev,cur=cur,F(6,5)*cur-prev;congruences+=1
    return {'native_matrix_product_checks':checks,'infinite_order_congruence_controls':congruences}


def finite_energy_controls():
    rows=[]
    for d in range(5):
        data=degree_data(d);G,P,Pi,S=(data[k] for k in ('G','P','Pi','S'))
        size=len(P);I=y49.identity(size);Z=y49.scale(I,F(0))
        assert mulmat(G,P)==mulmat(y49.transpose(P),G)
        assert mulmat(Pi,Pi)==Pi and mulmat(P,Pi)==Pi
        defect=mulmat(G,S);square_sum=Z
        for index in range(0,6,2):
            A=data['actions'][index];Ad=data['actions'][index+1]
            assert mulmat(A,Ad)==I and mulmat(mulmat(y49.transpose(A),G),A)==G
            diff=y49.add(I,A,-1)
            square_sum=y49.add(square_sum,y49.scale(
                mulmat(mulmat(y49.transpose(diff),G),diff),F(1,12)))
        assert defect==square_sum
        inertia=y37.inertia(defect)
        expected_null=int(d%2==0)
        assert inertia==(size-expected_null,0,expected_null)
        assert mulmat(y49.add(S,Pi),data['inverse'])==I
        rows.append({'degree':d,'dimension':size,'stationary_dimension':expected_null,
                     'defect_inertia':list(inertia)})
    return {'finite_degree_spaces':rows,'all_checked_fixed_spaces_are_radial':True}


def counted_readout_controls():
    x=y48.COORD
    samples=[y48.ONE,x[0],y48.mul(x[0],x[0]),y48.mul(x[1],x[2]),
             poly_power(x[0],4),y48.add(x[0],y48.scale(y48.mul(x[1],x[1]),F(2,3)))]
    checks=padding=0
    for p in samples:
        for k in range(4):
            transferred={}
            for d,part in split_degrees(p).items():
                v=apply(y49.power(degree_data(d)['P'],k),vector(part,d))
                transferred=y48.add(transferred,polynomial(v,d))
            assert direct_count(p,k)==evaluate(transferred,ONE);checks+=1
        for N in (1,2,4):
            independent=sum((direct_count(p,k) for k in range(N)),F(0))/N
            assert independent==evaluate(cesaro(p,N),ONE)
            # Count the padded atlas, retaining all multiplicities.
            padded=sum((12**(N-1-k)*sum((evaluate(p,x) for x in endpoints(k)),F(0))
                        for k in range(N)),F(0))
            assert padded/(N*12**(N-1))==independent;padding+=1
    p,q,r,s=x[0],x[1],x[2],y48.add(y48.ONE,x[3])
    paired_real=y48.add(y48.mul(p,r),y48.mul(q,s))
    paired_imag=y48.add(y48.mul(p,s),y48.scale(y48.mul(q,r),-1))
    real=imag=F(0);N=4
    for k in range(N):
        for point in endpoints(k):
            a=fabric.GQ(evaluate(p,point),evaluate(q,point))
            b=fabric.GQ(evaluate(r,point),evaluate(s,point))
            product=a.conj()*b
            real+=product.x/(N*12**k);imag+=product.y/(N*12**k)
    assert real==evaluate(cesaro(paired_real,N),ONE)
    assert imag==evaluate(cesaro(paired_imag,N),ONE)
    return {'independent_word_transfer_checks':checks,'padded_native_trace_checks':padding,
            'complex_native_pairing_checks':2,'largest_raw_word_count':len(endpoints(3))}


def tail_controls():
    x=y48.COORD
    samples=[x[0],y48.mul(x[0],x[0]),poly_power(x[0],3),poly_power(x[0],4),
             y48.add(x[0],y48.scale(y48.mul(x[1],x[2]),F(2,3)))]
    points=[ONE,*turns(),qtuple(fabric.rational_unit(F(2,3),F(1,5),F(-1,4)))]
    rows=[];checks=0
    for index,p in enumerate(samples):
        u,B=poisson(p);target=phi(p)
        for N in (1,3,9,27):
            averaged=cesaro(p,N)
            remainder=y48.add(averaged,y48.scale(y48.ONE,-target))
            pu={}
            for d,part in split_degrees(u).items():
                pu=y48.add(pu,polynomial(apply(y49.power(degree_data(d)['P'],N),
                                                      vector(part,d)),d))
            telescoped=y48.scale(y48.add(u,y48.scale(pu,-1)),F(1,N))
            assert y48.sphere_reduce(y48.add(remainder,y48.scale(telescoped,-1)))=={}
            for point in points:
                assert abs(evaluate(remainder,point))<=2*B/N;checks+=1
            rows.append({'source':index,'depth':N,'coefficient_tail_constant':str(B),
                         'uniform_error_upper':str(2*B/N),
                         'identity_readout_error':str(abs(evaluate(remainder,ONE)))})
    return {'pointwise_tail_checks':checks,'exact_telescoping_cases':len(rows),'cases':rows}


def moment_and_bridge_controls():
    projection_checks=legacy=rotation=multiplier=0
    for d in range(5):
        data=degree_data(d)
        for j,m in enumerate(basis(d)):
            column=[row[j] for row in data['Pi']]
            assert evaluate(polynomial(column,d),ONE)==native_moment(m)
            projection_checks+=1
    for d in range(9):
        for m in basis(d):
            assert native_moment(m)==y48.moment(m);legacy+=1
    unit=qtuple(fabric.rational_unit(F(1,3),F(-1,4),F(2,5)))
    for d in range(5):
        for m in basis(d):
            p={m:F(1)}
            for right in (False,True):
                assert phi(substitute(p,qmatrix(unit,right)))==phi(p);rotation+=1
    mons=[m for d in range(4) for m in basis(d)]
    gram=[[native_moment(tuple(a+b for a,b in zip(m,n))) for n in mons] for m in mons]
    assert y37.inertia(gram)==(30,0,5)
    zero=y48.add(RADIUS,y48.scale(y48.ONE,-1))
    assert phi(y48.mul(zero,zero))==0
    samples=[y48.ONE,*y48.COORD,
             y48.add(y48.COORD[0],y48.scale(y48.mul(y48.COORD[1],y48.COORD[2]),F(2,3)))]
    for p,f in itertools.product(samples,repeat=2):
        product=y48.mul(p,f);ceiling=sum(abs(c) for c in p.values())
        assert phi(y48.mul(product,product))<=ceiling**2*phi(y48.mul(f,f))
        multiplier+=1
    return {'independent_projected_moments':projection_checks,'legacy_moment_matches':legacy,
            'left_right_invariance_checks':rotation,'recognition_gram_inertia':[30,0,5],
            'bounded_multiplier_checks':multiplier,'unit_relation_null_energy':True}


def cosine_square(squared,terms=24):
    if squared<0 or squared>F(3,2):
        raise ValueError('certified cosine-square range is [0,3/2]')
    total=F(0)
    for k in range(terms+1):
        total+=(-1)**k*squared**k/factorial(2*k)
    nxt=(-1)**(terms+1)*squared**(terms+1)/factorial(2*terms+2)
    return _r(Iv(min(total,total+nxt),max(total,total+nxt)))


def heat_budget(d,t,n):
    if not isinstance(d,int) or d<0 or t<0 or not isinstance(n,int) or n<1:
        raise ValueError('nonnegative degree/time and positive integer refinement required')
    return F(3,8)*d**4*t*t/n


def heat_controls():
    identities=skew=0
    for d in range(5):
        data=degree_data(d);G=data['G']
        for a in range(1,4):
            D=y49.transpose([vector(native_generator({m:F(1)},a),d) for m in basis(d)])
            assert y49.add(mulmat(y49.transpose(D),G),mulmat(G,D))==y49.scale(G,0)
            skew+=1
        for m in basis(d):
            p={m:F(1)}
            euclidean=y48.add(*(y48.partial(y48.partial(p,j),j) for j in range(4)))
            expected=y48.scale(y48.add(y48.scale(p,d*(d+2)),
                                       y48.scale(y48.mul(RADIUS,euclidean),-1)),F(1,4))
            assert lap(p)==expected==y48.lap(p);identities+=1
    x0=y48.COORD[0]
    harmonics=[x0,
        y48.add(poly_power(x0,2),y48.scale(RADIUS,F(-1,4))),
        y48.add(poly_power(x0,3),y48.scale(y48.mul(x0,RADIUS),F(-1,2))),
        y48.add(poly_power(x0,4),y48.scale(y48.mul(poly_power(x0,2),RADIUS),F(-3,4)),
                 y48.scale(poly_power(RADIUS,2),F(1,16)))]
    for ell,p in enumerate(harmonics,1):
        assert lap(p)==y48.scale(p,F(ell*(ell+2),4))
    rows=[]
    for d,t,n in itertools.product((1,2),(F(1,4),F(1),F(2)),(8,32,128)):
        squared=6*t/n
        if d==1:
            sample=cosine_square(squared/4);energy=F(3,4)
        else:
            sample=_r((Iv(1)+Iv(2)*cosine_square(squared))/Iv(3));energy=F(2)
        approximate=y44.ivpower(sample,n);target=y48.exp_negative(energy*t)
        error=approximate-target;error_upper=max(abs(error.lo),abs(error.hi))
        bound=heat_budget(d,t,n)
        assert error_upper<=bound
        rows.append({'degree':d,'time':str(t),'refinement':n,
                     'independent_scalar_error_upper':str(error_upper),
                     'general_coefficient_bound':str(bound)})
    return {'skew_generator_checks':skew,'native_casimir_polynomial_checks':identities,
            'harmonic_eigenvalue_checks':4,'outward_heat_refinements':rows}


def refusal_controls():
    x0=y48.COORD[0];square=y48.mul(x0,x0)
    P=degree_data(2)['P']
    expected=y48.add(y48.scale(square,F(43,75)),y48.scale(RADIUS,F(8,75)))
    assert polynomial(apply(P,vector(square,2)),2)==expected
    # Nested raw-record energies, calculated from the exact conditional recursion.
    raw=[]
    for n in (1,4,16,64):
        variance=F(1,4)+F(3,4)*F(43,75)**n
        next_variance=F(1,4)+F(3,4)*F(43,75)**(n+1)
        difference=variance+next_variance-2*F(4,5)*variance
        assert difference==F(1,10)-F(1,50)*F(43,75)**n
        raw.append({'depth':n,'raw_record_step_energy':str(difference)})
    lost=y48.add(y48.mul(y48.COORD[2],y48.COORD[2]),
                  y48.mul(y48.COORD[3],y48.COORD[3]))
    one_axis=qmatrix(turns((1,))[0])
    assert substitute(lost,one_axis)==lost and evaluate(lost,ONE)==0 and phi(lost)==F(1,2)
    wrong=y48.exp_negative(F(1,8));right=y48.exp_negative(F(3,4))
    controls={
        'single_axis_does_not_select_full_reference':phi(lost)==F(1,2),
        'finite_depth_is_not_exact_invariance':evaluate(cesaro(x0,4),ONE)>0==phi(x0),
        'raw_record_sequence_is_not_recognition_cauchy':all(
            F(2,25)<=F(r['raw_record_step_energy'])<F(1,10) for r in raw),
        'endpoint_equality_does_not_erase_native_path_memory':
            qmul(turns()[0],turns()[1])==ONE and F(1)**2+F(-1)**2==2,
        'unnormalized_constant_readout_is_wrong':len(endpoints(3))!=1,
        'unweighted_depth_union_is_not_cesaro':
            sum((sum((evaluate(x0,p) for p in endpoints(k)),F(0)) for k in range(4)),F(0))
            /sum(12**k for k in range(4))!=evaluate(cesaro(x0,4),ONE),
        'wrong_heat_clock_factor_is_detected':wrong.lo>right.hi,
        'single_axis_heat_generator_is_different':
            y48.scale(native_generator(native_generator(x0,1),1),-1)!=lap(x0),
        'heat_tail_is_not_degree_uniform':heat_budget(4,F(1),16)>heat_budget(1,F(1),16),
        'zero_count_depth_rejected':False,
        'invalid_heat_refinement_rejected':False,
        'out_of_range_cosine_enclosure_rejected':False}
    for key,call in (
        ('zero_count_depth_rejected',lambda:cesaro(x0,0)),
        ('invalid_heat_refinement_rejected',lambda:heat_budget(1,F(1),0)),
        ('out_of_range_cosine_enclosure_rejected',lambda:cosine_square(F(2)))):
        try:
            call()
        except ValueError:
            controls[key]=True
    assert all(controls.values())
    return {'controls':controls,'raw_record_non_cauchy_samples':raw,'raw_record_energy_limit':'1/10'}


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
    return {'certificate_type':'YM50_NATIVE_COUNTED_REFERENCE_AND_HEAT_INTERTWINER',
        'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
        'origin_ledger':json.loads(ORIGIN.read_text()),
        'native_quaternion_controls':quaternion_controls(),
        'finite_coefficient_energy':finite_energy_controls(),
        'native_counted_readouts':counted_readout_controls(),
        'constructive_average_tails':tail_controls(),
        'reference_recognition_bridge':moment_and_bridge_controls(),
        'native_heat_bridge':heat_controls(),'refusals':refusal_controls(),
        'claim_status':'WRITTEN_NATIVE_COUNTED_RECORD_REFERENCE_CONSTRUCTION__'
            'COEFFICIENT_AND_UNIFORM_REFERENCE_IDENTITY__SYMMETRIC_TURN_HEAT_INTERTWINER__'
            'DECLARED_CHAIN_TRANSFERS_PRESERVED__GENERAL_UGD_PHYSICAL_SELECTION_4D_OPEN',
        'evidence_scope':{'proof':'written; not mechanically formalized',
            'native_source':'normalized Phi_Sigma on specified finite diagonal record kernels; rational counts and derived quaternion coefficients',
            'functional':'polynomial readouts and their uniform completion, with 2*B_p/N tail',
            'state_module':'new positive-form completion of readout limits; not convergence of raw endpoint records in the older refinement module',
            'heat':'specified symmetric native-turn protocol; finite-content error 3*d^4*t^2/(8*n)',
            'intertwining':'coefficient recognition completion, bounded multipliers, declared finite-chain heat/bridge transfers and A_unif-generated history sector',
            'open':'general RH E5C/E6; physical state/clock/trajectory selection; NCG quantum measure; actual row defect; full bounded-history density; 4D/AF/Clay/QG'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM50 fresh certificate/pin mismatch')
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
    print('YM50 PASS',sha)
    print(json.dumps({'word_checks':cert['native_counted_readouts']['independent_word_transfer_checks'],
        'projected_moments':cert['reference_recognition_bridge']['independent_projected_moments'],
        'tail_cases':cert['constructive_average_tails']['exact_telescoping_cases'],
        'heat_refinements':len(cert['native_heat_bridge']['outward_heat_refinements']),
        'refusal_groups':len(cert['refusals']['controls']),'general_UGD_and_4D':'OPEN'},sort_keys=True))


if __name__=='__main__':
    main()
