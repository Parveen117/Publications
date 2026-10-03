"""YM54 exact response/curvature/record controls. Python 3.12 only.

General results have written proofs. These finite controls reuse the existing
matrix, quaternion, polynomial and interval engines; no new operator engine.
Default/--check is read-only; --write touches only YM54 result and pin.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from math import factorial
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
import ym53_anisotropic_interaction as chain

z=chain.z
q=z.q
y=z.y
a=z.a
p=z.p
Iv=chain.Iv
exact=q.exact
mm=a.multiply
I=a.identity(2)
ZERO=a.scale(I,0)
K=a.diag((1,-1))
R=[[F(0),F(-1)],[F(1),F(0)]]
L=mm(R,K)
RESULT=HERE/'YM54_RESULT.json'
PIN=HERE/'EXPECTED_YM54.sha256'
SOURCES=HERE/'YM54_SOURCE_PINS.json'


def matrix2(M):
    if len(M)!=2 or any(len(row)!=2 for row in M):
        raise ValueError('two-mode exact matrix required')
    return [[exact(x) for x in row] for row in M]


def trace(M):return sum((M[i][i] for i in range(len(M))),F(0))


def det2(M):return M[0][0]*M[1][1]-M[0][1]*M[1][0]


def comm(A,B):return a.add(mm(A,B),mm(B,A),-1)


def conjugate(B,X):return mm(mm(B,X),a.inverse(B))


def realification(real,imag):
    # A faithful finite control of the already admitted central iota extension.
    return [real[i]+[-x for x in imag[i]] for i in range(2)]+[
        imag[i]+real[i] for i in range(2)]


def embed_quaternion(v):
    r=a.add(a.scale(I,v[0]),a.scale(R,v[3]))
    s=a.add(a.scale(K,v[1]),a.scale(L,v[2]))
    return realification(r,s)


def response_data(H,derivatives,B,factor_scale=F(1)):
    H=matrix2(H);B=matrix2(B);factor_scale=exact(factor_scale)
    if H!=a.transpose(H) or H[0][0]<=0 or det2(H)<=0:
        raise ValueError('positive invertible native response required')
    if len(derivatives)!=2:raise ValueError('two declared response directions required')
    ds=[matrix2(D) for D in derivatives]
    if any(D!=a.transpose(D) for D in ds):
        raise ValueError('response derivatives must preserve dagger')
    if factor_scale<=0 or mm(a.transpose(B),B)!=a.scale(H,factor_scale):
        raise ValueError('factor must satisfy B^dagger B = positive_scale H')
    invH=a.inverse(H)
    X=[mm(invH,D) for D in ds]
    sigma=[trace(x)/2 for x in X]
    shape=[a.add(x,a.scale(I,-s)) for x,s in zip(X,sigma)]
    Y=[conjugate(B,x) for x in shape]
    assert all(M==a.transpose(M) and trace(M)==0 for M in Y)
    vectors=[(M[0][0],M[0][1],F(0)) for M in Y]
    G=[[trace(mm(x,v)) for v in shape] for x in shape]
    curvature=a.scale(comm(*X),F(-1,4))
    Z=conjugate(B,curvature)
    orientation=1 if det2(B)>0 else -1
    marker=orientation*Z[1][0]
    E=trace(mm(a.transpose(Z),Z))
    protocol=q.protocol([(F(1,2),v) for v in vectors])
    C,m2,m4=q.moments(protocol)
    tau=trace(G)
    beta=4*marker*marker/tau if tau else F(0)
    return {'H':H,'derivatives':ds,'factor':B,'factor_scale':factor_scale,
            'X':X,'sigma':sigma,'shape':shape,'Y':Y,'vectors':vectors,
            'G':G,'tau':tau,'F':curvature,'Z':Z,'marker':marker,'E':E,
            'protocol':protocol,'C':C,'m2':m2,'m4':m4,'beta_floor':beta}


def fixture():
    return response_data(((2,1),(1,5)),(((0,1),(1,0)),((1,0),(0,0))),
                         ((2,1),(0,3)),F(2))


def from_shapes(vectors,B=I,slopes=(F(0),F(0))):
    B=matrix2(B);H=mm(a.transpose(B),B)
    Y=[a.add(a.scale(K,exact(v[0])),a.scale(L,exact(v[1]))) for v in vectors]
    ds=[a.add(mm(mm(a.transpose(B),M),B),a.scale(H,exact(s)))
        for M,s in zip(Y,slopes)]
    return response_data(H,ds,B)


def native_embedding_controls():
    assert mm(K,K)==I and mm(L,L)==I and mm(R,R)==a.scale(I,-1)
    assert comm(K,L)==a.scale(R,-2)
    basis=[tuple(F(i==j) for i in range(4)) for j in range(4)]
    products=0
    for x,w in itertools.product(basis,repeat=2):
        assert mm(embed_quaternion(x),embed_quaternion(w))==embed_quaternion(y.qmul(x,w))
        products+=1
    for x in basis:
        assert a.transpose(embed_quaternion(x))==embed_quaternion(y.qdagger(x))
    for v in ((F(1),F(0)),(F(2,3),F(-4,5)),(F(1,9),F(1,3))):
        Q=embed_quaternion((F(0),v[0],v[1],F(0)))
        assert a.transpose(Q)==a.scale(Q,-1)
        assert mm(Q,Q)==a.scale(a.identity(4),-sum(t*t for t in v))
    # A real shape word is not itself an imaginary compact-turn generator.
    assert mm(realification(K,ZERO),realification(K,ZERO))==a.identity(4)
    return {'quaternion_product_checks':products,'dagger_checks':4,
            'shape_lift_checks':3,'unlifted_shape_refused_as_skew_turn':True}


def response_controls():
    factors=(I,a.diag((2,3)),[[F(2),F(1)],[F(0),F(3)]],
             [[F(1),F(2)],[F(-1),F(1)]])
    pairs=(((1,0),(0,1)),((1,2),(3,-1)),((1,0),(0,F(1,10))),
           ((1,2),(2,4)),((0,0),(1,1)),((0,0),(0,0)))
    checks=positive=degenerate=0
    for B,pair,slopes in itertools.product(factors,pairs,((0,0),(1,-2))):
        d=from_shapes(pair,B,slopes);H=d['H'];G=d['G']
        for X,D in zip(d['X'],d['derivatives']):
            assert mm(a.transpose(X),H)==mm(H,X)==D
        full=[[trace(mm(x,w)) for w in d['X']] for x in d['X']]
        assert full==a.add(G,[[2*u*v for v in d['sigma']] for u in d['sigma']])
        V=a.transpose([list(v) for v in d['vectors']])
        assert G==a.scale(mm(a.transpose(V),V),2)
        assert d['C']==a.scale(mm(V,a.transpose(V)),F(1,2))
        assert det2(G)==8*d['E']==16*d['marker']**2
        assert a.psd(G) and a.psd(d['C'])
        assert d['m2']==d['tau']/4
        assert d['m4']==(G[0][0]**2+G[1][1]**2)/8
        if d['marker']:
            assert a.psd(a.add(a.scale(G,F(1,4)),a.scale(I,-d['beta_floor'])))
            assert a.y37.inertia(d['C'])==(2,0,1)
            assert not a.psd(a.add(d['C'],a.scale(a.identity(3),-d['beta_floor'])))
            positive+=1
        else:
            assert not det2(G) and d['beta_floor']==0
            degenerate+=1
        checks+=1
    H=a.diag((2,1));X=mm(a.inverse(H),L)
    assert trace(mm(X,X))==1 and trace(mm(a.transpose(X),X))==F(5,4)
    return {'response_jets':checks,'rank_two_floor_checks':positive,
            'degenerate_shape_cases':degenerate,
            'wrong_unweighted_adjoint_equality_rejected':True}


def counted_protocol_controls():
    examples=[fixture(),from_shapes(((1,0),(0,2))),
              from_shapes(((1,1),(2,-1))),from_shapes(((1,2),(2,4)))]
    literal=actions=positive=0
    for d in examples:
        words=[tuple(sign*t for t in v) for v in d['vectors'] for sign in (1,-1)]
        assert len(words)==4
        assert [sum(w[j] for w in words) for j in range(3)]==[0]*3
        actual=[[sum(w[i]*w[j] for w in words)/4 for j in range(3)] for i in range(3)]
        assert actual==d['C'];literal+=1
        for degree in range(4):
            D=q.generators(degree);size=len(D[0]);direct=a.scale(a.identity(size),0)
            for v in words:
                direction=a.scale(a.identity(size),0)
                for coefficient,base in zip(v,D):direction=a.add(direction,a.scale(base,coefficient))
                direct=a.add(direct,a.scale(mm(direction,direction),F(-1,4)))
            operator=q.operator(degree,d['C'])
            assert direct==operator
            gram=y.degree_data(degree)['G']
            assert a.psd(mm(gram,operator))
            actions+=1;positive+=1
    return {'four_label_moment_checks':literal,'independent_generator_checks':actions,
            'finite_coefficient_positive_energy_checks':positive,'degree_range':[0,3]}


def covariance_controls():
    transformations=[([[F(3,5),F(-4,5)],[F(4,5),F(3,5)]],(F(3,5),0,0,F(4,5))),
                     (K,(F(0),F(1),F(0),F(0)))]
    cases=0
    for d in (fixture(),from_shapes(((1,2),(3,-1)),[[2,1],[0,3]],(1,-2))):
        for U,g in transformations:
            g=tuple(map(F,g));rot=q.rotation(g)
            new=response_data(d['H'],d['derivatives'],mm(U,d['factor']),d['factor_scale'])
            assert new['G']==d['G'] and new['marker']==d['marker']
            assert new['C']==q.rotate_tensor(d['C'],rot)
            cases+=1
        O=transformations[0][0]
        ds=[a.add(a.scale(d['derivatives'][0],row[0]),a.scale(d['derivatives'][1],row[1])) for row in O]
        new=response_data(d['H'],ds,d['factor'],d['factor_scale'])
        assert new['C']==d['C']
        assert new['G']==mm(mm(O,d['G']),a.transpose(O));cases+=1
    d1=from_shapes(((1,0),(0,2)))
    d2=from_shapes(((F(3,5),F(8,5)),(F(-4,5),F(6,5))))
    assert d1['C']==d2['C'] and d1['m4']!=d2['m4']
    assert d1['m4']/384!=d2['m4']/384
    return {'factor_and_direction_covariances':cases,
            'same_heat_different_finite_records':{'first_m4':str(d1['m4']),
              'second_m4':str(d2['m4']),'scalar_fourth_coefficient_difference':str((d1['m4']-d2['m4'])/384)}}


def cos_sqrt_interval(x,N=20):
    x=exact(x);chain.integer(N,1)
    if not 0<=x<=1:raise ValueError('control requires 0<=cosine squared argument<=1')
    partial=sum(((-x)**k/factorial(2*k) for k in range(N+1)),F(0))
    tail=x**(N+1)/factorial(2*N+2)/(1-x/F((2*N+4)*(2*N+3)))
    return Iv(partial-tail,partial+tail)


def refinement_controls():
    rows=[]
    examples=[fixture(),from_shapes(((1,0),(0,2))),from_shapes(((1,1),(2,-1)))]
    for index,d in enumerate(examples):
        for t in (F(1,8),F(1,2),F(1)):
            for n in (4,16,64,256):
                # Independent actual four-label readout on every linear coefficient.
                one=Iv(0)
                for v in d['vectors']:
                    norm=sum(x*x for x in v)
                    one=one+cos_sqrt_interval(norm*t/(2*n))/Iv(2)
                finite=chain.refine.ivpower(one,n)
                target=chain.refine.ex(-d['m2']*t/4)
                error=finite-target
                actual=max(abs(error.lo),abs(error.hi))
                bound=q.budget(1,t,n,d['m2'],d['m4'])
                coarse=F(5,1536)*t*t*d['tau']**2/n
                assert actual<=bound<=coarse
                rows.append({'protocol':index,'time':str(t),'records_steps':n,
                    'linear_error_upper':chain.old.decimals(actual),
                    'written_degree_one_budget':chain.old.decimals(bound)})
    return {'independent_finite_turn_enclosures':len(rows),'cells':rows,
            'finite_degree_bound':'5*d^4*t^2*T^2/(1536*n)'}


def curvature_budget(floor_abs,T):
    floor_abs=exact(floor_abs);T=exact(T)
    if floor_abs<=0 or T<8*floor_abs:
        raise ValueError('require a positive curvature floor and nonempty trace budget T>=8*f0')
    return 4*floor_abs*floor_abs/T


def source_gate(data,floor_abs,T,theta,J,R):
    beta=curvature_budget(floor_abs,T)
    if abs(data['marker'])<floor_abs or data['tau']>T:
        raise ValueError('actual response fails the stated uniform source bounds')
    gate=chain.block_gate(beta,theta,J,R)
    return None if gate is None else dict(gate,beta_floor=beta,trace_ceiling=exact(T)/4)


def fixture_controls():
    d=fixture()
    assert d['vectors']==[(F(1,9),F(1,3),F(0)),(F(2,9),F(-1,6),F(0))]
    assert d['G']==[[F(20,81),F(-5,81)],[F(-5,81),F(25,162)]]
    assert d['marker']==F(-5,108) and d['tau']==F(65,162)
    assert d['C']==a.diag((F(5,162),F(5,72),F(0)))
    assert d['beta_floor']==F(5,234)
    free=z.rates((F(0),F(5,162),F(5,72)))
    assert free[0]==F(65,2592) and free[1]==F(5,162)
    assert free[0]!=d['tau']/2
    theta=F(1,262144);gate=source_gate(d,F(5,108),F(65,162),theta,6,F(5,4))
    assert gate is not None and gate['gamma'].lo>F(79467076,10**12)
    assert gate['beta_floor']/2400==F(1,112320)
    # Independently differentiate the original four-variable polynomial embedding.
    s,v=p.COORD[0],p.COORD[1]
    potential=p.add(p.scale(s,10),p.scale(v,-10),p.mul(s,s),p.mul(s,v),
                    p.scale(p.mul(v,v),F(5,2)),p.scale(p.mul(p.mul(s,s),v),F(1,2)))
    origin=(F(0),)*4
    H=[[y.evaluate(p.partial(p.partial(potential,i),j),origin) for j in range(2)] for i in range(2)]
    Ds=[[[y.evaluate(p.partial(p.partial(p.partial(potential,i),j),k),origin)
          for j in range(2)] for i in range(2)] for k in range(2)]
    assert H==d['H'] and Ds==d['derivatives']
    flat=response_data(H,[ZERO,ZERO],d['factor'],2)
    assert flat['H']==H and flat['C']==a.scale(q.I3,0) and not flat['marker']
    return {'H':[[str(v) for v in row] for row in H],
        'vectors':[list(map(str,v)) for v in d['vectors']],
        'shape_gram':[[str(v) for v in row] for row in d['G']],
        'marker':str(d['marker']),'tau':str(d['tau']),
        'C':[[str(v) for v in row] for row in d['C']],
        'exact_second_eigenvalue':'5/162','beta_from_curvature':'5/234',
        'sufficient_window':'abs(theta) < 1/112320','theta':str(theta),'J':6,'R':'5/4',
        'weighted_margin':str(gate['margin']),
        'gamma_lower':chain.old.decimals(gate['gamma'].lo,False),
        'gamma_upper':chain.old.decimals(gate['gamma'].hi),
        'free_full_rate':str(free[0]),'free_even_rate':str(free[1]),
        'point_H_does_not_select_protocol':True,'potential_derivative_replay':True}


def scale_and_failure_controls():
    d=fixture();checks=0
    for multiplier in (F(1,10),F(2),F(7,3)):
        # Local positive scalar scaling: choose a square multiplier for exact B.
        s=multiplier**2;slopes=(F(2),F(-3))
        H=a.scale(d['H'],s)
        Ds=[a.add(a.scale(D,s),a.scale(d['H'],s*k)) for D,k in zip(d['derivatives'],slopes)]
        new=response_data(H,Ds,a.scale(d['factor'],multiplier),2)
        assert new['C']==d['C'] and new['G']==d['G'] and new['marker']==d['marker']
        changed=response_data(d['H'],[a.scale(D,multiplier) for D in d['derivatives']],d['factor'],2)
        assert changed['C']==a.scale(d['C'],multiplier**2)
        assert changed['marker']==multiplier**2*d['marker']
        assert changed['beta_floor']==multiplier**2*d['beta_floor']
        checks+=1
    slow=[];fixed_curvature=[]
    for eps in (F(1),F(1,10),F(1,100)):
        small=from_shapes(((1,0),(0,eps)))
        assert small['C']==a.diag((F(1,2),eps**2/2,0))
        assert small['tau']/2==1+eps**2
        slow.append({'epsilon':str(eps),'second_eigenvalue':str(eps**2/2),
                     'half_gram_trace':str(small['tau']/2)})
        t=1/eps;fixed=from_shapes(((t,0),(0,1/t)))
        assert fixed['marker']==F(1,2)
        fixed_curvature.append({'t':str(t),'marker':'1/2','tau':str(fixed['tau']),
                                'second_eigenvalue':str(eps**2/2)})
    central=response_data(I,[I,K],I)
    assert central['G']==[[0,0],[0,2]] and not central['marker']
    idle=[]
    for rho in (F(1,10),F(1,100),F(1,1000)):
        proto=q.protocol([(1-rho,(0,0,0))]+[(rho*w,v) for w,v in d['protocol']])
        C,m2,m4=q.moments(proto)
        assert C==a.scale(d['C'],rho) and m2==rho*d['m2']
        idle.append({'active_weight':str(rho),'trace':str(m2)})
    return {'scalar_and_clock_scale_checks':checks,'independent_but_slow':slow,
            'fixed_curvature_without_budget':fixed_curvature,'fixed_response_idle_protocols':idle,
            'central_response_does_not_restore_shape_rank':True,
            'physical_protocol_and_clock_selected':False}


def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def source_checks():
    pins=json.loads(SOURCES.read_text());checked={}
    for path,expected in pins['upstream_sha256'].items():
        actual=digest(ROOT/path)
        if actual!=expected:raise ValueError('upstream source changed: '+path)
        checked[path]=actual
    for path in pins['local_inputs']:checked[path]=digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))]=digest(SOURCES)
    return checked


def run():
    if not __debug__:raise RuntimeError('Optimized Python is refused')
    return {'certificate_type':'YM54_NATIVE_RESPONSE_CURVATURE_TO_COUNTED_PROTOCOL',
        'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
        'native_embedding':native_embedding_controls(),'response_geometry':response_controls(),
        'counted_protocol':counted_protocol_controls(),'covariance':covariance_controls(),
        'refinement':refinement_controls(),'thermo_fixture':fixture_controls(),
        'scale_and_failures':scale_and_failure_controls(),
        'claim_status':'DECLARED_RESPONSE_PROTOCOL_AND_CURVATURE_BUDGET_GAP_BRIDGE__'
                       'PHYSICAL_SELECTION_AND_ANISOTROPIC_JOINT_LIMIT_OPEN',
        'evidence_scope':{'general_proof':'written, not mechanically formalized or expert-certified',
            'protocol':'two shape directions; equal signed counts; frozen response jets',
            'derived':'C from actual turn records; beta >= 4*f^2/tau; inherited YM53 gap',
            'external_sources':'lineage/context metadata; external proofs not recertified',
            'open':'physical response, directions, weights, clock, coupling and state; state-dependent '
                   'diffusion; anisotropic joint/volume limit; row closure; 4D/AF/Clay/QG'}}


def check(cert,result=RESULT,pin=PIN):
    sha=chain.old.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM54 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM54 PASS',sha)
    print(json.dumps({'response_jets':cert['response_geometry']['response_jets'],
        'fixture_beta':cert['thermo_fixture']['beta_from_curvature'],
        'fixture_gamma_lower':cert['thermo_fixture']['gamma_lower'],
        'physical_selection':'OPEN'},sort_keys=True))


if __name__=='__main__':main()
