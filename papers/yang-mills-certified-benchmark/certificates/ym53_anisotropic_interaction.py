"""YM53 exact controls for the written anisotropic interacting chain proof.

Python 3.12 only. Finite fixtures are not compact-carrier discretizations.
Default/--check is read-only; --write touches only YM53 result and pin.
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
import ym52_energy_observer as z
import ym45_temporal_blocks as old
import ym44_time_refinement as refine

a=z.a
p=z.p
Iv=old.Iv
exact=z.q.exact
RESULT=HERE/'YM53_RESULT.json'
PIN=HERE/'EXPECTED_YM53.sha256'
SOURCES=HERE/'YM53_SOURCE_PINS.json'
CELLS=((F(1,8192),8,F(3,2)),(F(1,4096),6,F(5,4)),
       (F(1,3072),5,F(11,10)),(F(1,2560),5,F(33,32)))


def integer(n,minimum=0):
    if isinstance(n,bool) or not isinstance(n,int) or n<minimum:
        raise ValueError('integer outside declared domain')
    return n


def geometric_tail(q,N=0):
    q=exact(q);N=integer(N)
    if not 0<=q<1:raise ValueError('geometric ratio must lie in [0,1)')
    return q**(N+1)*((N+2)**2/(1-q)+2*(N+2)*q/(1-q)**2
                       +q*(1+q)/(1-q)**3)


def heat_controls():
    factorial_floor=sum((F(9,2)**k/factorial(k) for k in range(9)),F(0))
    assert factorial_floor>84
    q=F(1,84);tail=geometric_tail(q)
    assert tail==(1+q)/(1-q)**3-1<F(1,20)
    assert F(19,21)**2>=old.EPS
    checks=0
    for q in (F(0),F(1,1000),F(1,84),F(1,2),F(9,10)):
        for N in (0,1,4,12,32):
            finite=sum(((n+1)**2*q**n for n in range(1,N+1)),F(0))
            assert finite+geometric_tail(q,N)==(1+q)/(1-q)**3-1
            assert geometric_tail(q,N)>=sum(((n+1)**2*q**n
                       for n in range(N+1,N+21)),F(0))
            checks+=1
    return {'factorial_order':8,'exp_9_over_2_lower':str(factorial_floor),
            'all_content_tail_upper':str(tail),'tail_identity_checks':checks,
            'coarse_time':'9/beta','kernel_floor':'19/20',
            'kernel_ceiling':'21/20','normalized_overlap':str(old.EPS)}


def spin_controls():
    checks=odd=rank_two=0
    for n in range(1,13):
        squares,G=z.spin_squares(n);j=F(n,2);I=a.identity(n+1)
        for vals in z.DIAGONALS:
            L=a.scale(I,0)
            for v,sq in zip(vals,squares):L=a.add(L,a.scale(sq,v))
            bound=vals[0]*j*j+vals[1]*j
            assert a.psd(a.multiply(G,a.add(L,a.scale(I,-bound))))
            assert bound>=vals[1]*j
            checks+=1;odd+=n%2;rank_two+=vals[0]==0<vals[1]
    return {'finite_spin_inequalities':checks,'half_integer_cases':odd,
            'rank_two_cases':rank_two,'degrees':[1,12],
            'infinite_extension':'YM53-T1 Casimir/operator-norm argument'}


def block_gate(beta,theta,J,R):
    beta=exact(beta);theta=exact(theta);R=exact(R);integer(J,1)
    if beta<=0 or R<=1:raise ValueError('positive beta and R>1 required')
    alpha=F(2)/(old.EPS*J)+240*abs(theta)/beta*J
    if R*alpha>=1:return None
    gamma=old.log_iv(Iv(R),old.LOG_TERMS)*Iv(beta)/Iv(10*J)
    assert 0<gamma.lo<=gamma.hi<beta/2
    return {'alpha':alpha,'margin':1-R*alpha,'gamma':gamma}


def skeleton(beta,h):
    beta=exact(beta);h=exact(h)
    if beta<=0 or not 0<h<=1/beta:raise ValueError('require 0<h<=1/beta')
    ratio=9/(beta*h)
    return -(-ratio.numerator//ratio.denominator)


def parameter_controls():
    cells=[];fine=identities=0
    for u,J,R in CELLS:
        gate=block_gate(1,u,J,R);assert gate is not None
        for beta in (F(1,100),F(1),F(3,2),F(7)):
            for step in (F(1),F(2,3),F(1,7),F(17,256),F(1,65536)):
                h=step/beta;r=skeleton(beta,h);ell=J*r
                assert 9/beta<=h*r<10/beta
                alpha=2*r/(old.EPS*ell)+24*h*u*beta*ell
                assert alpha<=gate['alpha'] and R*alpha<1
                fine+=1
        cells.append({'theta_over_beta':str(u),'J':J,'R':str(R),
            'alpha_upper':str(gate['alpha']),'weighted_margin':str(gate['margin']),
            'gamma_over_beta_lower':old.decimals(gate['gamma'].lo,False),
            'gamma_over_beta_upper':old.decimals(gate['gamma'].hi)})
    for u in (F(0),F(1,4096),F(1,2400),F(1,64),F(1,16)):
        for J in range(1,31):
            alpha=F(5,2*J)+240*u*J
            assert J*(alpha-1)==F((J-5)**2,10)+(240*u-F(1,10))*J*J
            if u>=F(1,2400):assert alpha>=1
            identities+=1
    for u in (F(0),F(1,100000),F(1,3000),F(999,2400000)):
        alpha=F(1,2)+1200*u;R=(1+1/alpha)/2
        assert R*alpha==(1+alpha)/2<1
        assert block_gate(1,u,5,R) is not None
    sample=block_gate(F(3,2),F(1,4096),8,F(4,3))
    assert sample['alpha']==F(5,8) and sample['margin']==F(1,6)
    return {'window':'abs(theta) < beta/2400','cells':cells,
            'fine_step_cases':fine,'endpoint_identity_checks':identities,
            'rank_two_example':{'tensor':['0','3/2','3/2'],'beta':'3/2',
                'theta':'1/4096','J':8,'R':'4/3','alpha':'5/8','margin':'1/6',
                'gamma_lower':old.decimals(sample['gamma'].lo,False),
                'gamma_upper':old.decimals(sample['gamma'].hi)}}


def path_kernel():
    return [[F(3,4),F(1,4),F(0),F(0)],
            [F(1,4),F(1,2),F(1,4),F(0)],
            [F(0),F(1,4),F(1,2),F(1,4)],
            [F(0),F(0),F(1,4),F(3,4)]]


def short_bridge(left,right):
    K=path_kernel()
    return old.normalized([K[left][i]*K[i][right] for i in range(4)])


def zero_bridge_controls():
    K=path_kernel();K3=a.multiply(a.multiply(K,K),K)
    assert K==a.transpose(K) and a.psd(K)
    assert all(sum(row)==1 for row in K)
    assert all(x>0 for row in K3 for x in row)
    rejected=False
    try:short_bridge(0,3)
    except ValueError:rejected=True
    assert rejected
    left=short_bridge(0,1);right=short_bridge(2,3)
    joint=old.coupling(left,right)
    assert [sum(row) for row in joint]==left
    assert [sum(row[j] for row in joint) for j in range(4)]==right
    cost=sum(joint[i][j] for i,j in itertools.product(range(4),repeat=2) if i!=j)
    assert cost<=1
    # Positive spatial tilts cannot change the support of an admissible law.
    tilted=old.normalized([v*w for v,w in zip(left,(1,2,3,4))])
    assert [v>0 for v in tilted]==[v>0 for v in left]
    return {'zero_normalizer_rejected':rejected,'coarse_power_positive':3,
            'admissible_pairs':[[0,1],[2,3]],'invalid_intermediate':[0,3],
            'grouped_endpoint_hamming_cost':str(cost),
            'positive_tilt_preserves_support':True}


COORD8=tuple({tuple(int(i==j) for i in range(8)):F(1)} for j in range(8))


def site_derivative(f,site,axis):
    # Reuse the native four-coordinate generator, lifted to a product.
    out={}
    for j in range(4):
        native=z.y.native_generator(p.COORD[j],axis+1)
        field={((0,)*4+m if site else m+(0,)*4):v for m,v in native.items()}
        out=p.add(out,p.mul(field,p.partial(f,4*site+j)))
    return out


def site_lap(f,site,C):
    return p.scale(p.add(*(p.scale(site_derivative(site_derivative(f,site,j),site,i),
                    C[i][j]) for i,j in itertools.product(range(3),repeat=2))),-1)


def bilinear_energy(df,C):
    return p.add(*(p.scale(p.mul(df[i],df[j]),C[i][j])
                   for i,j in itertools.product(range(3),repeat=2)))


def bond_controls():
    v=p.add(*(p.mul(COORD8[j],COORD8[j+4]) for j in range(4)))
    radii=[p.add(*(p.mul(COORD8[4*i+j],COORD8[4*i+j]) for j in range(4)))
           for i in range(2)]
    R=p.mul(*radii);v2=p.mul(v,v)
    derivatives=[[site_derivative(v,i,j) for j in range(3)] for i in range(2)]
    for ds in derivatives:
        assert p.add(*(p.mul(d,d) for d in ds))==p.scale(p.add(R,p.scale(v2,-1)),F(1,4))
    tensors=[a.diag(vals) for vals in z.DIAGONALS[2:]]
    tensors += [z.q.rotate_tensor(a.diag((0,1,2)),z.ROT2)]
    checks=0
    for C,D in zip(tensors,tensors[1:]+tensors[:1]):
        traces=[sum(T[i][i] for i in range(3)) for T in (C,D)]
        lam=sum(traces)/4
        assert p.add(site_lap(v,0,C),site_lap(v,1,D))==p.scale(v,lam)
        E=p.add(*(bilinear_energy(ds,T) for ds,T in zip(derivatives,(C,D))))
        complements=[a.add(a.scale(z.q.I3,t),T,-1) for T,t in zip((C,D),traces)]
        assert all(a.psd(T) for T in complements)
        sos=p.add(p.scale(v2,lam),*(bilinear_energy(ds,T)
                         for ds,T in zip(derivatives,complements)))
        assert p.add(p.scale(R,lam),p.scale(E,-1))==sos
        checks+=1
    return {'site_gradient_identities':2,'inhomogeneous_tensor_pairs':checks,
            'bond_eigen_and_positive_complement_identities':checks,
            'coordinate_count':8,'constant':'lambda=(tr(C)+tr(D))/4'}


def refinement_budget(width,theta,trace_max,t,mesh):
    integer(width,1)
    theta=exact(theta);trace_max=exact(trace_max);t=exact(t);mesh=exact(mesh)
    if trace_max<=0 or t<0 or not 0<=mesh<=t:
        raise ValueError('positive trace and 0<=mesh<=time required')
    b=abs(theta)*(width-1)
    if not b or not t or not mesh:return F(0)
    return b*t*refine.ex(b*t).hi*refine.iv_sqrt(Iv(3*trace_max*mesh)).hi


def refinement_controls():
    moments=0;rows=[]
    for lam in (F(1,100),F(3,2),F(3),F(9)):
        for s in (F(0),F(1,65536),F(1,16),F(1),F(4)):
            q=refine.ex(-lam*s)
            variance=Iv(1)-q*q;drift=(Iv(1)-q)*(Iv(1)-q)
            direct=Iv(2)*(Iv(1)-q)
            assert not (variance+drift).separated_from(direct)
            # Scalar analytic bound is checked with outward enclosures.
            assert direct.lo<=2*lam*s and direct.hi>=0
            if s:assert direct.hi<=2*lam*s
            moments+=1
    theta=F(1,4096);gate=block_gate(F(3,2),theta,8,F(4,3))
    for width in (1,2,4,9):
        for t in (F(1,4),F(1),F(4)):
            b=theta*(width-1);floor=refine.ex(-b*t).lo
            ceiling=refine.ex(-gate['gamma'].lo*t).hi
            n=16;steps=0;previous=None
            while True:
                error=refinement_budget(width,theta,3,t,t/n)
                if previous is not None:assert error<previous or not error
                ratio=refine.gap_gate(floor,ceiling,error)
                steps+=1
                if ratio is not None:break
                previous=error;n*=4
                assert n<=2**26
            assert ratio<1
            rows.append({'width':width,'time':str(t),'steps':n,
                'error_upper':old.decimals(error),'normalized_gap_ceiling_upper':old.decimals(ratio),
                'mesh_checks':steps})
    return {'moment_enclosures':moments,'gap_transport_cells':rows,
            'error':'b*time*exp(b*time)*sqrt(3*trace_max*mesh)',
            'uniform_in_unbounded_trace':False}


def ground_energy_controls():
    # Independent four-state identity, not a replacement for native calculus.
    L=a.scale(a.identity(4),0)
    weights=(F(1),F(2),F(3));h=list(map(F,(1,2,3,2)))
    for i,w in enumerate(weights):
        L[i][i]+=w;L[i+1][i+1]+=w;L[i][i+1]-=w;L[i+1][i]-=w
    Lh=p.apply(L,h);B=[v/w for v,w in zip(Lh,h)]
    H=a.add(L,a.diag(B),-1);weighted=a.multiply(a.multiply(a.diag(h),H),a.diag(h))
    edge=a.scale(a.identity(4),0)
    for i,w in enumerate(weights):
        v=w*h[i]*h[i+1]
        edge[i][i]+=v;edge[i+1][i+1]+=v;edge[i][i+1]-=v;edge[i+1][i]-=v
    assert weighted==edge and a.psd(H) and p.apply(H,h)==[0]*4
    assert a.y37.inertia(H)==(3,0,1)
    norm=sum(x*x for x in h);checks=0;free_mismatch=False
    for f in ([1,1,1,1],[1,0,0,0],[1,-1,2,-2],[0,1,1,0]):
        fh=[x*y for x,y in zip(f,h)]
        energy=p.dot(fh,p.apply(H,fh))/norm
        direct=sum(w*h[i]*h[i+1]*(f[i]-f[i+1])**2
                   for i,w in enumerate(weights))/norm
        free=p.dot(f,p.apply(L,f))/4
        assert energy==direct>=0
        free_mismatch|=energy!=free;checks+=1
    assert free_mismatch
    return {'weighted_form_matrix_entries':16,'source_checks':checks,
            'ground_kernel_dimension':1,'interaction_potential':list(map(str,B)),
            'omitting_ground_weight_rejected':free_mismatch}


def scale_controls():
    checks=0
    for beta in (F(1,100),F(1),F(3,2)):
        theta=beta/4096;h=F(1,7)/beta
        base=block_gate(beta,theta,6,F(5,4))
        for k in (F(1,10),F(2),F(7,3)):
            scaled=block_gate(k*beta,k*theta,6,F(5,4))
            assert scaled['alpha']==base['alpha'] and scaled['margin']==base['margin']
            assert scaled['gamma'].lo==k*base['gamma'].lo
            assert scaled['gamma'].hi==k*base['gamma'].hi
            assert skeleton(beta,h)==skeleton(k*beta,h/k)
            # The dimensionless arguments of the norm budget stay identical.
            b=3*theta;t=F(2);M=3*beta;mesh=t/100
            assert (k*b)*(t/k)==b*t and 3*(k*M)*(mesh/k)==3*M*mesh
            checks+=1
    return {'clock_covariant_cells':checks,'physical_clock_selected':False}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    return {'certificate_type':'YM53_ANISOTROPIC_INTERACTING_CHAIN_GAP',
        'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
        'all_content_heat':heat_controls(),'spin_controls':spin_controls(),
        'parameters':parameter_controls(),'admissible_bridges':zero_bridge_controls(),
        'bond_energy':bond_controls(),'norm_refinement':refinement_controls(),
        'ground_source_energy':ground_energy_controls(),'clock':scale_controls(),
        'reused_YM45_controls':{'bridges':old.bridge_controls(),'tilts':old.tilt_controls(),
            'incidence':old.incidence_controls(),'joint_updates':old.joint_update_controls()},
        'claim_status':'DECLARED_ANISOTROPIC_CHAIN_UNIFORM_GAP_AND_FIXED_WIDTH_TIME_LIMIT__'
                       'ANISOTROPIC_JOINT_LIMIT_AND_PHYSICAL_SELECTION_OPEN',
        'evidence_scope':{'general_proof':'written; not mechanically formalized or expert-certified',
            'finite_fixtures':'algebraic controls, not compact-carrier discretizations',
            'carrier':'native quaternion coefficient completion, YM50 representation bridge',
            'window':'sufficient, not optimized or a phase boundary',
            'external_sources':'method attribution only; no external theorem recertified',
            'open':'physical protocol/state/action/clock; anisotropic volume/joint limit; actual '
                   'row closure; four-dimensional QFT, asymptotic freedom and Clay'}}


def check(cert,result=RESULT,pin=PIN):
    sha=old.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM53 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM53 PASS',sha)
    print(json.dumps({'window':cert['parameters']['window'],
        'rank_two_example':cert['parameters']['rank_two_example'],
        'all_content_proof':'written, finite controls reproduced','physical_selection':'OPEN'},sort_keys=True))


if __name__=='__main__':main()
