"""YM57 exact neighbour-turn controls. Python 3.12 only; no new one-site engine.

--check is read-only; --write updates only YM57 evidence.
Finite controls support, but do not replace, the all-content written proof.
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
import ym56_lambda_curvature_gap as previous

a,p,q,y=previous.a,previous.p,previous.q,previous.y
exact=previous.exact
RESULT=HERE/'YM57_RESULT.json'
PIN=HERE/'EXPECTED_YM57.sha256'
SOURCES=HERE/'YM57_SOURCE_PINS.json'


def width(n, minimum=1):
    if isinstance(n,bool) or not isinstance(n,int) or n<minimum:
        raise ValueError('Invalid chain width')
    return n


def lower_bound(n,k):
    width(n,2)
    k=exact(k)
    if k<=0:
        raise ValueError('Positive edge rate required')
    return min(F(1,630),k/420)


def lift(poly,site,n):
    width(n)
    if isinstance(site,bool) or not isinstance(site,int) or not 0<=site<n:
        raise ValueError('Invalid site')
    return {(0,)*(4*site)+m+(0,)*(4*(n-site-1)):c for m,c in poly.items()}


def derivative(poly,site,axis):
    out={}
    for m,c in poly.items():
        local=y.native_generator({m[4*site:4*site+4]:c},axis)
        for term,v in local.items():
            full=m[:4*site]+term+m[4*site+4:]
            out[full]=out.get(full,F(0))+v
    return p.clean(out)


def bond(poly,site):
    return p.add(derivative(poly,site,2),p.scale(derivative(poly,site+1,2),-1))


def operator(poly,n,k,lam=F(0)):
    width(n)
    k,lam=exact(k),exact(lam)
    if k<0:
        raise ValueError('Nonnegative edge rate required')
    terms=[]
    for i in range(n):
        terms.append(derivative(derivative(poly,i,1),i,1))
        terms.append(p.scale(derivative(derivative(poly,i,2),i,2),lam*lam))
    terms.extend(p.scale(bond(bond(poly,i),i),k) for i in range(n-1))
    return p.scale(p.add(*terms),F(-1,2))


def phi(poly,n):
    return sum((c*y.prod(y.phi({m[4*i:4*i+4]:F(1)}) for i in range(n))
                for m,c in poly.items()),F(0))


def variance(poly,n):
    return phi(p.mul(poly,poly),n)-phi(poly,n)**2


def paired_product(x,z):
    return tuple(y.qmul(v,w) for v,w in zip(x,z))


def echo(turn,half_turn=None):
    h=half_turn if half_turn is not None else (F(0),F(1),F(0),F(0))
    b=(turn,y.qdagger(turn))
    u=(y.ONE,h)
    return paired_product(paired_product(paired_product(b,u),b),(y.ONE,y.qdagger(h)))


def echo_controls():
    turns=[y.ONE,tuple(-v for v in y.ONE),(F(0),F(0),F(1),F(0)),
           (F(3,5),F(0),F(4,5),F(0)),(F(5,13),F(0),F(-12,13),F(0))]
    for turn in turns:
        assert y.qmul(turn,y.qdagger(turn))==y.ONE
        assert echo(turn)==(y.qmul(turn,turn),y.ONE)
    assert echo(turns[3],y.ONE)!=(y.qmul(turns[3],turns[3]),y.ONE)
    # Cauchy budget in the Euler/echo estimate and rational pi bound.
    assert 4**2+1**2+2**2==21
    assert F(22,7)**2<10
    return dict(exact_quaternion_words=len(turns),wrong_half_turn_refused=True,
                cauchy_coefficient=21,pi_upper_used='22/7',
                all_angles_proof='YM57-T3 conjugation and commuting second-axis flows')


def obstruction_controls():
    n=2
    w=lift(previous.witness(1),0,n)
    coords=[[lift(f,i,n) for f in p.COORD] for i in range(n)]
    V=p.add(*(p.mul(coords[0][j],coords[1][j]) for j in range(4)))
    checks=0
    for f in (lift(p.ONE,0,n),coords[0][0],coords[1][1],w,p.mul(*[c[0] for c in coords])):
        for theta in (F(-2),F(0),F(1,3)):
            def original(g):
                return p.add(operator(g,n,0),p.scale(p.mul(V,g),-theta))
            assert original(p.mul(w,f))==p.mul(w,original(f))
            checks+=1
    assert not derivative(w,0,1)
    assert bond(w,0)
    dw=previous.y.native_generator(previous.witness(1),2)
    # D2(w) consists of two signed disjoint-coordinate products.
    assert len(dw)==2 and all(abs(v)==1 for v in dw.values())
    support=[{i for i,e in enumerate(m) if e} for m in dw]
    assert len(support[0])==len(support[1])==2 and not support[0]&support[1]
    assert phi(w,n)==0 and variance(w,n)==F(1,12)
    return dict(potential_commutation_checks=checks,rank_one_invariant=True,
                new_bond_acts_nontrivially=True,weighted_variance_must_not_be_replaced=True)


@lru_cache(maxsize=None)
def block(degrees):
    mons=[tuple(v for local in row for v in local)
          for row in itertools.product(*(y.basis(d) for d in degrees))]
    n=len(degrees)
    polys=[{m:F(1)} for m in mons]
    G=[[phi(p.mul(f,g),n) for g in polys] for f in polys]
    means=[phi(f,n) for f in polys]
    centered=[[G[i][j]-means[i]*means[j] for j in range(len(mons))]
              for i in range(len(mons))]
    return mons,polys,G,centered


def matrix_controls():
    count=positive=0
    degrees_list=((1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(1,2),(1,0,1))
    for degrees in degrees_list:
        n=len(degrees)
        mons,polys,G,centered=block(degrees)
        for k in (F(1,10),F(1),F(3)):
            columns=[operator(f,n,k) for f in polys]
            H=[[col.get(m,F(0)) for col in columns] for m in mons]
            form=a.multiply(G,H)
            assert form==a.transpose(form)
            assert a.psd(a.add(form,a.scale(centered,-lower_bound(n,k))))
            # Independent factorization by the native derivative matrices.
            factored=a.scale(G,0)
            for site in range(n):
                images=[derivative(f,site,1) for f in polys]
                factored=a.add(factored,[[phi(p.mul(v,w),n)/2 for w in images] for v in images])
            for site in range(n-1):
                images=[bond(f,site) for f in polys]
                factored=a.add(factored,[[k*phi(p.mul(v,w),n)/2 for w in images] for v in images])
            assert factored==form
            count+=1;positive+=1
    return dict(tensor_blocks=[list(d) for d in degrees_list],
                rational_gap_inequalities=count,independent_square_factorizations=positive,
                max_block_dimension=max(len(block(d)[0]) for d in degrees_list),
                all_content_and_width_proof='YM57-T3--T4 finite words and variance tensorization')


def graph_controls():
    rows=[]
    for n in (2,3,4,5,16,64,257):
        neighbour=[i+1 if i<n-1 else n-2 for i in range(n)]
        indegrees=[neighbour.count(i) for i in range(n)]
        edge_counts=[0]*(n-1)
        for i,j in enumerate(neighbour):edge_counts[min(i,j)]+=1
        assert max(indegrees)<=2 and max(edge_counts)<=2
        assert max(1+c for c in indegrees)<=3
        for k in (F(1,100),F(1),F(5)):
            g=lower_bound(n,k)
            # Replace pi^2 by its strict upper budget 10 in (8).
            assert g*F(630,2)<=F(1,2)
            assert g*210<=k/2
        rows.append(dict(width=n,max_indegree=max(indegrees),max_edge_uses=max(edge_counts)))
    return dict(rows=rows,uniform_rational_rate_at_unit_edge='1/630')


def record_controls():
    checks=0
    for n in (2,3,5):
        for k in (F(1,10),F(1),F(3)):
            W=n+k*(n-1)
            directions=[]
            for i in range(n):
                v=[F(0)]*(2*n);v[2*i]=1
                directions.append((F(1),v))
            for i in range(n-1):
                v=[F(0)]*(2*n);v[2*i+1]=1;v[2*(i+1)+1]=-1
                directions.append((k,v))
            records=[(weight/(2*W),[sign*x for x in v])
                     for weight,v in directions for sign in (-1,1)]
            assert sum(weight for weight,v in records)==1
            assert all(sum(weight*v[j] for weight,v in records)==0 for j in range(2*n))
            covariance=[[sum(W*weight*v[i]*v[j] for weight,v in records)
                         for j in range(2*n)] for i in range(2*n)]
            expected=a.scale(a.identity(2*n),0)
            for i in range(n):expected[2*i][2*i]=1
            for i in range(n-1):
                u,v=2*i+1,2*(i+1)+1
                expected[u][u]+=k;expected[v][v]+=k
                expected[u][v]-=k;expected[v][u]-=k
            assert covariance==expected and a.psd(covariance)
            checks+=1
    return dict(exact_signed_record_second_moments=checks,
                spatial_covariance_singular=True,
                singular_covariance_is_not_a_gaplessness_test=True)


def symmetry_and_clock_controls():
    identities=0
    rows=[]
    for n in (2,3,6):
        w=lift(previous.witness(1),0,n)
        for k in (F(1,10),F(1),F(2)):
            E=phi(p.mul(w,operator(w,n,k)),n)
            quotient=E/variance(w,n)
            assert quotient==k/2
            assert operator(w,n,k)==p.scale(w,k/2)
            W=n+k*(n-1)
            rows.append(dict(width=n,k=k,extensive_endpoint_rayleigh=quotient,
                             normalized_endpoint_rayleigh=quotient/W))
            for lam in (F(0),F(1,2),F(-2)):
                assert phi(p.mul(w,operator(w,n,k,lam)),n)/variance(w,n)==(k+lam*lam)/2
        for site in range(n):
            f=p.add(lift(p.COORD[0],site,n),lift(previous.witness(2),site,n))
            def flip(g):
                return {m:c*(-1)**sum(m[4*i+2]+m[4*i+3] for i in range(n)) for m,c in g.items()}
            assert flip(operator(f,n,1))==operator(flip(f),n,1)
            identities+=1
    assert not operator(lift(previous.witness(1),0,1),1,1)
    assert not operator(lift(previous.witness(1),0,3),3,0)
    # Fixed-rate whole-chain update cannot retain a width-uniform floor.
    normalized=[]
    for n in (2,10,100,1000):normalized.append(F(1,2*(2*n-1)))
    assert normalized==sorted(normalized,reverse=True)
    return dict(local_symmetry_checks=identities,endpoint_witnesses=rows,
                fixed_rate_clock_upper_bounds=normalized,
                scalar_interaction_imported=False,infinite_volume_constructed=False)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_checks():
    manifest=json.loads(SOURCES.read_text());checked={}
    for path,expected in manifest['upstream_sha256'].items():
        actual=digest(ROOT/path)
        if actual!=expected:raise ValueError('upstream source changed: '+path)
        checked[path]=actual
    for path in manifest['local_inputs']:checked[path]=digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))]=digest(SOURCES)
    return checked


def run():
    previous.runtime_check()
    return previous.encode(dict(certificate_type='YM57_NEIGHBOUR_ECHO_UNIFORM_COMPACT_GAP',
        verdict='PASS',runtime_policy='Python 3.12 only',input_sha256=source_checks(),
        echo=echo_controls(),obstruction=obstruction_controls(),matrices=matrix_controls(),
        graph=graph_controls(),records=record_controls(),symmetry_clock=symmetry_and_clock_controls(),
        claim_status='SELECTED_CORRELATED_KINETIC_FINITE_CHAIN_UNIFORM_GAP__4D_CLAY_OPEN',
        proof_scope='Written all-content proof; finite exact controls are not formal verification',
        open='Correlated infinite-volume histories; scalar interaction; physical gauge action and clock; spatial continuum'))


def check(cert,result=RESULT,pin=PIN):
    sha=previous.native.chain.old.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM57 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=previous.native.chain.old.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM57 PASS',sha)
    print(json.dumps(dict(runtime=sys.version.split()[0],rate_at_unit_edge='1/630',
                          finite_tensor_controls=cert['matrices']['rational_gap_inequalities'])))


if __name__=='__main__':main()
