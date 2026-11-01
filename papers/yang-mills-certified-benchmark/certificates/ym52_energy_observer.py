"""YM52 exact energy/record/observer controls. Python 3.12 only.

Written all-content proofs are distinct from these finite controls.
Default and --check are read-only; --write touches YM52 evidence only.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from math import comb
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE))
import ym51_heat_selection as q

y=q.y
p=q.p
a=q.a
RESULT=HERE/'YM52_RESULT.json'
PIN=HERE/'EXPECTED_YM52.sha256'
SOURCES=HERE/'YM52_SOURCE_PINS.json'
EVEN=p.add(p.mul(p.COORD[1],p.COORD[1]),p.mul(p.COORD[2],p.COORD[2]),
           p.scale(y.RADIUS,F(-1,2)))
DIAGONALS=tuple(tuple(map(F,row)) for row in
               ((0,0,0),(0,0,3),(0,1,2),(0,F(3,2),F(3,2)),
                (F(1,100),F(1,100),F(149,50)),(F(1,4),F(1,2),F(9,4)),
                (F(1,2),1,F(3,2)),(1,1,1),(1,2,9)))
ROT=q.rotation((F(1,2),)*4)
ROT2=q.rotation((F(3,5),F(4,5),F(0),F(0)))


def gamma(f,g,C):
    C=q.tensor(C)
    df=[y.native_generator(f,i+1) for i in range(3)]
    dg=[y.native_generator(g,i+1) for i in range(3)]
    return p.add(*(p.scale(p.mul(df[i],dg[j]),C[i][j])
                   for i,j in itertools.product(range(3),repeat=2)))


def energy(f,C):
    return y.phi(gamma(f,f,C))


def variance(f):
    return y.phi(p.mul(f,f))-y.phi(f)**2


def antipodal(f):
    return {m:c*(-1)**sum(m) for m,c in f.items()}


def parity(f,even=True):
    return p.scale(p.add(f,p.scale(antipodal(f),1 if even else -1)),F(1,2))


def cross(u,v):
    return tuple(u[(i+1)%3]*v[(i+2)%3]-u[(i+2)%3]*v[(i+1)%3]
                 for i in range(3))


def cofactor(C):
    C=q.tensor(C)
    out=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rows=[r for r in range(3) if r!=i]
            cols=[s for s in range(3) if s!=j]
            row.append((-1)**(i+j)*(C[rows[0]][cols[0]]*C[rows[1]][cols[1]]
                                    -C[rows[0]][cols[1]]*C[rows[1]][cols[0]]))
        out.append(row)
    assert a.psd(out)
    return out


def bracket_tensor(records):
    records=q.protocol(records)
    out=a.scale(q.I3,0)
    for (w,u),(z,v) in itertools.combinations(records,2):
        c=cross(u,v)
        out=a.add(out,[[w*z*c[i]*c[j] for j in range(3)] for i in range(3)])
    return out


def rates(eigenvalues):
    vals=tuple(q.exact(v) for v in eigenvalues)
    if len(vals)!=3 or vals!=tuple(sorted(vals)) or vals[0]<0:
        raise ValueError('three ordered nonnegative exact eigenvalues required')
    x,z,w=vals
    odd=(x+z+w)/4
    even=x+z
    return min(odd,even),even,odd


def energy_controls():
    basis=[p.ONE,*p.COORD,EVEN,
           p.add(p.COORD[0],p.mul(p.COORD[1],p.COORD[3]))]
    tensors=[a.diag(v) for v in DIAGONALS[1:]]
    tensors.append(q.rotate_tensor(a.diag(DIAGONALS[2]),ROT2))
    products=balances=0
    for C in tensors:
        for f,g in itertools.product(basis,repeat=2):
            defect=p.add(q.lap(p.mul(f,g),C),p.scale(p.mul(f,q.lap(g,C)),-1),
                         p.scale(p.mul(g,q.lap(f,C)),-1),p.scale(gamma(f,g,C),2))
            assert not defect
            assert y.phi(p.mul(f,q.lap(g,C)))==y.phi(gamma(f,g,C))
            products+=1
        for f in basis:
            h=q.lap(f,C)
            assert energy(f,C)>=0
            assert -2*y.phi(p.mul(f,h))==-2*energy(f,C)
            assert -2*y.phi(gamma(f,h,C))==-2*y.phi(p.mul(h,h))
            balances+=1
    f=p.COORD[0];C=q.I3
    assert y.phi(p.mul(f,f))==F(1,4)
    assert energy(f,C)==F(3,16)
    assert -2*energy(f,C)==F(-3,8)
    return {'leibniz_and_integration_by_parts':products,
            'norm_and_energy_derivatives':balances,
            'missing_factor_two_rejected':True}


def average_step(f,words):
    return p.scale(p.add(*(y.substitute(f,y.qmatrix(w)) for w in words)),
                   F(1,len(words)))


def record_controls():
    words=tuple(t for t in y.turns() if not t[3])
    inputs=[p.COORD[0],p.COORD[1],EVEN,
            p.add(p.COORD[0],p.scale(EVEN,2))]
    checks=0;lost=[]
    for f in inputs:
        current=f;counts=(y.ONE,)
        for n in range(1,4):
            counts=tuple(y.qmul(w,z) for z in counts for w in words)
            current=average_step(current,words)
            direct=average_step(f,counts)
            assert current==direct
            second=average_step(p.mul(f,f),counts)
            noise=p.add(second,p.scale(p.mul(current,current),-1))
            explicit=p.scale(p.add(*(p.mul(
                p.add(y.substitute(f,y.qmatrix(w)),p.scale(current,-1)),
                p.add(y.substitute(f,y.qmatrix(w)),p.scale(current,-1)))
                for w in counts)),F(1,len(counts)))
            assert noise==explicit
            assert y.phi(second)==y.phi(p.mul(f,f))
            assert y.phi(noise)>=0
            assert y.phi(p.mul(f,f))==y.phi(p.mul(current,current))+y.phi(noise)
            at_one=[y.evaluate(f,y.qmul(w,y.ONE)) for w in counts]
            raw=sum(z*z for z in at_one)/len(counts)-(sum(at_one)/len(counts))**2
            assert raw==y.evaluate(noise,y.ONE)
            checks+=1
            lost.append(str(y.phi(noise)))
    assert any(F(v)>0 for v in lost)
    return {'independent_word_variance_checks':checks,
            'largest_literal_word_bank':len(counts),'integrated_noise':lost}


def bracket_controls():
    protocols=[q.eta_protocol(F(0)),q.eta_protocol(F(1,100)),
               q.protocol([(F(1,2),(2,0,0)),(F(1,2),(0,2,0))]),
               q.protocol([(F(1,4),v) for v in
                           ((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))]),
               q.protocol([(F(1,3),v) for v in ((1,2,0),(0,1,2),(2,0,1))])]
    pair_checks=0;identities=0
    for records in protocols:
        C=q.moments(records)[0];B=bracket_tensor(records)
        assert B==cofactor(C)
        split=tuple((w/2,v) for w,v in records for _ in range(2))
        assert bracket_tensor(split)==B
        assert cofactor(q.rotate_tensor(C,ROT2))==q.rotate_tensor(B,ROT2)
        for f in [*p.COORD,EVEN]:
            rhs={}
            for (w,u),(z,v) in itertools.combinations(records,2):
                comm=p.add(q.derivative(q.derivative(f,v),u),
                           p.scale(q.derivative(q.derivative(f,u),v),-1))
                assert comm==q.derivative(f,cross(u,v))
                rhs=p.add(rhs,p.scale(p.mul(comm,comm),w*z))
                pair_checks+=1
            assert rhs==gamma(f,f,B)
            identities+=1
    ranks=[]
    for vals in DIAGONALS:
        C=a.diag(vals);B=cofactor(C)
        # Rank of this dimensionless diagnostic only; not a physical sum of units.
        rank=a.y37.inertia(a.add(C,B))[0]
        active=sum(v>0 for v in vals)
        assert rank==(3 if active>=2 else active)
        ranks.append([active,rank])
    return {'cofactor_factorizations':len(protocols),'bracket_actions':pair_checks,
            'bracket_energy_identities':identities,'rank_controls':ranks}


@lru_cache(maxsize=None)
def spin_squares(n):
    if isinstance(n,bool) or not isinstance(n,int) or n<0:
        raise ValueError('nonnegative integer tensor degree required')
    size=n+1
    plus=a.scale(a.identity(size),0);minus=a.scale(a.identity(size),0)
    for k in range(size):
        if k:plus[k-1][k]=F(k)
        if k<n:minus[k+1][k]=F(n-k)
    G=a.diag([F(1,comb(n,k)) for k in range(size)])
    assert a.multiply(G,plus)==a.multiply(a.transpose(minus),G)
    J1=a.scale(a.add(plus,minus),F(1,2))
    iJ2=a.scale(a.add(plus,minus,-1),F(1,2))
    squares=(a.multiply(J1,J1),a.scale(a.multiply(iJ2,iJ2),-1),
             a.diag([(F(n,2)-k)**2 for k in range(size)]))
    return squares,G


def spin_controls():
    inequalities=casimirs=0;odd_floors=0
    for n in range(1,13):
        squares,G=spin_squares(n);j=F(n,2);I=a.identity(n+1)
        total=a.add(a.add(squares[0],squares[1]),squares[2])
        assert total==a.scale(I,j*(j+1));casimirs+=1
        for sq in squares:
            assert a.psd(a.multiply(G,a.add(a.scale(I,j*j),sq,-1)))
            if n%2:
                assert a.psd(a.multiply(G,a.add(sq,a.scale(I,F(-1,4)))))
                odd_floors+=1
        for vals in DIAGONALS:
            L=a.scale(I,0)
            for v,sq in zip(vals,squares):L=a.add(L,a.scale(sq,v))
            full,even,odd=rates(vals)
            floor=odd if n%2 else even
            assert a.psd(a.multiply(G,a.add(L,a.scale(I,-floor))))
            if not n%2:
                bound=vals[0]*j*j+vals[1]*j
                assert a.psd(a.multiply(G,a.add(L,a.scale(I,-bound))))
                assert bound>=even
            assert floor>=full;inequalities+=1
    return {'finite_spin_blocks':casimirs,'degree_range':[1,12],
            'gap_matrix_inequalities':inequalities,
            'half_integer_axis_floor_checks':odd_floors,
            'all_degrees_supported_by':'YM52-T4 written block proof, not finite cutoff'}


@lru_cache(maxsize=None)
def gram(d):
    mons=[{m:F(1)} for m in y.basis(d)]
    G=[[y.phi(p.mul(f,g)) for g in mons] for f in mons]
    means=[y.phi(f) for f in mons]
    centered=[[G[i][j]-means[i]*means[j] for j in range(len(mons))]
              for i in range(len(mons))]
    return G,centered


def polynomial_gap_controls():
    controls=0;rows=[]
    for vals in DIAGONALS:
        C=a.diag(vals);full,even,odd=rates(vals)
        assert q.lap(p.COORD[0],C)==p.scale(p.COORD[0],odd)
        assert q.lap(EVEN,C)==p.scale(EVEN,even)
        assert variance(EVEN)==F(1,12) and y.phi(EVEN)==0
        for d in range(1,5):
            G,V=gram(d);L=q.operator(d,C);bound=odd if d%2 else even
            matrix=a.add(a.multiply(G,L),a.scale(V,-bound))
            assert a.psd(matrix);controls+=1
        rows.append({'eigenvalues':list(map(str,vals)),
                     'full':str(full),'even':str(even),'odd':str(odd)})
    C=q.rotate_tensor(a.diag(DIAGONALS[2]),ROT2)
    for d in range(1,5):
        G,V=gram(d);full,even,odd=rates(DIAGONALS[2])
        bound=odd if d%2 else even
        assert a.psd(a.add(a.multiply(G,q.operator(d,C)),a.scale(V,-bound)))
        controls+=1
    for eta in [F(0),F(1,1000),F(1,4),F(3,8),F(1,2),F(1)]:
        assert rates((eta,eta,3-2*eta))[0]==min(F(3,4),2*eta)
    assert rates((0,F(3,2),F(3,2)))[0]==F(3,4)
    return {'all_source_polynomial_matrix_checks':controls,
            'degree_range':[1,4],'sharp_mode_table':rows,
            'rank_two_positive_rate':'3/4','rank_one_centered_rate':'0'}


def observer_controls():
    f=p.add(p.COORD[0],p.scale(EVEN,2));n=0
    assert parity(f)==p.scale(EVEN,2)
    assert parity(f,False)==p.COORD[0]
    assert y.phi(p.mul(parity(f),parity(f,False)))==0
    for vals in DIAGONALS:
        C=a.diag(vals)
        assert parity(q.lap(f,C))==q.lap(parity(f),C)
        assert parity(q.lap(f,C),False)==q.lap(parity(f,False),C)
        assert variance(f)==variance(parity(f))+variance(parity(f,False))
        assert energy(f,C)==energy(parity(f),C)+energy(parity(f,False),C)
        n+=1
    plateau=unique=stable=0
    for denominator in (4,8,12):
        s=F(3)
        for ia in range(3*denominator+1):
            for ib in range(ia,3*denominator-ia+1):
                ic=3*denominator-ia-ib
                if ic<ib:continue
                vals=tuple(F(k,denominator) for k in (ia,ib,ic))
                full,even,_=rates(vals)
                assert full<=s/4 and (full==s/4)==(vals[2]<=3*s/4)
                assert even<=2*s/3 and (even==2*s/3)==(vals==(1,1,1))
                epsilon=2*s/3-even
                assert max(abs(v-s/3) for v in vals)<=2*epsilon
                plateau+=1;unique+=1;stable+=1
    for epsilon in (F(0),F(1,10),F(1,4),F(1,2)):
        vals=(1-2*epsilon,1+epsilon,1+epsilon)
        assert 2-rates(vals)[1]==epsilon
        assert max(abs(v-1) for v in vals)==2*epsilon
    assert rates((0,F(3,2),F(3,2)))[0]==rates((1,1,1))[0]
    assert rates((0,F(3,2),F(3,2)))[1]<rates((1,1,1))[1]
    return {'parity_energy_checks':n,'full_objective_plateau_checks':plateau,
            'even_unique_maximizer_checks':unique,'stability_checks':stable,
            'stability_constant_two_is_sharp':True,
            'full_gap_optimization_does_not_select_isotropy':True}


def quadratic_uniform_bound(f):
    """Rational sphere bound: |constant|+linear l1+max quadratic row l1."""
    linear=[F(0)]*4;A=a.scale(a.identity(4),0);constant=F(0)
    for m,value in f.items():
        degree=sum(m)
        if degree==0:constant+=value
        elif degree==1:linear[m.index(1)]+=value
        elif degree==2:
            slots=[i for i in range(4) for _ in range(m[i])]
            i,j=slots
            if i==j:A[i][i]+=value
            else:A[i][j]+=value/2;A[j][i]+=value/2
        else:raise ValueError('enclosure helper requires degree at most two')
    return abs(constant)+sum(map(abs,linear))+max(sum(map(abs,row)) for row in A)


def entropy_series(f,C,epsilon,order=12):
    """Prove a positive density range before making the scalar enclosures."""
    epsilon=q.exact(epsilon)
    r=abs(epsilon)*quadratic_uniform_bound(f)
    if not 0<r<1 or isinstance(order,bool) or not isinstance(order,int) or order<2:
        raise ValueError('certified |epsilon f|<1 and integer order>=2 required')
    if y.phi(f)!=0:raise ValueError('centered density perturbation required')
    powers=[y.poly_power(f,k) for k in range(order+1)]
    H=sum((F((-1)**k,k*(k-1))*epsilon**k*y.phi(powers[k])
           for k in range(2,order+1)),F(0))
    htail=r**(order+1)/F(order*(order+1))/(1-r)
    G=gamma(f,f,C);E=y.phi(G)
    I=epsilon**2*sum(((-epsilon)**k*y.phi(p.mul(powers[k],G))
                     for k in range(order+1)),F(0))
    itail=epsilon**2*r**(order+1)/(1-r)*E
    return (H-htail,H+htail),(I-itail,I+itail)


def entropy_controls():
    identities=enclosures=0;rows=[]
    functions=[p.COORD[0],p.scale(EVEN,2)]
    for vals in (DIAGONALS[2],DIAGONALS[4],DIAGONALS[7]):
        C=a.diag(vals)
        for f in functions:
            for k in range(2,9):
                lhs=y.phi(p.mul(y.poly_power(f,k-1),q.lap(f,C)))
                rhs=(k-1)*y.phi(p.mul(y.poly_power(f,k-2),gamma(f,f,C)))
                assert lhs==rhs;identities+=1
            epsilon=F(1,4);m=1-epsilon;M=1+epsilon
            H,I=entropy_series(f,C,epsilon)
            V=epsilon**2*variance(f);E=epsilon**2*energy(f,C)
            assert V/(2*M)<=H[0]<=H[1]<=V/(2*m)
            assert E/M<=I[0]<=I[1]<=E/m
            rate=rates(vals)[1 if parity(f)==f else 0]
            assert I[0]>=2*m*rate/M*H[1]
            rows.append({'tensor':list(map(str,vals)),
                'observer':'even' if parity(f)==f else 'full',
                'entropy_interval':list(map(str,H)),
                'production_interval':list(map(str,I)),
                'bounded_density_decay_lower':str(2*m*rate/M)})
            enclosures+=1
    return {'entropy_differentiation_coefficients':identities,
            'rational_entropy_production_enclosures':enclosures,
            'density_bounds':['3/4','5/4'],'cases':rows}


def scale_controls():
    f=p.add(p.COORD[0],EVEN);checks=0
    for vals in DIAGONALS[2:]:
        C=a.diag(vals)
        for k in (F(1,10),F(2),F(7,3)):
            scaled=a.scale(C,k)
            assert energy(f,scaled)==k*energy(f,C)
            assert rates(tuple(k*v for v in vals))==tuple(k*v for v in rates(vals))
            assert cofactor(scaled)==a.scale(cofactor(C),k*k)
            checks+=1
    return {'linear_clock_and_quadratic_bracket_scalings':checks,
            'physical_energy_and_clock_units_selected':False}


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
    return {'certificate_type':'YM52_NATIVE_ENERGY_NOISE_AND_OBSERVER_RELAXATION',
        'verdict':'PASS','runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
        'energy':energy_controls(),'discarded_records':record_controls(),
        'bracket_tensor':bracket_controls(),'spin_blocks':spin_controls(),
        'polynomial_gaps':polynomial_gap_controls(),'observer_selection':observer_controls(),
        'entropy':entropy_controls(),'clock':scale_controls(),
        'claim_status':'COMPACT_FREE_PROTOCOL_ENERGY_AND_OBSERVER_RATES_PROVED__'
                       'PHYSICAL_SELECTION_AND_INTERACTING_TRANSFER_OPEN',
        'evidence_scope':{
            'general_proof':'written, not mechanically formalized or expert-certified',
            'gap':'optimal above constants in the specified compact coefficient completion',
            'lineage':'positive-definite full/even spectral formula credited to Lauret',
            'selection':'even-observer objective and fixed trace budget are declared',
            'entropy':'defined relative entropy for bounded positive density; no physical units',
            'external_sources':'metadata pinned; CI does not recertify external theorems',
            'open':'physical energy/time/action; anisotropic interacting estimates; actual row '
                   'closure; four-dimensional QFT, asymptotic freedom and Clay'}}


def check(cert,result=RESULT,pin=PIN):
    sha=y.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM52 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=y.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM52 PASS',sha)
    print(json.dumps({'spin_gap_checks':cert['spin_blocks']['gap_matrix_inequalities'],
        'polynomial_gap_checks':cert['polynomial_gaps']['all_source_polynomial_matrix_checks'],
        'entropy_enclosures':cert['entropy']['rational_entropy_production_enclosures'],
        'physical_selection':'OPEN'},sort_keys=True))


if __name__=='__main__':main()
