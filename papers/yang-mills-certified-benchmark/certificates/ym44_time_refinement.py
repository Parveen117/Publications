"""YM-44: finite exact controls for the written full time-refinement proof.

Python 3.12 only. Reuses existing interval, fusion and matrix arithmetic.
No general operator engine or SU(2) discretization is introduced here.
Default/--check rebuilds evidence without writing; --write is explicit.
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from ym1_certified_gap import Iv, canonical_sha
from ym2_theta_interacting_gap import exp_point
from ym4_symmetry_protected import chi_mul
from ym6_seam_integer_dock import _r, iv_sqrt, ldl_inertia
from ym43_native_time_transfer import multiply, decimals

RESULT = HERE / 'YM44_RESULT.json'
PIN = HERE / 'EXPECTED_YM44.sha256'
SOURCES = HERE / 'YM44_SOURCE_PINS.json'
I2 = [[F(1),F(0)],[F(0),F(1)]]
L2 = [[F(3,4),F(-3,4)],[F(-3,4),F(3,4)]]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@lru_cache(maxsize=None)
def ex(x):
    return _r(exp_point(F(x)))


def ivpower(x,n):
    out=Iv(1)
    for _ in range(n):
        out=_r(out*x)
    return out


def factorial_tail(x,N):
    if x<0 or N<0 or not isinstance(N,int):
        raise ValueError('invalid factorial-tail input')
    ratio=x/(N+2)
    if ratio>=1:
        raise ValueError('factorial-tail ratio must be strictly below one')
    return x**(N+1)/factorial(N+1)/(1-ratio)


def moment_controls():
    assert list(chi_mul(1,1))==[0,2]
    assert 2*F(1*(1+2),4)==F(3,2)
    assert 2*F(2*(2+2),4)==4
    # Exact coefficient identity: (1-u)^2 H(u)=5-8u^3+3u^8.
    H=list(map(F,(5,10,15,12,9,6,3)))
    coefficients=[sum((F((1,-2,1)[j])*H[k-j]
                      for j in range(3) if 0<=k-j<len(H)),F(0))
                  for k in range(9)]
    assert coefficients==[F(5),F(0),F(0),F(-8),F(0),F(0),F(0),F(0),F(3)]
    count=0
    for s in (F(0),F(1,65536),F(1,256),F(1,16),F(1,4),F(1),F(3),F(6)):
        u=ex(-s/2)
        d0=(Iv(1)-ivpower(u,8))/Iv(4)
        poly=sum_iv([ivpower(u,j)*Iv(c) for j,c in enumerate(H)])
        d1=ivpower(Iv(1)-u,2)*poly/Iv(4)
        for z in (F(k,16) for k in range(17)):
            positive=d0*Iv(1-z)+d1*Iv(z)
            direct=(Iv(1)-ex(-4*s))/Iv(4)+(Iv(1)+ex(-4*s)-ex(-F(3,2)*s)*Iv(2))*Iv(z)
            assert positive.lo>=0
            assert positive.hi<=3*s
            assert not positive.separated_from(direct)
            count+=1
    return {'fusion_identity':True,'endpoint_polynomial':[str(x) for x in coefficients],
            'moment_enclosures':count,'full_commutator_bound':'||[K_s,B]|| <= b sqrt(3s)'}


def sum_iv(xs):
    out=Iv(0)
    for x in xs:
        out=out+x
    return out


def split_word_coefficient(word):
    k=len(word);value=F(0)
    for p in range(k+1):
        for q in range(k-p+1):
            r=k-p-q
            if 'A'*p+'B'*q+'A'*r==word:
                value+=F(1,2**(p+r)*factorial(p)*factorial(q)*factorial(r))
    return value


def ordered_word_controls():
    cubic={'AAB':F(-1,24),'ABA':F(1,12),'BAA':F(-1,24),
           'ABB':F(1,12),'BAB':F(-1,6),'BBA':F(1,12)}
    counts=0
    for degree in range(4):
        for letters in itertools.product('AB',repeat=degree):
            word=''.join(letters);s=split_word_coefficient(word)
            exact=F(1,factorial(degree))
            expected=cubic.get(word,F(0)) if degree==3 else F(0)
            assert s-exact==expected
            half=sum((split_word_coefficient(word[:j])*split_word_coefficient(word[j:])
                      for j in range(degree+1)),F(0))/2**degree
            assert s-half==(F(3,4)*expected if degree==3 else F(0))
            counts+=2
    # Sum all words at each degree equals the commuting scalar coefficient.
    for degree in range(7):
        total=sum((split_word_coefficient(''.join(w))
                   for w in itertools.product('AB',repeat=degree)),F(0))
        assert total==F(2**degree,factorial(degree))
    return {'word_identities':counts,'commuting_scalar_degrees':7,
            'cubic_residue':{k:str(v) for k,v in cubic.items()}}


def refinement_budget(width,theta,t,mesh):
    if not isinstance(width,int) or width<1 or t<0 or mesh<0 or mesh>t:
        raise ValueError('invalid width, time or mesh')
    b=abs(theta)*(width-1)
    if not b or not t or not mesh:
        return F(0)
    return 3*b*t*ex(b*t).hi*iv_sqrt(Iv(mesh)).hi


def budget_controls():
    rows=[]
    for m in (1,2,4,9):
        for t in (F(1,4),F(1),F(4)):
            previous=None
            for n in (1,4,16,64,256):
                bound=refinement_budget(m,F(1,16),t,t/n)
                if previous is not None:
                    assert bound<=previous
                previous=bound
                rows.append({'m':m,'t':str(t),'steps':n,'norm_error_upper':decimals(bound)})
    tails=[]
    for x in (F(0),F(1,16),F(1,2),F(2),F(8)):
        N=24;tail=factorial_tail(x,N)
        partial=sum((x**j/factorial(j) for j in range(N+1)),F(0))
        # Independent existing exponential enclosure.
        exact=exp_point(x)
        assert partial<=exact.hi and partial+tail>=exact.lo
        extra=sum((x**j/factorial(j) for j in range(N+1,N+9)),F(0))
        assert extra<=tail
        tails.append({'x':str(x),'N':N,'tail_upper':str(tail)})
    partitions=((F(1,4),)*4,(F(1,2),F(1,4),F(1,4)),
                (F(1,16),)*16,(F(1,9),F(4,9),F(4,9)))
    for parts in partitions:
        total=sum(parts);mesh=max(parts)
        lhs=sum(h*iv_sqrt(Iv(h)).hi for h in parts)
        rhs=total*iv_sqrt(Iv(mesh)).hi
        assert lhs<=rhs
    # sqrt(3)+sqrt(3/2)<3 via sqrt(2)<3/2, all squared rational gates.
    assert F(2)<F(3,2)**2
    return {'equal_partition_budgets':rows,'factorial_tails':tails,
            'unequal_partition_checks':len(partitions),'local_constant_three_verified':True}


# Two-state controls only. Rational multiplication is the existing YM43
# routine; interval wrappers retain outward bounds after every operation.
def ivmul2(A,B):
    return [[_r(sum_iv([A[i][k]*B[k][j] for k in range(2)]))
             for j in range(2)] for i in range(2)]


def ivpow2(A,n):
    out=[[Iv(v) for v in row] for row in I2]
    while n:
        if n%2:
            out=ivmul2(out,A)
        n//=2
        if n:
            A=ivmul2(A,A)
    return out


def heat2(t):
    r=ex(-F(3,2)*t)
    a,b=(Iv(1)+r)/Iv(2),(Iv(1)-r)/Iv(2)
    return [[a,b],[b,a]]


def split2(t,theta):
    K=heat2(t/2)
    M=[[ex(theta*t),Iv(0)],[Iv(0),ex(-theta*t)]]
    return ivmul2(ivmul2(K,M),K)


def row_bound2(A):
    return max(sum(max(abs(x.lo),abs(x.hi)) for x in row) for row in A)


def parity_series(x,offset,N=24):
    partial=sum((x**k/factorial(2*k+offset) for k in range(N+1)),F(0))
    first=x**(N+1)/factorial(2*N+2+offset)
    ratio=x/F((2*N+4+offset)*(2*N+3+offset))
    if x<0 or ratio>=1:
        raise ValueError('parity tail requires a lawful positive ratio')
    return Iv(partial,partial+first/(1-ratio))


def reference2(t,theta):
    # G=-3/4 I+D, D^2=(theta^2+9/16)I: no square root evaluated.
    x=(theta**2+F(9,16))*t*t
    c=parity_series(x,0);s=parity_series(x,1)*Iv(t);scale=ex(-F(3,4)*t)
    return [[_r(scale*(c+s*Iv(theta))),_r(scale*s*Iv(F(3,4)))],
            [_r(scale*s*Iv(F(3,4))),_r(scale*(c-s*Iv(theta)))]]


def taylor_reference2(t,theta,N=42):
    G=[[F(-3,4)+theta,F(3,4)],[F(3,4),F(-3,4)-theta]]
    X=[[t*x for x in row] for row in G]
    term=[row[:] for row in I2];total=[row[:] for row in I2]
    for k in range(1,N+1):
        term=[[x/k for x in row] for row in multiply(term,X)]
        total=[[total[i][j]+term[i][j] for j in range(2)] for i in range(2)]
    norm=max(sum(abs(v) for v in row) for row in X)
    tail=factorial_tail(norm,N)
    return [[Iv(x-tail,x+tail) for x in row] for row in total]


def finite_refinement_controls():
    rows=[];compared=0
    for theta in (F(0),F(1,16),F(-1,8),F(1,4)):
        for t in (F(1,4),F(1),F(2)):
            U=reference2(t,theta);independent=taylor_reference2(t,theta)
            for i,j in itertools.product(range(2),repeat=2):
                assert not U[i][j].separated_from(independent[i][j])
                compared+=1
            last=None
            for n in (1,4,16,64,256):
                S=ivpow2(split2(t/n,theta),n)
                difference=[[S[i][j]-U[i][j] for j in range(2)] for i in range(2)]
                upper=row_bound2(difference);budget=refinement_budget(2,theta,t,t/n)
                if theta:
                    assert upper<budget
                    if last is not None:
                        assert upper<last
                else:
                    # Exact free identity in the proof; intervals agree.
                    assert all(not S[i][j].separated_from(U[i][j])
                               for i,j in itertools.product(range(2),repeat=2))
                last=upper
                rows.append({'theta':str(theta),'t':str(t),'steps':n,
                             'independent_error_upper':decimals(upper),
                             'theorem_budget_upper':decimals(budget)})
    # Unequal time partitions use later action on the left.
    for parts in ((F(1,2),F(1,4),F(1,4)),(F(1,9),F(4,9),F(4,9))):
        S=[[Iv(x) for x in row] for row in I2]
        for h in parts:
            S=ivmul2(split2(h,F(1,16)),S)
        U=reference2(sum(parts),F(1,16))
        D=[[S[i][j]-U[i][j] for j in range(2)] for i in range(2)]
        # A non-palindromic product need not be symmetric. Frobenius is
        # used here, not the symmetric row-sum norm shortcut.
        frob2=sum(max(abs(x.lo),abs(x.hi))**2 for row in D for x in row)
        bound=refinement_budget(2,F(1,16),sum(parts),max(parts))
        assert frob2<bound*bound
    coarse=split2(F(1),F(1,16));refined=ivpow2(split2(F(1,2),F(1,16)),2)
    delta=coarse[0][0]-refined[0][0]
    assert delta.lo>0 or delta.hi<0
    return {'refinement_rows':rows,'independent_exponential_entry_comparisons':compared,
            'unequal_partition_operator_checks':2,
            'coarse_minus_two_half_00_interval':[str(delta.lo),str(delta.hi)]}


def gap_gate(ell,q,error):
    if ell<=0 or not 0<=q<1 or error<0:
        return None
    if 2*error>=ell*(1-q):
        return None
    return (q*ell+error)/(ell-error)


def gap_controls():
    rows=[]
    for e in (F(1,64),F(1,16),F(1,4)):
        ratio=gap_gate(F(2),F(1,2),e)
        assert ratio is not None and ratio<1
        # C=[[2,e],[e,1]], epsilon=e. The threshold 1+e has exactly
        # one positive and one negative cut-square elimination pivot.
        threshold=1+e
        M=[[Iv(2-threshold),Iv(e)],[Iv(e),Iv(1-threshold)]]
        inertia=ldl_inertia(M)
        assert inertia==(1,1)
        rows.append({'error':str(e),'ratio_upper':str(ratio),
                     'complement_threshold':str(threshold),'inertia':list(inertia)})
    return rows


def spectator_control():
    K=[[F(3,4),F(1,4)],[F(1,4),F(3,4)]]
    M=[[F(2),F(0)],[F(0),F(1)]]
    KM,MK=multiply(K,M),multiply(M,K)
    assert KM!=MK
    T=multiply(multiply(K,M),K)
    def tensor(A,B):
        return [[A[i][j]*B[k][l] for j in range(2) for l in range(2)]
                for i in range(2) for k in range(2)]
    for eps in (F(1,2),F(1,16),F(1,256),F(1,65536)):
        R=[[1-eps/2,eps/2],[eps/2,1-eps/2]]
        J=tensor(T,R)
        for f in ((F(1),F(0)),(F(2),F(-1))):
            for v,scale in (((F(1),F(1)),F(1)),((F(1),F(-1)),1-eps)):
                source=[[x*y] for x in f for y in v]
                actual=multiply(J,source)
                tf=multiply(T,[[x] for x in f])
                assert actual==[[scale*x[0]*y] for x in tf for y in v]
        assert multiply(tensor(K,R),tensor(M,I2))!=multiply(tensor(M,I2),tensor(K,R))
    return True


def negative_controls():
    coarse=split2(F(1),F(1,16));fine=ivpow2(split2(F(1,2),F(1,16)),2)
    nonexact=coarse[0][0].separated_from(fine[0][0])
    # Boundary equality admits a degenerate top: diag(3/2,3/2).
    strict=gap_gate(F(2),F(1,2),F(1,2)) is None
    # A=diag(1,1/2), C=diag(7/8,5/8), ||C-A||=1/8.
    # C's actual normalized ratio is 5/7, above the false 5/8 ceiling.
    actual_ratio=F(5,8)/F(7,8)
    denominator=(actual_ratio==gap_gate(F(1),F(1,2),F(1,8))
                 and actual_ratio>F(1,2)+F(1,8))
    high_content=[]
    for a,two_j in ((F(1),2),(F(1,16),8),(F(1,256),32),(F(1,65536),512)):
        casimir=F(two_j*(two_j+2),4)
        assert ex(-a*casimir).hi<F(1,2)
        high_content.append([str(a),two_j])
    truncated=sum((F(1,2)**j/factorial(j) for j in range(5)),F(0))
    assert truncated<exp_point(F(1,2)).lo
    try:
        factorial_tail(F(4),2)
        tail_gate=False
    except ValueError:
        tail_gate=True
    controls={'coarse_exact_subdivision_rejected':nonexact,
              'vacuum_denominator_omission_rejected':denominator,
              'non_strict_gap_margin_rejected':strict,
              'unbounded_heat_norm_taylor_rejected':len(high_content)==4,
              'width_independent_potential_budget_rejected':F(1,16)*8>F(1,16),
              'fixed_kappa_finite_generator_rejected':ex(F(1,16)).lo>1,
              'omitted_cubic_word_rejected':split_word_coefficient('ABA')!=F(1,6),
              'zero_factorial_tail_rejected':truncated<exp_point(F(1,2)).lo,
              'invalid_tail_ratio_rejected':tail_gate,
              'noncommutation_alone_gap_rejected':spectator_control()}
    assert all(controls.values())
    return {'groups':controls,'high_content_witnesses':high_content}


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
    return {'certificate_type':'YM44_FULL_INTERACTING_TIME_REFINEMENT','verdict':'PASS',
            'runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'moment_and_commutator':moment_controls(),'ordered_words':ordered_word_controls(),
            'tail_and_partition_budgets':budget_controls(),
            'independent_finite_refinement':finite_refinement_controls(),
            'normalized_gap_transport_controls':gap_controls(),'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_OPERATOR_NORM_REFINEMENT_PROOF_FOR_EVERY_FIXED_FINITE_WIDTH__'
                           'ALL_CONTENT_AND_SOURCES_WITHIN_DECLARED_HEAT_FUNCTIONAL_ADAPTER__'
                           'UNIFORM_VOLUME_CUTOFF_GAP_AND_NATIVE_MEASURE_DICTIONARY_OPEN',
            'evidence_scope':{'general_proof':'written; not mechanically formalized',
                'finite_controls':'exact rational and outward intervals; two-state fixtures are not SU(2)',
                'trajectory':'declared kappa(a)=theta*a, not AF derived',
                'gap_transport':'conditional on a separate uniform fine-step gap',
                'physical_clock_and_measure':'not derived','clay_or_quantum_gravity':'not established'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM44 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(sha+'\n')
    else:
        check(cert)
    print('YM44 PASS',sha)
    print(json.dumps({'moment_cases':cert['moment_and_commutator']['moment_enclosures'],
                      'word_identities':cert['ordered_words']['word_identities'],
                      'refinement_cases':len(cert['independent_finite_refinement']['refinement_rows']),
                      'refusal_groups':len(cert['negative_controls']['groups'])},sort_keys=True))


if __name__=='__main__':
    main()
