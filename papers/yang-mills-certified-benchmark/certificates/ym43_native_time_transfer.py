"""YM-43: exact controls for the written native time-transfer bound.

Python 3.12 only in CI. Finite arithmetic checks the hypotheses and examples;
YM43_NATIVE_TIME_TRANSFER_BOUND.md contains the general proofs. The full
declared heat-kernel functional is not derived from the EMK primitives here.
Default/--check is read-only. --write explicitly generates the new record.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from ym1_certified_gap import Iv, log_iv, canonical_sha, LOG_TERMS
import ym19_dobrushin_dock as y19
import ym37_space_transfer as y37

RESULT = HERE / 'YM43_RESULT.json'
PIN = HERE / 'EXPECTED_YM43.sha256'
SOURCES = HERE / 'YM43_SOURCE_PINS.json'
SCALE = 10**12
R_SPACE, V_TIME = F(33,32), F(1,128)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def up(x, scale=SCALE):
    return F(-((-x.numerator * scale)//x.denominator), scale)


def down(x, scale=SCALE):
    return F((x.numerator * scale)//x.denominator, scale)


def decimals(x, upward=True, places=12):
    z = up(x, 10**places) if upward else down(x, 10**places)
    k = z.numerator * (10**places)//z.denominator
    sign = '-' if k < 0 else ''
    k = abs(k)
    return f'{sign}{k//10**places}.{k%10**places:0{places}d}'


def gate(cs, ct, q):
    if cs < 0 or ct < 0 or not 0 < q < 1:
        return None
    eta = 2*cs + ct*(q+1/q)
    return eta if eta < 1 else None


def coefficient_cell(a, kappa):
    S = y19.S_a(a)
    if S.hi >= 1:
        return {'status':'REFUSED_KERNEL_FLOOR'}
    ct = up(y19.coeff(y19.delta_t(a)).hi)
    cs = up(y19.coeff(y19.delta_s(kappa)).hi)
    alpha = 2*(cs+ct)
    choice = next(((F(1,2**n), n) for n in range(10,0,-1)
                   if gate(cs,ct,F(1,2**n)) is not None), None)
    common = {'a':str(a), 'kappa':str(kappa), 'c_space_upper':str(cs),
              'c_time_upper':str(ct), 'alpha_upper':str(alpha),
              'kernel_floor_lower':str(down(1-S.hi)),
              'heat_tail_included':True}
    if choice is None:
        assert alpha >= 1  # these refusal cells fail every q with this budget
        return {**common, 'status':'REFUSED_OVERLAP_BUDGET'}
    q,n = choice
    eta = gate(cs,ct,q)
    free = y19.exp_neg(F(3,4)*a)
    assert free.hi < q
    log2 = log_iv(Iv(F(2)), LOG_TERMS)
    rate = log2*Iv(F(n))
    # Independent native odd-series enclosure on q itself.
    direct = -log_iv(Iv(q), LOG_TERMS)
    assert rate.lo <= direct.hi and direct.lo <= rate.hi
    return {**common, 'status':'CERTIFIED_WITHIN_DECLARED_FUNCTIONAL_CARRIER',
            'q':str(q), 'eta_exact':str(eta),
            'eta_upper_decimal':decimals(eta),
            'gap_step_lower':decimals(rate.lo,False),
            'gap_per_a_lower':decimals(rate.lo/a,False),
            'free_B_ratio_upper':decimals(free.hi), 'free_B_below_q':True}


def barrier_controls(cs,ct,q):
    assert gate(cs,ct,q) is not None
    count = 0
    for m in (1,2,5,11):
        for L in (2,3,8,17):
            def v(r): return q**r + q**(L-r)
            assert v(0)>=1 and v(L)>=1
            for i in range(m):
                for r in range(1,L):
                    interior = cs*((i>0)+(i+1<m))*v(r)
                    boundary = F(0)
                    for rr in (r-1,r+1):
                        if rr in (0,L): boundary += ct
                        else: interior += ct*v(rr)
                    assert interior+boundary <= v(r)
                    assert v(r-1)+v(r+1) == (q+1/q)*v(r)
                    count += 2
    return count


def coupling(p,q):
    common = [min(a,b) for a,b in zip(p,q)]
    r = 1-sum(common)
    n = len(p)
    J = [[common[i] if i==j else F(0) for j in range(n)] for i in range(n)]
    if r:
        for i in range(n):
            for j in range(n):
                J[i][j] += (p[i]-common[i])*(q[j]-common[j])/r
    assert [sum(row) for row in J] == p
    assert [sum(J[i][j] for i in range(n)) for j in range(n)] == q
    assert all(x>=0 for row in J for x in row)
    assert sum(J[i][j] for i in range(n) for j in range(n) if i!=j) == r
    return J


def ws(x,y): return R_SPACE**2 if x==y else F(1)
def wt(x,y): return 1+V_TIME*x*y


def normalized(values):
    z = sum(values)
    assert z > 0 and all(v>0 for v in values)
    return [v/z for v in values]


def conditional(neighbours, bias=(F(2),F(3))):
    # First two neighbours are spatial; last two temporal.
    values=[]
    for x,b in zip((-1,1),bias):
        values.append(b*ws(x,neighbours[0])*ws(x,neighbours[1])
                      *wt(x,neighbours[2])*wt(x,neighbours[3]))
    return normalized(values)


def overlap_controls():
    cs,ct = 1-R_SPACE**-4, 1-((1-V_TIME)/(1+V_TIME))**2
    count=0
    for neighbours in itertools.product((-1,1),repeat=4):
        p=conditional(neighbours)
        for j in range(4):
            changed=list(neighbours);changed[j]*=-1
            q=conditional(changed)
            J=coupling(p,q)
            tv=sum(abs(a-b) for a,b in zip(p,q))/2
            assert tv <= (cs if j<2 else ct)
            assert J[0][1]+J[1][0] == tv
            count+=1
    return count


def update_controls():
    states=list(itertools.product((-1,1),repeat=2))
    def law(left,right):
        return normalized([wt(left,x)*wt(x,y)*wt(y,right) for x,y in states])
    mu,nu=law(-1,-1),law(1,1)
    joint=[[a*b for b in nu] for a in mu]
    ct=1-((1-V_TIME)/(1+V_TIME))**2
    counts=0
    for step in range(5):
        old=[sum(joint[i][j] for i,x in enumerate(states) for j,y in enumerate(states)
                 if x[k]!=y[k]) for k in (0,1)]
        next_joint=[[F(0)]*4 for _ in range(4)]
        for i,x in enumerate(states):
            for j,y in enumerate(states):
                for site in (0,1):
                    def cond(s,left,right):
                        return normalized([wt(left,z)*wt(z,s[1]) if site==0
                                           else wt(s[0],z)*wt(z,right) for z in (-1,1)])
                    J=coupling(cond(x,-1,-1),cond(y,1,1))
                    for a,sx in enumerate((-1,1)):
                        for b,sy in enumerate((-1,1)):
                            xx=list(x);yy=list(y);xx[site]=sx;yy[site]=sy
                            next_joint[states.index(tuple(xx))][states.index(tuple(yy))] += joint[i][j]*J[a][b]/2
        assert [sum(row) for row in next_joint]==mu
        assert [sum(next_joint[i][j] for i in range(4)) for j in range(4)]==nu
        for k in (0,1):
            now=sum(next_joint[i][j] for i,x in enumerate(states) for j,y in enumerate(states) if x[k]!=y[k])
            assert now <= old[k]/2 + ct*old[1-k]/2 + ct/2
            counts+=1
        joint=next_joint
    return {'joint_updates':5,'mismatch_recurrences':counts,'stationary_marginals_exact':True}


def transfer(m):
    states=list(itertools.product((-1,1),repeat=m))
    root=[R_SPACE**sum(x[i]==x[i+1] for i in range(m-1)) for x in states]
    T=[]
    for i,x in enumerate(states):
        row=[]
        for j,y in enumerate(states):
            v=root[i]*root[j]
            for a,b in zip(x,y):v*=wt(a,b)/2
            row.append(v)
        T.append(row)
    return states,T


def multiply(A,B):
    return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]


def marginal_controls():
    m,L=2,4
    states,T=transfer(m); n=len(states)
    powers=[[[F(i==j) for j in range(n)] for i in range(n)]]
    for _ in range(L):powers.append(multiply(powers[-1],T))
    all_marginals={r:[] for r in range(1,L)}
    identities=0
    for a,x in enumerate(states):
        for b,z in enumerate(states):
            marg=[[F(0)]*n for _ in range(L-1)]
            Z=F(0)
            for middle in itertools.product(states,repeat=L-1):
                rows=(x,)+middle+(z,);value=F(1)
                for row in middle:value*=ws(row[0],row[1])
                for u,v in zip(rows,rows[1:]):
                    for s in range(m):value*=wt(u[s],v[s])
                Z+=value
                for r,row in enumerate(middle):marg[r][states.index(row)]+=value
            for r in range(1,L):
                direct=[v/Z for v in marg[r-1]]
                matrix=[powers[r][a][j]*powers[L-r][j][b]/powers[L][a][b] for j in range(n)]
                assert direct==matrix
                identities+=n
                all_marginals[r].append(direct)
    cs,ct=1-R_SPACE**-4,1-((1-V_TIME)/(1+V_TIME))**2
    q=F(1,16);assert gate(cs,ct,q) is not None
    comparisons=0
    for r,laws in all_marginals.items():
        for p,p2 in itertools.combinations(laws,2):
            tv=sum(abs(a-b) for a,b in zip(p,p2))/2
            assert tv<=m*(q**r+q**(L-r))
            comparisons+=1
    return {'width':m,'time_length':L,'path_to_transfer_equalities':identities,
            'all_row_observable_TV_comparisons':comparisons,'q':str(q)}


def positive_kernel(T):
    n=len(T)
    return n>0 and all(len(row)==n for row in T) and all(
        T[i][j]>0 and T[i][j]==T[j][i] for i in range(n) for j in range(n))


def vacuum_controls(T,steps=7):
    if not positive_kernel(T):raise ValueError('positive symmetric kernel required')
    n=len(T);ell=min(x for row in T for x in row);U=max(x for row in T for x in row)
    zeta=1-(ell/U)**2
    f=[F(1)]*n; prev=None; count=0
    for k in range(steps):
        fn=y37.matvec(T,f);r=[a/b for a,b in zip(fn,f)]
        lo,hi=min(r),max(r)
        h=[v/(sum(f)/n) for v in f]; hn=[v/(sum(fn)/n) for v in fn]
        if k==0: initial_lo,initial_width=lo,hi-lo
        else:
            assert prev[0]<=lo<=hi<=prev[1]
            assert hi-lo<=zeta*(prev[1]-prev[0])
            assert max(abs(a-b) for a,b in zip(hn,h)) <= U/ell*initial_width/initial_lo*zeta**k
            count+=3
        assert all(ell/U<=v<=U/ell for v in hn)
        count+=1;prev=(lo,hi);f=fn
    return count


def transfer_controls():
    q=F(1,16);rows=[]
    for m in (1,2,3,4):
        _,T=transfer(m);n=len(T)
        assert y37.inertia(T)==(n,0,0)
        lower=min(sum(row) for row in T)
        ceiling=q*lower
        signs=y37.inertia([[T[i][j]-(ceiling if i==j else 0) for j in range(n)] for i in range(n)])
        assert signs==(1,n-1,0)
        # An all-space q bound would incorrectly remove the vacuum.
        assert signs[0]!=0
        checks=vacuum_controls(T)
        # Cut-square log convexity, tested on an unselected mixed source.
        v=[F((i%5)-2) for i in range(n)];mom=[]
        for _ in range(7):
            mom.append(sum(x*x for x in v));v=y37.matvec(T,v)
        assert all(mom[k]**2<=mom[k-1]*mom[k+1] for k in range(1,6))
        rows.append({'m':m,'dimension':n,'q':str(q),'vacuum_floor':str(lower),
                     'threshold_inertia':list(signs),'vacuum_iteration_checks':checks,
                     'log_convexity_checks':5})
    return rows


def negative_controls():
    c=F(1,8);q=F(1,16)
    assert 2*c+c*q<1 and gate(c,c,q) is None  # missing return-time term
    assert gate(F(1,3),F(1,3),F(1,2)) is None
    assert all(gate(F(0),F(0),q) is None for q in (F(0),F(1),F(2)))
    p,p2=[F(1,3),F(2,3)],[F(3,4),F(1,4)]
    a=[min(x,y) for x,y in zip(p,p2)];r=1-sum(a)
    bad=[[a[i] if i==j else F(0) for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):bad[i][j]+=(p[i]-a[i])*(p2[j]-a[j]) # omitted division by r
    assert [sum(row) for row in bad]!=p
    # A blind probe vanishes while a different mean-zero source is slow.
    eps=F(1,10)
    P=[[(1-eps)*F(int(i//2==j//2),2)+eps/4 for j in range(4)] for i in range(4)]
    f,g=[F(1),F(-1),F(0),F(0)],[F(1),F(1),F(-1),F(-1)]
    assert y37.matvec(P,f)==[0]*4
    assert y37.matvec(P,g)==[(1-eps)*v for v in g]
    assert sum(v*v for v in y37.matvec(P,g))>F(1,4)*sum(v*v for v in g)
    for bad_kernel in ([[F(1),F(-1)],[F(-1),F(1)]],[[F(2),F(1)],[F(2),F(2)]]):
        try:vacuum_controls(bad_kernel)
        except ValueError:pass
        else:raise AssertionError('invalid kernel was accepted')
    return {'omitted_return_time_term':True,'noncontracting_budget':True,
            'invalid_q':True,'wrong_common_part_normalization':True,
            'hidden_slow_source':True,'nonpositive_kernel':True,'nonsymmetric_kernel':True,
            'vacuum_projection_omitted':True}


def source_checks():
    pins=json.loads(SOURCES.read_text());checked={}
    for path,expected in pins['upstream_sha256'].items():
        actual=digest(ROOT/path)
        if actual!=expected:raise ValueError('upstream source changed: '+path)
        checked[path]=actual
    for path in pins['local_inputs']:checked[path]=digest(ROOT/path)
    checked[str(SOURCES.relative_to(ROOT))]=digest(SOURCES)
    old42=json.loads((HERE/'YM42_RESULT.json').read_text());old42.pop('input_sha256')
    if canonical_sha(old42)!=pins['ym42_runtime_only']['mathematical_payload_sha256']:
        raise ValueError('YM42 mathematical output changed during runtime migration')
    return checked


def run():
    if not __debug__:raise RuntimeError('Optimized Python is refused')
    inputs=source_checks();grid={};barriers=0
    for a,k in y19.GRID:
        row=coefficient_cell(a,k);assert row['status'].startswith('CERTIFIED')
        barriers+=barrier_controls(F(row['c_space_upper']),F(row['c_time_upper']),F(row['q']))
        grid[f'a={a},kappa={k}']=row
    refused={f'a={a},kappa={k}':coefficient_cell(a,k) for a,k in y19.FAIL_CELLS}
    assert all(r['status'].startswith('REFUSED') for r in refused.values())
    return {'certificate_type':'YM43_NATIVE_TIME_TRANSFER_BOUND','verdict':'PASS',
            'runtime_policy':'Python 3.12 only','input_sha256':inputs,
            'grid':grid,'refused':refused,'strip_barrier_equalities_and_bounds':barriers,
            'common_part_coupling_cases':overlap_controls(),'conditional_updates':update_controls(),
            'independent_row_marginals':marginal_controls(),'finite_transfer_controls':transfer_controls(),
            'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_GENERAL_PROOF_WITHIN_DECLARED_POSITIVE_FUNCTIONAL_CARRIER__'
                           'FULL_TIME_OPERATOR_BOUND_ALL_FINITE_WIDTHS_AT_CERTIFIED_COARSE_CELLS__'
                           'NOT_PRIMITIVE_MEASURE_DERIVATION__CONTINUUM_NG_ALL_A_CLAY_OPEN',
            'evidence_scope':{'general_theorem':'written proof, not formalized',
                'machine_certificate':'exact parameter hypotheses and finite identity/refusal controls',
                'YM19_historical_record':'preserved; external comparison premise replaced on this scoped route',
                'YM42_original_Wilson_grid':'no transfer of this heat-kernel result claimed',
                'native_measure_dictionary':'OPEN','cutoff_uniformity':'OPEN','physical_measurement':False}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM43 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group();group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM43 PASS',sha)
    print(json.dumps({k:{name:r[name] for name in ('q','eta_upper_decimal','gap_step_lower')}
                      for k,r in cert['grid'].items()},sort_keys=True))


if __name__=='__main__':main()
