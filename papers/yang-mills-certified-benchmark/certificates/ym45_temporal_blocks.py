"""YM-45: exact controls for the written temporal-block uniform-gap proof.

Python 3.12 only. Finite binary fixtures test coupling and normalization;
they are not a discretization or replacement of the admitted SU(2) carrier.
Default/--check is read-only. No shared operator engine is introduced.
"""
import argparse
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
from ym1_certified_gap import Iv, LOG_TERMS, canonical_sha, log_iv
from ym19_dobrushin_dock import S_a
from ym43_native_time_transfer import coupling, decimals, multiply, vacuum_controls
from ym44_time_refinement import ex, gap_gate, refinement_budget

RESULT = HERE / 'YM45_RESULT.json'
PIN = HERE / 'EXPECTED_YM45.sha256'
SOURCES = HERE / 'YM45_SOURCE_PINS.json'
EPS = F(4,5)
CELLS = ((F(1,8192),11,F(2)), (F(1,4096),8,F(4,3)),
         (F(1,2048),6,F(21,20)), (F(1,1792),5,F(33,32)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def block_gate(theta,J,R):
    if not isinstance(J,int) or J<1 or R<=1:
        return None
    alpha=F(2)/(EPS*J)+168*abs(theta)*J
    if R*alpha>=1:
        return None
    gamma=log_iv(Iv(R),LOG_TERMS)/Iv(7*J)
    # The independent free B-sector must obey the same ceiling.
    if gamma.lo<=0 or gamma.hi>=F(3,4):
        return None
    return {'alpha':alpha,'margin':1-R*alpha,'gamma':gamma}


def parameter_controls():
    heat=S_a(F(6))
    assert heat.hi<F(1,20)
    assert F(19,21)**2>=EPS
    cells=[]
    for theta,J,R in CELLS:
        gate=block_gate(theta,J,R)
        assert gate is not None
        cases=0
        for a in (F(1),F(2,3),F(1,7),F(17,256),F(1,1024),F(1,65536)):
            r=-((-6*a.denominator)//a.numerator)
            ell=J*r
            assert 6<=a*r<7
            alpha=2*r/(EPS*ell)+24*a*abs(theta)*ell
            assert alpha<=gate['alpha'] and R*alpha<1
            assert a*ell<=7*J
            cases+=1
        cells.append({'theta_abs':str(theta),'J':J,'R':str(R),
                      'alpha_upper':str(gate['alpha']),'weighted_margin':str(gate['margin']),
                      'gamma_lower':decimals(gate['gamma'].lo,False),
                      'gamma_upper':decimals(gate['gamma'].hi),
                      'fine_step_budget_cases':cases})
    # Polynomial identity: this budget cannot pass at/above the endpoint.
    # J*(alpha-1) = (J-5)^2/10 + (168*theta-1/10)*J^2.
    for theta in (F(0),F(1,4096),F(1,1680),F(1,64),F(1,16)):
        for J in range(1,31):
            alpha=F(5,2*J)+168*theta*J
            assert J*(alpha-1)==F((J-5)**2,10)+(168*theta-F(1,10))*J*J
    # An interval statement, not an inference from these samples:
    # J=5, alpha=1/2+840|theta|<1 for |theta|<1/1680.
    for theta in (F(0),F(1,100000),F(1,2000),F(999,1680000)):
        alpha=F(1,2)+840*theta
        R=(1+1/alpha)/2
        assert R>1 and R*alpha==(1+alpha)/2<1
        assert block_gate(theta,5,R) is not None
    return {'full_heat_tail_at_six_upper':decimals(heat.hi),
            'heat_floor':'19/20','heat_ceiling':'21/20','skeleton_overlap':str(EPS),
            'open_window':'abs(theta) < 1/1680','cells':cells,
            'endpoint_budget_identity_checks':150}


def normalized(values):
    z=sum(values,F(0))
    if z<=0 or any(x<0 for x in values):
        raise ValueError('invalid positive normalization')
    return [x/z for x in values]


@lru_cache(maxsize=None)
def power_kernel(n,q=F(1,2)):
    if not isinstance(n,int) or n<0 or not 0<=q<1:
        raise ValueError('invalid finite heat fixture')
    return [[(1+q**n)/2,(1-q**n)/2],[(1-q**n)/2,(1+q**n)/2]]


def bridge_transition(x,z,step,remaining,q=F(1,2)):
    A,B=power_kernel(step,q),power_kernel(remaining,q)
    return normalized([A[x][y]*B[y][z] for y in range(2)])


def bridge_law(k,left,right,q=F(1,2),delta=F(0),env=None):
    paths=list(itertools.product(range(2),repeat=k));weights=[]
    K=power_kernel(1,q)
    for path in paths:
        full=(left,)+path+(right,);w=F(1)
        for x,y in zip(full,full[1:]):
            w*=K[x][y]
        if env is not None:
            for u,x in enumerate(path):
                for y in env[u]:
                    w*=1+delta*(2*x-1)*(2*y-1)
        weights.append(w)
    return paths,normalized(weights)


def hamming(x,y):
    return sum(a!=b for a,b in zip(x,y))


def coupling_cost(paths,p,q):
    joint=coupling(p,q)
    return sum((joint[i][j]*hamming(x,y)
                for i,x in enumerate(paths) for j,y in enumerate(paths)),F(0))


def bridge_controls():
    r=5;q=F(1,2)
    assert ((1-q**r)/(1+q**r))**2>=EPS
    K=power_kernel(1,q);P=[[F(1),F(0)],[F(0),F(1)]]
    for n in range(1,26):
        P=multiply(P,K)
        assert P==power_kernel(n,q)
    cases=0;max_cost=F(0);enumerated=0
    for k in (1,4,5,6,9,10,14,20,24):
        for right in (0,1):
            joint=[[F(0),F(1)],[F(0),F(0)]]
            t=0;cost=F(0);rounds=0
            while t+r<=k:
                mismatch=joint[0][1]+joint[1][0]
                cost+=r*mismatch
                nxt=[[F(0),F(0)],[F(0),F(0)]]
                for x in range(2):
                    for y in range(2):
                        px=bridge_transition(x,right,r,k+1-t-r,q)
                        py=bridge_transition(y,right,r,k+1-t-r,q)
                        C=coupling(px,py)
                        if x!=y:
                            assert C[0][0]+C[1][1]>=EPS
                        for u in range(2):
                            for v in range(2):
                                nxt[u][v]+=joint[x][y]*C[u][v]
                joint=nxt;t+=r;rounds+=1
                assert joint[0][1]+joint[1][0]<=(1-EPS)**rounds
                for left,axis in ((0,0),(1,1)):
                    actual=([sum(row) for row in joint] if axis==0 else
                            [sum(joint[i][j] for i in range(2)) for j in range(2)])
                    expect=bridge_transition(left,right,t,k+1-t,q)
                    assert actual==expect
                    if k<=9:
                        paths,law=bridge_law(k,left,right,q)
                        independent=[sum((p for path,p in zip(paths,law) if path[t-1]==u),F(0))
                                     for u in range(2)]
                        assert actual==independent
                        enumerated+=1
            cost+=(k-t)*(joint[0][1]+joint[1][0])
            assert cost<=r/EPS
            max_cost=max(max_cost,cost);cases+=1
    return {'skeleton_bridge_cases':cases,'direct_path_marginal_checks':enumerated,
            'matrix_power_identities':25,'segment_cost_upper':str(max_cost),
            'uniform_budget':str(r/EPS)}


def tilt_controls():
    delta=F(1,4096);eta=delta/(1-delta);cases=0;flips=0
    for k in range(1,6):
        env=tuple((u%2,(u+1)%2) for u in range(k))
        for left,right in itertools.product(range(2),repeat=2):
            paths,p=bridge_law(k,left,right)
            _,tilt=bridge_law(k,left,right,delta=delta,env=env)
            tv=sum(abs(x-y) for x,y in zip(p,tilt))/2
            assert tv<=4*eta*k
            assert coupling_cost(paths,p,tilt)<=4*eta*k*k
            # Direct density-ratio gate includes the changed normalizer.
            ratio=((1-delta)/(1+delta))**(2*k)
            assert all(x/y>=ratio for x,y in zip(tilt,p))
            for u in range(k):
                changed=list(env);changed[u]=(1-env[u][0],env[u][1])
                _,other=bridge_law(k,left,right,delta=delta,env=tuple(changed))
                assert sum(abs(x-y) for x,y in zip(tilt,other))/2<=4*eta
                assert coupling_cost(paths,tilt,other)<=4*eta*k
                flips+=1
            cases+=1
    return {'tilted_block_cases':cases,'one_spatial_coordinate_changes':flips,
            'rational_tilt_log_budget':str(eta),'normalizers_checked':True}


def blocks(m,N,ell):
    if min(m,N,ell)<1 or any(not isinstance(v,int) for v in (m,N,ell)):
        raise ValueError('positive integer block dimensions required')
    return [(i,tuple(range(max(1,s),min(N,s+ell-1)+1)))
            for i in range(m) for s in range(2-ell,N+1)]


def incidence_controls():
    layouts=0;points=0;weighted=0
    eta=F(1,100000)
    for m,N,ell in itertools.product((1,2,4),(1,2,4,7,12),(1,2,3,5,8)):
        labels=blocks(m,N,ell)
        assert len(labels)==m*(N+ell-1)
        D=1/EPS+8*eta*ell*ell
        base=1+F(1,10*ell);R=base**ell
        w=lambda u:base**min(u,N+1-u)
        for i in range(m):
            for u in range(1,N+1):
                removal=0;temporal=0;spatial=0;cost=F(0)
                for j,B in labels:
                    if i==j and u in B:
                        removal+=1
                    if i==j and u in (B[0]-1,B[-1]+1):
                        temporal+=1;cost+=max(w(v) for v in B)*D
                    if abs(i-j)==1 and u in B:
                        spatial+=1;cost+=max(w(v) for v in B)*4*eta*ell
                assert removal==ell and temporal<=2 and spatial<=2*ell
                assert cost<=w(u)*R*(2/EPS+24*eta*ell*ell)
                points+=1;weighted+=1
            for u in (0,N+1):
                hits=[B for j,B in labels if j==i and u in (B[0]-1,B[-1]+1)]
                assert len(hits)==ell
                assert sum(max(w(v) for v in B)*D for B in hits)<=ell*R*D
                weighted+=1
        layouts+=1
    return {'clipped_layouts':layouts,'interior_incidence_checks':points,
            'weighted_outgoing_budget_checks':weighted}


def strip_law(m,N,left,right,q,delta):
    states=list(itertools.product(range(2),repeat=m*N));weights=[]
    K=power_kernel(1,q)
    for x in states:
        w=F(1)
        for i in range(m):
            col=(left[i],)+tuple(x[(u-1)*m+i] for u in range(1,N+1))+(right[i],)
            for a,b in zip(col,col[1:]):
                w*=K[a][b]
        for u in range(N):
            for i in range(m-1):
                w*=1+delta*(2*x[u*m+i]-1)*(2*x[u*m+i+1]-1)
        weights.append(w)
    return states,normalized(weights)


def conditional_paths(states,law,x,indices):
    outside=[j for j in range(len(x)) if j not in indices]
    ids=[k for k,s in enumerate(states) if all(s[j]==x[j] for j in outside)]
    return ids,normalized([law[k] for k in ids])


def joint_update_controls():
    m,N,ell=2,2,8;q=F(1,100);delta=F(1,2000000)
    eta=delta/(1-delta);R=F(21,20)
    alpha=2/(EPS*ell)+24*eta*ell
    assert R*alpha<1 and ((1-q)/(1+q))**2>=EPS
    labels=blocks(m,N,ell);M=len(labels);D=1/EPS+8*eta*ell*ell
    states,p=strip_law(m,N,(0,0),(0,0),q,delta)
    rows=[]
    for different in (False,True):
        _,other=strip_law(m,N,((1,1) if different else (0,0)),(0,0),q,delta)
        n=len(states);out=[[F(0) for _ in states] for _ in states]
        cache={};pair_checks=0
        for i,x in enumerate(states):
            for j,y in enumerate(states):
                weight=p[i]*other[j]/M
                for col,B in labels:
                    indices=tuple((u-1)*m+col for u in B)
                    key=(i,j,indices)
                    if key not in cache:
                        ix,px=conditional_paths(states,p,x,indices)
                        iy,py=conditional_paths(states,other,y,indices)
                        cache[key]=(ix,iy,coupling(px,py))
                    ix,iy,C=cache[key]
                    for a,ia in enumerate(ix):
                        for b,ib in enumerate(iy):
                            out[ia][ib]+=weight*C[a][b]
                    pair_checks+=1
        assert [sum(row) for row in out]==p
        assert [sum(out[i][j] for i in range(n)) for j in range(n)]==other
        old=sum((p[i]*other[j]*R*hamming(x,y)
                 for i,x in enumerate(states) for j,y in enumerate(states)),F(0))
        new=sum((out[i][j]*R*hamming(x,y)
                 for i,x in enumerate(states) for j,y in enumerate(states)),F(0))
        forcing=(m*ell*R*D if different else F(0))
        ceiling=(1-ell*(1-R*alpha)/M)*old+forcing/M
        assert new<=ceiling and new<old
        rows.append({'different_left_boundary':different,'pair_label_updates':pair_checks,
                     'both_marginals_exact':True,'weighted_mismatch_before':decimals(old),
                     'weighted_mismatch_after':decimals(new),'recurrence_ceiling':decimals(ceiling)})
    return {'states':len(states),'update_labels':M,'cases':rows}


def vacuum_and_limit_controls():
    T=[[F(9,8),F(1,8)],[F(1,8),F(3,4)]]
    power=[[F(1),F(0)],[F(0),F(1)]]
    for _ in range(5):
        power=multiply(power,T)
    assert multiply(T,power)==multiply(power,T)
    count=vacuum_controls(power,steps=9)
    h=[F(1),F(1)];ranges=[]
    for _ in range(12):
        h=[sum(a*b for a,b in zip(row,h)) for row in power]
        h=[x/sum(h) for x in h]
        ratios=[sum(a*b for a,b in zip(row,h))/h[i] for i,row in enumerate(T)]
        ranges.append(max(ratios)-min(ratios))
    assert all(b<a for a,b in zip(ranges,ranges[1:]))
    rows=[]
    for theta,J,R in CELLS:
        gamma=block_gate(theta,J,R)['gamma']
        for m,t in itertools.product((1,4,16),(F(1),F(6))):
            b=abs(theta)*(m-1);ell=ex(-b*t).lo
            q=ex(-gamma.lo*t).hi
            ratios=[]
            for n in (64,4096,262144):
                error=refinement_budget(m,theta,t,t/n)
                bound=gap_gate(ell,q,error)
                ratios.append(None if bound is None else decimals(bound))
            assert ratios[-1] is not None
            rows.append({'theta_abs':str(theta),'m':m,'t':str(t),
                         'limiting_ratio_upper':decimals(q),'transported_ratio_ceilings':ratios})
    return {'coarse_power_vacuum_checks':count,'fine_vacuum_ratio_ranges_decrease':True,
            'fixed_width_limit_controls':rows,'refinement_counts':[64,4096,262144]}


def negative_controls():
    # Omitting the future factor changes an actual conditioned path marginal.
    true=bridge_transition(0,1,5,1)
    false=power_kernel(5)[0]
    # A normalized tilted law can fall below the unnormalized tilt floor.
    base=[F(9,10),F(1,10)];w=[F(2),F(1)];tilt=normalized([x*y for x,y in zip(base,w)])
    # A fixed number of fine slices has no cutoff-independent physical length.
    small=F(1,1024);r=-((-6*small.denominator)//small.numerator)
    controls={
        'future_bridge_factor_omission_rejected':true!=false,
        'tilt_normalizer_omission_rejected':tilt[1]/base[1]<min(w),
        'non_strict_weighted_gate_rejected':block_gate(F(1,4096),8,F(64,41)) is None,
        'endpoint_window_has_no_strict_unweighted_margin':F(1,2)+840*F(1,1680)==1,
        'old_theta_one_sixteenth_not_promoted':all(block_gate(F(1,16),J,F(101,100)) is None for J in range(1,101)),
        'missing_temporal_tilt_would_admit_uncertified_cell':F(101,100)*(F(1,2)+280*F(1,1000))<1
            and block_gate(F(1,1000),5,F(101,100)) is None,
        'unclipped_block_coverage_rejected':sum(1 for s in range(1,3) if 1 in range(s,s+3))!=3,
        'single_spatial_neighbor_count_rejected':sum(abs(i-1)==1 for i in range(3))==2,
        'fixed_fine_slice_count_rejected':2*r/(EPS*40)>1,
        'constant_per_step_ratio_rejected':ex(-F(1,100)/1024).lo>F(99,100),
        'vacuum_denominator_omission_rejected':gap_gate(F(1),F(1,2),F(1,8))==F(5,7)>F(5,8),
        'invalid_block_dimensions_rejected':False}
    try:
        blocks(2,0,8)
    except ValueError:
        controls['invalid_block_dimensions_rejected']=True
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
    return {'certificate_type':'YM45_TEMPORAL_BLOCK_UNIFORM_GAP','verdict':'PASS',
            'runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'full_heat_and_parameters':parameter_controls(),'free_bridge_coupling':bridge_controls(),
            'tilted_conditionals':tilt_controls(),'clipped_block_incidence':incidence_controls(),
            'independent_joint_updates':joint_update_controls(),
            'vacuum_and_time_limit':vacuum_and_limit_controls(),'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_ALL_SOURCE_GAP_UNIFORM_IN_FINITE_WIDTH_AND_FINE_TIME_STEP__'
                           'SMALL_BRIDGE_WINDOW_WITHIN_DECLARED_HEAT_FUNCTIONAL_CHAIN__'
                           'FIXED_WIDTH_TIME_LIMIT_GAP__NATIVE_DICTIONARY_AND_4D_CONTINUUM_OPEN',
            'evidence_scope':{'general_proof':'written; not mechanically formalized',
                'finite_controls':'exact rational and outward intervals; binary fixtures are not SU(2)',
                'trajectory':'declared kappa(a)=theta*a; abs(theta)<1/1680; not AF derived',
                'uniformity':'rate gamma independent of finite width m and 0<a<=1; prefactors need not be',
                'time_limit':'all t>0, fixed finite width; no infinite-volume construction',
                'physical_clock_measure_and_4D_Clay':'not established'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM45 fresh certificate/pin mismatch')
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
    print('YM45 PASS',sha)
    print(json.dumps({'parameter_cells':len(cert['full_heat_and_parameters']['cells']),
                      'block_layouts':cert['clipped_block_incidence']['clipped_layouts'],
                      'refusal_groups':len(cert['negative_controls'])},sort_keys=True))


if __name__=='__main__':
    main()
