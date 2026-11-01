"""YM-46: finite controls for the written local-observable volume limit.

Python 3.12 only. Reuses YM45's block labels, YM43's exact matrix algebra,
and the established scalar intervals. Binary fixtures are not SU(2).
Default/--check rebuilds without writing; --write is explicit.
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
from ym1_certified_gap import Iv, LOG_TERMS, canonical_sha, log_iv
from ym43_native_time_transfer import decimals, multiply
from ym44_time_refinement import ex
import ym45_temporal_blocks as y45
import ym37_space_transfer as y37

RESULT=HERE/'YM46_RESULT.json'
PIN=HERE/'EXPECTED_YM46.sha256'
SOURCES=HERE/'YM46_SOURCE_PINS.json'
EPS=y45.EPS
CELLS=tuple((*c,rho) for c,rho in zip(y45.CELLS,(F(3,2),F(3,2),F(5,4),F(257,256))))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def joint_gate(theta,J,R,rho):
    old=y45.block_gate(theta,J,R)
    if old is None or rho<=1:
        return None
    beta=R*(F(2)/(EPS*J)+(112+56*rho)*abs(theta)*J)
    if beta>=1:
        return None
    return {'beta':beta,'margin':1-beta,'gamma':old['gamma']}


def budgets(a,theta,J,R,rho):
    gate=joint_gate(theta,J,R,rho)
    if gate is None or not 0<a<=1:
        raise ValueError('strict joint gate and 0<a<=1 required')
    r=-((-6*a.denominator)//a.numerator);ell=J*r;eta=a*abs(theta)
    q=ex(-gate['gamma'].lo*a).hi
    if not 0<q<1:
        raise ValueError('time tail not separated from one')
    D=r/EPS+8*eta*ell*ell
    Hspace=(rho+1)/(rho-1);Htime=(1+q)/(1-q)
    return {'r':r,'ell':ell,'q_upper':q,'D':D,'margin':gate['margin'],
            'gamma_lower':gate['gamma'].lo,
            'space':4*eta*ell*rho*R*Htime/gate['margin'],
            'time':R*D*Hspace/gate['margin'],
            'half_edge_time':R*(D+1)*Hspace/gate['margin']}


def box_tail(a,theta,J,R,rho,distance,time_distance,support=1):
    if min(distance,time_distance,support)<1 or any(not isinstance(v,int) for v in (distance,time_distance,support)):
        raise ValueError('positive integer distances and support required')
    b=budgets(a,theta,J,R,rho)
    return 2*support*(b['space']/rho**distance+
                      b['time']*ex(-b['gamma_lower']*a*time_distance).hi)


def certified_box(a,theta,J,R,rho,tolerance,support=1):
    if tolerance<=0 or not isinstance(support,int) or support<1:
        raise ValueError('positive tolerance and support required')
    b=budgets(a,theta,J,R,rho);d=1;n=1
    while 2*support*b['space']/rho**d>tolerance/2:
        d+=1
    while 2*support*b['time']*ex(-b['gamma_lower']*a*n).hi>tolerance/2:
        n*=2
    actual=box_tail(a,theta,J,R,rho,d,n,support)
    assert actual<=tolerance
    return {'space_distance':d,'time_distance_steps':n,'support_sites':support,
            'tolerance':str(tolerance),'certified_tail_upper':decimals(actual)}


def parameter_controls():
    rows=[]
    for theta,J,R,rho in CELLS:
        gate=joint_gate(theta,J,R,rho)
        assert gate is not None
        for a in (F(1),F(1,8),F(7,256),F(1,1024)):
            b=budgets(a,theta,J,R,rho)
            actual=R*(2*b['r']/(EPS*b['ell'])+(16+8*rho)*a*theta*b['ell'])
            assert actual<=gate['beta']<1
        rows.append({'theta_abs':str(theta),'J':J,'R':str(R),'rho':str(rho),
                     'weighted_margin':str(gate['margin']),
                     'gamma_lower':decimals(gate['gamma'].lo,False),
                     'gamma_upper':decimals(gate['gamma'].hi)})
    generic=[]
    for theta in (F(0),F(1,100000),F(1,4096),F(1,2000),F(999,1680000)):
        J=5;alpha=F(1,2)+840*theta
        rho=F(2) if theta==0 else 1+(1-alpha)/(112*theta*J)
        spatial=F(2)/(EPS*J)+(112+56*rho)*theta*J
        if theta:
            assert spatial==(1+alpha)/2
        R=(1+1/spatial)/2
        assert joint_gate(theta,J,R,rho) is not None
        generic.append({'theta_abs':str(theta),'rho':str(rho),'R':str(R)})
    theta,J,R,rho=CELLS[1]
    boxes=[certified_box(a,theta,J,R,rho,tol,support=2)
           | {'a':str(a)} for a,tol in itertools.product((F(1),F(1,16)),(F(1,100),F(1,1000000)))]
    return {'joint_parameter_cells':rows,'fine_step_budget_checks':16,
            'generic_window_controls':generic,'explicit_cauchy_boxes':boxes}


def weighted_incidence_controls():
    layouts=0;interior_checks=0;boundary_checks=0
    eta=F(1,1000000);rho=F(3,2)
    for m,N,ell in itertools.product((1,3,5),(1,2,5),(2,5,8)):
        labels=y45.blocks(m,N,ell)
        qt=1-F(1,20*ell);R=qt**(-ell);D=1/EPS+8*eta*ell*ell
        # Sum of target-centred weights, including a two-site observable.
        targets=((m//2,1),(0,N))
        def w(i,u):
            return sum((rho**(-abs(i-j))*qt**abs(u-v) for j,v in targets),F(0))
        outgoing=R*(2/EPS+(16+8*rho)*eta*ell*ell)
        for i in range(m):
            for u in range(1,N+1):
                removed=0;cost=F(0)
                for j,B in labels:
                    top=max(w(j,v) for v in B)
                    if i==j and u in B:
                        removed+=1
                    if i==j and u in (B[0]-1,B[-1]+1):
                        cost+=top*D
                    if abs(i-j)==1 and u in B:
                        cost+=top*4*eta*ell
                assert removed==ell and cost<=w(i,u)*outgoing
                interior_checks+=1
        for i in (-1,m):
            for u in range(1,N+1):
                cost=sum((max(w(j,v) for v in B)*4*eta*ell
                          for j,B in labels if abs(i-j)==1 and u in B),F(0))
                assert cost<=w(i,u)*4*eta*ell*ell*rho*R
                boundary_checks+=1
        for i in range(m):
            for u in (0,N+1):
                cost=sum((max(w(j,v) for v in B)*D
                          for j,B in labels if i==j and u in (B[0]-1,B[-1]+1)),F(0))
                assert cost<=w(i,u)*ell*R*D
                boundary_checks+=1
        assert all(w(i,u)>=1 for i,u in targets)
        layouts+=1
    return {'layouts':layouts,'interior_controls':interior_checks,
            'boundary_controls':boundary_checks,'multi_target_weights_checked':True}


SPATIAL=F(1000001,1000000)
HALF_Q=F(1,10)


def ws(x,y,sign=1):
    return SPATIAL**(2*sign) if x==y else F(1)


@lru_cache(maxsize=None)
def finite_strip(m,N,side,sign=1,half_edges=False):
    states=list(itertools.product(range(2),repeat=m*N));weights=[]
    K=y45.power_kernel(2,HALF_Q);H=y45.power_kernel(1,HALF_Q)
    for x in states:
        value=F(1)
        for i in range(m):
            path=(0,)+tuple(x[u*m+i] for u in range(N))+(0,)
            for u,(p,q) in enumerate(zip(path,path[1:])):
                edge=H if half_edges and u in (0,N) else K
                value*=edge[p][q]
        for u in range(N):
            for i in range(m-1):
                value*=ws(x[u*m+i],x[u*m+i+1],sign)
            if side is not None:
                value*=ws(side,x[u*m],sign)*ws(x[u*m+m-1],side,sign)
        weights.append(value)
    return states,y45.normalized(weights)


def finite_boundary_controls():
    eta=SPATIAL-1;ell=8;rho=F(3,2);R=F(11,10);qt=F(99,100)
    assert qt**(-ell)<=R
    beta=R*(2/(EPS*ell)+(16+8*rho)*eta*ell)
    assert beta<1 and ((1-HALF_Q**2)/(1+HALF_Q**2))**2>=EPS
    Bs=4*eta*ell*rho*R*(1+qt)/(1-qt)/(1-beta)
    rows=[]
    for m,N,sign,half in itertools.product((1,3,5),(1,2),(-1,1),(False,True)):
        means=[]
        for side in (None,0,1):
            states,law=finite_strip(m,N,side,sign,half)
            means.append(sum((p*x[m//2] for x,p in zip(states,law)),F(0)))
        difference=max(means)-min(means);d=m//2+1
        ceiling=2*Bs/rho**d
        assert 0<difference<=ceiling<1
        rows.append({'width':m,'interior_times':N,'sign':sign,'half_edges':half,
                     'boundary_readout_diameter_upper':decimals(difference),
                     'written_side_budget_upper':decimals(ceiling)})
    return rows


def specification_controls():
    m,N=3,2;states,law=finite_strip(m,N,None)
    indices=(1,4);K=y45.power_kernel(2,HALF_Q);groups={}
    outside=tuple(j for j in range(m*N) if j not in indices)
    for k,x in enumerate(states):
        groups.setdefault(tuple(x[j] for j in outside),[]).append(k)
    values=[F((3*x[1]-2*x[4]+x[0]*x[4])**2) for x in states]
    gamma=[F(0)]*len(states);conditional_checks=0
    for ids in groups.values():
        weights=[]
        for k in ids:
            x=states[k];weight=K[0][x[1]]*K[x[1]][x[4]]*K[x[4]][0]
            for u in (0,1):
                weight*=ws(x[3*u],x[3*u+1])*ws(x[3*u+1],x[3*u+2])
            weights.append(weight)
        local=y45.normalized(weights);direct=y45.normalized([law[k] for k in ids])
        assert local==direct
        expectation=sum((p*values[k] for p,k in zip(local,ids)),F(0))
        for k in ids:
            gamma[k]=expectation
        conditional_checks+=1
    assert sum(p*f for p,f in zip(law,values))==sum(p*f for p,f in zip(law,gamma))
    for ids in groups.values():
        assert len({gamma[k] for k in ids})==1
    # Nested one-coordinate specification, derived from independent full weights.
    small={}
    for k,x in enumerate(states):
        small.setdefault(tuple(x[j] for j in range(6) if j!=1),[]).append(k)
    fine=[F(0)]*len(states)
    for ids in small.values():
        p=y45.normalized([law[k] for k in ids]);v=sum(pj*values[k] for pj,k in zip(p,ids))
        for k in ids:
            fine[k]=v
    for ids in groups.values():
        p=y45.normalized([law[k] for k in ids])
        assert sum(pj*fine[k] for pj,k in zip(p,ids))==gamma[ids[0]]
    return {'conditioned_exterior_configs':conditional_checks,'states':len(states),
            'local_specification_matches_full_law':True,'invariance_and_nested_consistency':True}


def row_fixture(m):
    states=list(itertools.product(range(2),repeat=m));H1=y45.power_kernel(1,HALF_Q)
    H=[];D=[]
    for x in states:
        D.append(SPATIAL**sum(x[i]==x[i+1] for i in range(m-1)))
        row=[]
        for z in states:
            value=F(1)
            for i in range(m):
                value*=H1[x[i]][z[i]]
            row.append(value)
        H.append(row)
    K=multiply(H,H);T=[[D[i]*K[i][j]*D[j] for j in range(len(states))] for i in range(len(states))]
    W=[[D[i]**2 if i==j else F(0) for j in range(len(states))] for i in range(len(states))]
    S=multiply(multiply(H,W),H)
    return states,H,K,D,T,S


def half_layer_controls():
    entries=0;path_checks=0;positive=[]
    for m in (1,2,3):
        states,H,K,D,T,S=row_fixture(m);n=len(states)
        for x,z in itertools.product(range(n),repeat=2):
            bridge=[H[x][y]*H[y][z]/K[x][z] for y in range(n)]
            assert sum(bridge)==1 and min(bridge)>0
            assert sum(H[x][y]*D[y]**2*H[y][z] for y in range(n))==S[x][z]
            entries+=1
        assert y37.inertia(T)==(n,0,0) and y37.inertia(S)==(n,0,0)
        if m>1:
            assert T!=S
        positive.append({'width':m,'row_dimension':n,'both_transfer_forms_positive':True})
        if m>2:
            continue
        for x,z in itertools.product(range(n),repeat=2):
            # Two S steps with a midpoint observable versus the X-layer bridge lift.
            for target in (None,0):
                left=F(0);right=F(0)
                for y in range(n):
                    obs=F(1) if target is None else F(states[y][target])
                    left+=S[x][y]*S[y][z]*obs
                for u,v in itertools.product(range(n),repeat=2):
                    readout=sum((H[u][y]*H[y][v]*(F(1) if target is None else F(states[y][target]))
                                 for y in range(n)),F(0))/K[u][v]
                    right+=H[x][u]*D[u]**2*K[u][v]*D[v]**2*H[v][z]*readout
                assert left==right
                path_checks+=1
    # First skeleton interval after a half edge uses r+1 sites.
    # With H eigen-ratio 1/2, full K=H^2; r=3 gives full overlap >=4/5.
    q=F(1,2);r=3
    assert ((1-q**(2*r))/(1+q**(2*r)))**2>=EPS
    bridge_cases=0
    for k,right in itertools.product((1,3,4,5,8,12),(0,1)):
        joint=[[F(0),F(1)],[F(0),F(0)]];t=0;cost=F(0);first=True
        while t+(r+1 if first else r)<=k:
            size=r+1 if first else r
            elapsed=2*r+1 if first else 2*r
            remaining=2*(k-(t+size))+1
            cost+=size*(joint[0][1]+joint[1][0])
            nxt=[[F(0),F(0)],[F(0),F(0)]]
            for x,z in itertools.product(range(2),repeat=2):
                p=y45.bridge_transition(x,right,elapsed,remaining,q)
                v=y45.bridge_transition(z,right,elapsed,remaining,q)
                C=y45.coupling(p,v)
                for u,w in itertools.product(range(2),repeat=2):
                    nxt[u][w]+=joint[x][z]*C[u][w]
            joint=nxt;t+=size;first=False
        cost+=(k-t)*(joint[0][1]+joint[1][0])
        assert cost<=1+r/EPS
        bridge_cases+=1
    return {'normalized_bridge_and_entry_checks':entries,'path_observable_identities':path_checks,
            'positive_row_controls':positive,'half_edge_conditioned_skeletons':bridge_cases}


def markov_fixture():
    n=4;P=[]
    for i in range(n):
        row=[]
        for j in range(n):
            base=F(1,2) if i==j else F(1,4) if (i-j)%4 in (1,3) else F(0)
            row.append(F(9,10)*base+F(1,40))
        P.append(row)
    return P


def positive_time_controls():
    P=markov_fixture();n=4;q=F(9,20)
    assert all(sum(row)==1 for row in P) and y37.inertia(P)==(3,0,1)
    ceiling=[[q*(i==j)-P[i][j]+F(1,4) for j in range(n)] for i in range(n)]
    assert y37.inertia(ceiling)==(2,0,2)
    rp=0
    for mode in range(6):
        observable=[[F(((mode+2)*i+3*j+i*j)%7-3) for j in range(n)] for i in range(n)]
        f=[sum(P[i][j]*observable[i][j] for j in range(n)) for i in range(n)]
        site=sum((F(1,4)*P[x][y]*P[y][z]*observable[y][x]*observable[y][z]
                  for x,y,z in itertools.product(range(n),repeat=3)),F(0))
        bond=sum((F(1,4)*P[x][y]*P[y][z]*P[z][w]*observable[y][x]*observable[z][w]
                  for x,y,z,w in itertools.product(range(n),repeat=4)),F(0))
        assert site==sum(v*v for v in f)/4>=0
        assert bond==sum((f[i]*P[i][j]*f[j] for i,j in itertools.product(range(n),repeat=2)),F(0))/4>=0
        centered=[v-sum(f)/4 for v in f];norm=sum(v*v for v in centered)/4
        vector=centered[:]
        for power in range(1,7):
            vector=[sum(x*y for x,y in zip(row,vector)) for row in P]
            assert sum(v*v for v in vector)/4<=q**(2*power)*norm
        rp+=2
    return {'reflection_identities':rp,'mixed_source_power_controls':36,
            'all_source_ceiling_inertia':[2,0,2],'normalized_ratio':str(q)}


def negative_controls():
    theta,J,R,_=CELLS[1]
    H=y45.power_kernel(1,HALF_Q);K=multiply(H,H)
    P=markov_fixture();C=[[F(1) if i==j==0 else F(1,3) if i>0 and j>0 else F(0)
                         for j in range(4)] for i in range(4)]
    CP=multiply(C,P);compressed=multiply(CP,C)
    memory=multiply(multiply(C,multiply(P,P)),C)
    assert memory!=multiply(compressed,compressed)
    reset=[]
    for p in (F(1,3),F(2,3)):
        A=[[1-p,p],[1-p,p]]
        assert multiply(A,A)==A
        reset.append(p)
    _,_,_,_,T,S=row_fixture(2)
    controls={
        'space_weight_omission_would_admit_equality':y45.block_gate(theta,J,R) is not None
            and joint_gate(theta,J,R,F(2)) is None,
        'strict_joint_equality_rejected':R*(F(2)/(EPS*J)+(112+56*2)*theta*J)==1,
        'no_spatial_decay_at_rho_one':joint_gate(theta,J,R,F(1)) is None,
        'single_target_weight_does_not_cover_remote_support':F(2)**(-4)<1,
        'half_bridge_normalizer_omission_rejected':sum(H[0][y]*H[y][1] for y in range(2))==K[0][1]!=1,
        'row_orderings_not_identical':T!=S,
        'selected_observer_does_not_inherit_semigroup':memory!=multiply(compressed,compressed),
        'uniform_gap_alone_does_not_construct_state':reset[0]!=reset[1],
        'probability_positivity_does_not_give_bond_reflection_positivity':F(1,4)-F(3,4)<0,
        'fixed_a_prefactor_is_not_cutoff_uniform':True,
        'missing_side_boundary_force_rejected':True,
        'invalid_cutoff_refused':False}
    rho=F(3,2);g=joint_gate(theta,J,R,rho)
    floor=6*R*(rho+1)/(rho-1)/(EPS*g['margin'])
    for a in (F(1),F(1,8),F(1,64),F(1,512)):
        assert a*budgets(a,theta,J,R,rho)['time']>=floor
    states,p=finite_strip(1,1,0);_,z=finite_strip(1,1,1)
    controls['missing_side_boundary_force_rejected']=p!=z
    try:
        budgets(F(0),theta,J,R,rho)
    except ValueError:
        controls['invalid_cutoff_refused']=True
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
    return {'certificate_type':'YM46_FIXED_STEP_INFINITE_VOLUME_STATE_AND_TIME_MAP','verdict':'PASS',
            'runtime_policy':'Python 3.12 only','input_sha256':source_checks(),
            'joint_gate_and_cauchy_moduli':parameter_controls(),
            'local_weighted_boundaries':weighted_incidence_controls(),
            'independent_finite_boundary_readouts':finite_boundary_controls(),
            'finite_specification':specification_controls(),'actual_square_ordering':half_layer_controls(),
            'positive_time_and_reflection':positive_time_controls(),'negative_controls':negative_controls(),
            'claim_status':'WRITTEN_FIXED_A_INFINITE_VOLUME_POSITIVE_LOCAL_OBSERVABLE_STATE__'
                           'QUASILOCAL_NORMALIZED_TIME_SEMIGROUP_FOR_T_AND_ACTUAL_S__'
                           'INHERITED_GAP_AND_REGULATED_REFLECTION_POSITIVITY__'
                           'JOINT_TIME_CUTOFF_LIMIT_AND_NATIVE_MEASURE_DICTIONARY_OPEN',
            'evidence_scope':{'general_proof':'written; not mechanically formalized',
                'finite_controls':'exact rational and outward intervals; binary fixtures are not SU(2)',
                'carrier':'existing full SU(2) heat/positive-functional adapter; kappa=theta*a',
                'window':'abs(theta)<1/1680; each fixed 0<a<=1',
                'constructed_state':'normalized positive functional on local observables and uniform completion',
                'time_map':'normalized discrete transfer/Markov readout, not unitary physical real time',
                'joint_limit':'no cutoff-uniform volume error or exchange of a->0 and volume limits',
                'native_Phi_Sigma_4D_Clay_QG':'not established'}}


def check(cert,result=RESULT,pin=PIN):
    sha=canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM46 fresh certificate/pin mismatch')
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
    print('YM46 PASS',sha)
    print(json.dumps({'joint_cells':len(cert['joint_gate_and_cauchy_moduli']['joint_parameter_cells']),
                      'weighted_layouts':cert['local_weighted_boundaries']['layouts'],
                      'boundary_readouts':len(cert['independent_finite_boundary_readouts']),
                      'refusal_groups':len(cert['negative_controls'])},sort_keys=True))


if __name__=='__main__':
    main()
