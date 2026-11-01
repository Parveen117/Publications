"""YM51: dependency audit and native heat-selection controls.

Python 3.12 only. Written proofs and exact controls have separate scopes.
Default/--check is read-only; --write changes this chapter's evidence only.
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
import ym50_native_reference as y

a=y.y49
p=y.y48
RESULT=HERE/'YM51_RESULT.json'
PIN=HERE/'EXPECTED_YM51.sha256'
SOURCES=HERE/'YM51_SOURCE_PINS.json'
LEDGER=HERE/'YM51_DEPENDENCY_LEDGER.json'
I3=a.identity(3)
ZERO3=(F(0),)*3
AXES=tuple(tuple(F(i==j) for i in range(3)) for j in range(3))
SLOTS=tuple((i,j) for i in range(3) for j in range(i,3))
SLOW=p.add(p.mul(p.COORD[2],p.COORD[2]),p.mul(p.COORD[3],p.COORD[3]),
           p.scale(y.RADIUS,F(-1,2)))


def exact(value):
    if isinstance(value,bool) or not isinstance(value,(int,F)):
        raise ValueError('exact rational input required')
    return F(value)


def protocol(records):
    out=[]
    for weight,vector in records:
        weight=exact(weight)
        if weight<0 or len(vector)!=3:
            raise ValueError('nonnegative weight and three coefficients required')
        out.append((weight,tuple(exact(v) for v in vector)))
    if not out or sum(w for w,_ in out)!=1:
        raise ValueError('protocol weights must sum to one')
    return tuple(out)


def moments(records):
    records=protocol(records)
    C=[[sum((w*v[i]*v[j] for w,v in records),F(0)) for j in range(3)]
       for i in range(3)]
    m2=sum((w*a.dot(v,v) for w,v in records),F(0))
    m4=sum((w*a.dot(v,v)**2 for w,v in records),F(0))
    assert a.psd(C) and sum(C[i][i] for i in range(3))==m2 and m4>=m2*m2
    return C,m2,m4


def eta_protocol(eta):
    eta=exact(eta)
    if not 0<=eta<=1:
        raise ValueError('eta must lie in [0,1]')
    weights=((3-2*eta)/9,eta/9,eta/9)
    return protocol([(F(2,3),ZERO3)]+[(w,tuple(3*v for v in axis))
                                    for w,axis in zip(weights,AXES)])


def tensor(C):
    if len(C)!=3 or any(len(row)!=3 for row in C):
        raise ValueError('three by three tensor required')
    C=[[exact(v) for v in row] for row in C]
    if C!=a.transpose(C) or not a.psd(C):
        raise ValueError('symmetric positive semidefinite tensor required')
    return C


def derivative(poly,v):
    return p.add(*(p.scale(y.native_generator(poly,i+1),v[i]) for i in range(3)))


def lap(poly,C):
    C=tensor(C)
    return p.scale(p.add(*(p.scale(y.native_generator(
        y.native_generator(poly,j+1),i+1),C[i][j])
        for i in range(3) for j in range(3))),-1)


@lru_cache(maxsize=None)
def generators(d):
    return tuple(a.transpose([y.vector(y.native_generator({m:F(1)},j),d)
                              for m in y.basis(d)]) for j in range(1,4))


def operator(d,C):
    C=tensor(C);D=generators(d);out=a.scale(a.identity(len(y.basis(d))),0)
    for i,j in itertools.product(range(3),repeat=2):
        out=a.add(out,a.scale(a.multiply(D[i],D[j]),-C[i][j]))
    return out


def rotation(g):
    if y.qmul(g,y.qdagger(g))!=y.ONE:
        raise ValueError('unit native quaternion required')
    cols=[y.qmul(y.qmul(g,(F(0),)+v),y.qdagger(g))[1:] for v in AXES]
    R=a.transpose(cols)
    assert a.multiply(a.transpose(R),R)==I3
    return R


def rotate_tensor(C,R):
    return a.multiply(a.multiply(R,tensor(C)),a.transpose(R))


def budget(d,t,n,m2,m4):
    if isinstance(d,bool) or not isinstance(d,int) or d<0:
        raise ValueError('nonnegative integer degree required')
    if isinstance(n,bool) or not isinstance(n,int) or n<1:
        raise ValueError('positive integer refinement required')
    t,m2,m4=map(exact,(t,m2,m4))
    if min(t,m2,m4)<0 or m4<m2*m2:
        raise ValueError('nonnegative time and consistent moments required')
    return F(d**4,96*n)*t*t*(m4+3*m2*m2)


def protocol_controls():
    tetra=protocol([(F(1,4),tuple(map(F,v))) for v in
                    ((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))])
    rows=[]
    records=[tetra]+[eta_protocol(eta) for eta in (F(0),F(1,100),F(1,4),F(1))]
    factor_checks=energy_checks=0
    for index,rec in enumerate(records):
        C,m2,m4=moments(rec)
        for d in range(5):
            L=operator(d,C);G=y.degree_data(d)['G']
            factored=a.scale(a.identity(len(L)),0)
            for w,v in rec:
                D=a.scale(factored,0)
                for coefficient,base in zip(v,generators(d)):
                    D=a.add(D,a.scale(base,coefficient))
                factored=a.add(factored,a.scale(a.multiply(D,D),-w))
            assert L==factored;factor_checks+=1
            GL=a.multiply(G,L)
            assert GL==a.transpose(GL) and a.psd(GL);energy_checks+=1
        rows.append({'protocol':index,'C':[[str(v) for v in row] for row in C],
                     'm2':str(m2),'m4':str(m4)})
    assert moments(tetra)[0]==I3
    assert budget(4,F(2),32,F(3),F(9))==y.heat_budget(4,F(2),32)
    return {'protocols':rows,'factorization_checks':factor_checks,
            'all_source_energy_checks':energy_checks,'YM50_bound_recovered':True}


def invariance_controls():
    # The native reference is unchanged for these exact finite unit turns.
    checks=0
    g=y.qtuple(y.fabric.rational_unit(F(2,3),F(-1,4),F(1,5)))
    for d in range(5):
        for m in y.basis(d):
            poly={m:F(1)}
            assert y.phi(y.substitute(poly,y.qmatrix(g)))==y.phi(poly);checks+=1
    bracket=0
    for d in range(1,5):
        D=generators(d)
        for i,j,k in ((0,1,2),(1,2,0),(2,0,1)):
            assert a.add(a.multiply(D[i],D[j]),a.multiply(D[j],D[i]),-1)==D[k]
            bracket+=1
    return {'same_reference_turn_checks':checks,'unchanged_native_bracket_checks':bracket}


def covariance_controls():
    qcycle=(F(1,2),)*4
    gs=[(F(0),)+v for v in AXES]+[qcycle,
        y.qtuple(y.fabric.rational_unit(F(1,3),F(2,5),F(-1,4)))]
    Cs=[I3,a.diag((F(5,2),F(1,4),F(1,4))),
        moments(protocol([(F(1,2),(1,2,0)),(F(1,2),(0,1,1))]))[0]]
    checks=0
    for C,g,d in itertools.product(Cs,gs,range(4)):
        R=rotation(g);rotated=rotate_tensor(C,R)
        U=y.action_matrix(d,y.qmatrix(g))
        inverse=y.action_matrix(d,y.qmatrix(y.qdagger(g)))
        assert a.multiply(a.multiply(inverse,operator(d,C)),U)==operator(d,rotated)
        checks+=1
    # Nullspace of the four fixed-law symmetry equations, on six symmetric slots.
    bases=[]
    for i,j in SLOTS:
        B=a.scale(I3,0);B[i][j]=B[j][i]=F(1);bases.append(B)
    equations=[]
    for g in gs[:4]:
        R=rotation(g)
        columns=[]
        for B in bases:
            difference=a.add(a.multiply(a.multiply(R,B),a.transpose(R)),B,-1)
            columns.append([difference[i][j] for i,j in SLOTS])
        equations.extend(a.transpose(columns))
    gram=a.multiply(a.transpose(equations),equations)
    assert y.y37.inertia(gram)==(5,0,1)
    unit_slots=[F(i==j) for i,j in SLOTS]
    assert a.apply(equations,unit_slots)==[F(0)]*len(equations)
    C=Cs[1];R=rotation(qcycle)
    assert rotate_tensor(C,R)!=C
    assert operator(2,rotate_tensor(C,R))!=operator(2,C)
    return {'covariant_operator_checks':checks,'fixed_law_constraint_inertia':[5,0,1],
            'covariance_does_not_imply_isotropy':True}


def selection_controls():
    rows=[];linear=quadratic=0
    for eta in (F(0),F(1,1000),F(1,100),F(1,4),F(1,2),F(1)):
        C,m2,m4=moments(eta_protocol(eta))
        for x in p.COORD:
            assert lap(x,C)==p.scale(x,F(3,4));linear+=1
        assert lap(SLOW,C)==p.scale(SLOW,2*eta);quadratic+=1
        assert y.phi(SLOW)==0 and y.phi(p.mul(SLOW,SLOW))==F(1,12)
        inertia=y.y37.inertia(a.multiply(y.degree_data(2)['G'],operator(2,C)))
        assert inertia[1]==0
        rows.append({'eta':str(eta),'same_linear_rate':'3/4','quadratic_rate':str(2*eta),
                     'quadratic_source_norm_squared':'1/12',
                     'full_axis_support':eta>0,'degree_two_energy_inertia':list(inertia)})
    assert rows[1]['full_axis_support'] and F(rows[1]['quadratic_rate'])<F(1,100)
    return {'linear_blindness_checks':linear,'quadratic_witness_checks':quadratic,
            'family':rows,'uniform_rate_over_this_family':'NO; upper bound 2*eta'}


def recovery_controls():
    recs=[eta_protocol(F(1,5)),
          protocol([(F(1,3),(1,2,3)),(F(2,3),(2,-1,1))]),
          protocol([(F(1,4),(1,0,0)),(F(1,4),(0,1,0)),
                    (F(1,4),(0,0,1)),(F(1,4),(1,1,1))])]
    checks=0
    for rec in recs:
        C,_,_=moments(rec)
        for i,j in SLOTS:
            poly=p.mul(p.COORD[i+1],p.COORD[j+1])
            assert -2*y.evaluate(lap(poly,C),y.ONE)==C[i][j];checks+=1
    # A positive off-diagonal perturbation is invisible to diagonal-only probes.
    altered=a.add(I3,[[F(0),F(1,4),F(0)],[F(1,4),F(0),F(0)],[F(0)]*3])
    assert a.psd(altered)
    for i in range(3):
        poly=p.mul(p.COORD[i+1],p.COORD[i+1])
        assert y.evaluate(lap(poly,altered),y.ONE)==y.evaluate(lap(poly,I3),y.ONE)
    cross=p.mul(p.COORD[1],p.COORD[2])
    assert y.evaluate(lap(cross,altered),y.ONE)==F(-1,8)
    assert y.evaluate(lap(cross,I3),y.ONE)==0
    return {'quadratic_tensor_recovery_checks':checks,'minimum_fixed_linear_channels':6,
            'diagonal_only_bank_has_positive_counterexample':True}


def count_controls():
    # eta=1/4: denominator 72 after splitting each non-idle weight into two signs.
    eta=F(1,4);rec=eta_protocol(eta)
    letters=[(rec[0][0],y.ONE)]
    for index,(weight,_) in enumerate(rec[1:]):
        for turn in y.turns((index+1,)):
            letters.append((weight/2,turn))
    labels=[q for w,q in letters for _ in range(int(72*w))]
    assert len(labels)==72 and sum(w for w,_ in letters)==1
    samples=[p.ONE,*p.COORD,SLOW,p.mul(p.COORD[1],p.COORD[2])]
    count=0
    for poly in samples:
        weighted=sum((w*y.evaluate(poly,q) for w,q in letters),F(0))
        raw=sum((y.evaluate(poly,q) for q in labels),F(0))/72
        assert weighted==raw;count+=1
        second=sum((u*v*y.evaluate(poly,y.qmul(g,h))
                    for u,g in letters for v,h in letters),F(0))
        raw2=sum((y.evaluate(poly,y.qmul(g,h)) for g in labels for h in labels),F(0))/72**2
        assert second==raw2;count+=1
    return {'independent_labelled_count_checks':count,'one_step_labels':72,
            'two_step_labels':72**2}


def refinement_controls():
    rows=[]
    for eta,t,n,d in itertools.product((F(0),F(1,100),F(1,4),F(1)),
                                       (F(1,4),F(1),F(2)),(32,128,512),(1,2)):
        if d==1:
            cosine=y.cosine_square(F(9,2)*t/n)
            step=(y.Iv(F(2,3))+y.Iv(F(1,3))*cosine);rate=F(3,4)
        else:
            weight=2*eta/9;cosine=y.cosine_square(18*t/n)
            step=y.Iv(1-weight)+y.Iv(weight)*cosine;rate=2*eta
        approximate=y.y44.ivpower(y._r(step),n);target=p.exp_negative(rate*t)
        difference=approximate-target;error=max(abs(difference.lo),abs(difference.hi))
        bound=budget(d,t,n,F(3),F(27))
        assert error<=bound
        rows.append({'eta':str(eta),'time':str(t),'n':n,'degree':d,
                     'error_upper':str(error),'bound':str(bound)})
    return {'outward_scalar_refinement_checks':len(rows),'cases':rows}


def clock_and_micro_controls():
    first=protocol([(1,(1,0,0))])
    second=protocol([(F(3,4),ZERO3),(F(1,4),(2,0,0))])
    C,m2,m4=moments(first);B,n2,n4=moments(second)
    assert C==B and m2==n2==1 and (m4,n4)==(1,4)
    assert (n4-m4)/384==F(1,128)
    scale_checks=0
    for c,d in itertools.product((F(1,4),F(1),F(7,3)),range(4)):
        assert operator(d,a.scale(I3,c))==a.scale(operator(d,I3),c);scale_checks+=1
    ratios=[]
    for c,theta,step in ((F(1,4),F(1,10000),F(1,8)),
                         (F(2),F(1,4096),F(1,16)),(F(3),F(-1,9000),F(1,10))):
        b=c*step;relative=theta/c
        assert b*relative==step*theta
        ratios.append({'c':str(c),'theta':str(theta),'a':str(step),
                       'b':str(b),'relative_interaction':str(relative)})
    return {'same_second_moment_different_microsteps':True,
            'microstep_fourth_coefficient_difference':'1/128',
            'clock_rescaling_checks':scale_checks,'interaction_rescalings':ratios,
            'physical_clock_selection':'OPEN'}


def validate_ledger(ledger):
    statuses={'NATIVE_DERIVED_UNDER_CONTRACT','PROVED_REPRESENTATION',
              'DECLARED_SELECTION','OPEN_EXTENSION'}
    nodes=ledger.get('nodes',[]);by_id={row['id']:row for row in nodes}
    if len(by_id)!=len(nodes):
        raise ValueError('duplicate dependency id')
    sources=ledger.get('sources',{})
    for row in nodes:
        if row.get('status') not in statuses or not row.get('contract') or not row.get('boundary'):
            raise ValueError('missing status or contract boundary')
        if not row.get('sources') or any(s not in sources for s in row['sources']):
            raise ValueError('missing source')
        if any(parent not in by_id for parent in row['parents']):
            raise ValueError('missing parent')
    visiting=set();visited=set()
    def visit(key):
        if key in visiting:
            raise ValueError('cyclic dependency')
        if key in visited:return
        visiting.add(key)
        for parent in by_id[key]['parents']:visit(parent)
        visiting.remove(key);visited.add(key)
    for key in by_id:visit(key)
    for key in ledger['must_remain_selections']:
        if key not in by_id or by_id[key]['status']!='DECLARED_SELECTION':
            raise ValueError('physical or protocol selection promoted: '+key)
    for key in ledger['must_remain_open']:
        if key not in by_id or by_id[key]['status']!='OPEN_EXTENSION':
            raise ValueError('open obligation promoted: '+key)
    for key in ('raw_order_defect','connection_curvature','observer_readout'):
        if not by_id[key].get('target_type'):
            raise ValueError('curvature target type absent')
    if len({by_id[key]['target_type'] for key in
            ('raw_order_defect','connection_curvature','observer_readout')})!=3:
        raise ValueError('distinct curvature targets collapsed')
    return {'nodes':len(nodes),'source_groups':len(sources),
            'status_counts':{s:sum(row['status']==s for row in nodes) for s in sorted(statuses)},
            'acyclic':True}


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
    ledger=json.loads(LEDGER.read_text())
    return {'certificate_type':'YM51_NATIVE_TOOL_AUDIT_AND_HEAT_SELECTION',
            'verdict':'PASS','runtime_policy':'Python 3.12 only',
            'input_sha256':source_checks(),'dependency_audit':validate_ledger(ledger),
            'protocol_energy':protocol_controls(),'native_frame':invariance_controls(),
            'covariance_and_isotropy':covariance_controls(),
            'nonselection_witness':selection_controls(),'tensor_observer':recovery_controls(),
            'counted_protocol':count_controls(),'refinement':refinement_controls(),
            'clock_and_microstructure':clock_and_micro_controls(),
            'claim_status':'SPECIFIED_NATIVE_PROTOCOL_CLASS_AND_SELECTION_BOUNDARY_PROVED__'
                'PHYSICAL_ISOTROPY_CLOCK_INTERACTION_AND_4D_NOT_SELECTED',
            'evidence_scope':{
                'proof':'written general arguments; not mechanically formalized or expert-certified',
                'audit':'current YM43-50 dependencies; not the entire framework',
                'isotropy':'conditional on fixed-law active invariance; covariance alone is insufficient',
                'decay':'free alternative protocols; no new interacting gap theorem',
                'open':'physical state/action/clock, general UGD, actual row closure, NCG quantum measure, 4D/AF/Clay/QG'}}


def check(cert,result=RESULT,pin=PIN):
    sha=y.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=sha:
        raise ValueError('YM51 fresh certificate/pin mismatch')
    return sha


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true');group.add_argument('--check',action='store_true')
    args=parser.parse_args();cert=run();sha=y.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n');PIN.write_text(sha+'\n')
    else:check(cert)
    print('YM51 PASS',sha)
    print(json.dumps({'dependency_nodes':cert['dependency_audit']['nodes'],
        'covariance_checks':cert['covariance_and_isotropy']['covariant_operator_checks'],
        'linear_blindness_checks':cert['nonselection_witness']['linear_blindness_checks'],
        'refinement_checks':cert['refinement']['outward_scalar_refinement_checks'],
        'physical_selection':'OPEN'},sort_keys=True))


if __name__=='__main__':main()
