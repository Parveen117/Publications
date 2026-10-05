"""YM86: smooth trace measurements and energy-constrained response, Python 3.12.

Exact finite controls supplement written proofs. Weighted-Haar and matrix
fixtures are not interacting vacuum estimates. Default/--check are read-only;
--write changes only this stage's result and pin.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import isqrt
from pathlib import Path
import sys

import ym85_exterior_frames_energy as previous

y, p, poly = previous, previous.p, previous.previous
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULT = HERE/'YM86_RESULT.json'
PIN = HERE/'EXPECTED_YM86.sha256'
SOURCES = HERE/'YM86_SOURCE_MANIFEST.json'
BASE_COMMIT = '6088c2344325829949409d4418229a436604b11a'
BASE_TREE = 'd704459d5a829c830fb0ecac0cee2733003c73cd'
PARENT_HASH = '4bce553a56aaceb74da3968c681ead609b6a1e56cfea974ad8f678b9a55c353d'
MODEL = dict(
    smooth_seven_readings_separate_three_holonomy_orbits=True,
    seven_linear_markers_span_all_physical_observables=False,
    pair_trace_data_alone_determine_orientation=False,
    complete_full_electric_Gram_includes_reference_loops=True,
    actual_bare_energy_bounds_uniform_in_total_volume_for_fixed_support=True,
    actual_normalized_state_energy_or_variance_floor_proved=False,
    actual_energy_tail_bound_requires_complete_native_cutoff=True,
    actual_projection_can_be_replaced_by_Haar_projection_without_comparison=False,
    response_Gram_contains_full_outputs=True,
    energy_constrained_bound_is_unrestricted_operator_norm=False,
    actual_weak_coupling_response_matrix_bound_supplied=False,
    density_ratio_bound_uniform_in_coupling_or_growing_support=False,
    comparison_fixtures_are_interacting_vacua=False,
    positive_four_dimensional_continuum_mass_gap=False,
    continuum_state_survival_proved=False,
    analytic_proofs_formally_machine_verified=False)


def invariants(a,b,x):
    a,b,x = y.unit(a),y.unit(b),y.unit(x)
    return (a[0],b[0],x[0],y.dot(a[1:],b[1:]),y.dot(a[1:],x[1:]),
            y.dot(b[1:],x[1:]),y.dot(y.cross(a[1:],b[1:]),x[1:]))


def determinant3(a):
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))


def trace_control(a,b,x):
    a,b,x = y.unit(a),y.unit(b),y.unit(x)
    a0,b0,x0,u,d1,d2,d3 = invariants(a,b,x)
    ab,ax,bx = y.qmul(a,b)[0],y.qmul(a,x)[0],y.qmul(b,x)[0]
    abx,bax = y.qmul(y.qmul(a,b),x)[0],y.qmul(y.qmul(b,a),x)[0]
    assert (u,d1,d2)==(a0*b0-ab,a0*x0-ax,b0*x0-bx)
    assert d3==(bax-abx)/2==a0*b0*x0-x0*u-a0*d2-b0*d1-abx
    gram=[[1-a0*a0,u,d1],[u,1-b0*b0,d2],[d1,d2,1-x0*x0]]
    assert determinant3(gram)==d3*d3 and p.psd(gram)
    return dict(readings=(a0,b0,x0,u,d1,d2,d3), ordered_traces=(abx,bax),
                imaginary_Gram=gram, determinant=d3*d3)


def relational_derivative(a,b,x,da,db,dx):
    av,bv,xv=a[1:],b[1:],x[1:]
    dav,dbv,dxv=da[1:],db[1:],dx[1:]
    return (y.dot(dav,xv)+y.dot(av,dxv), y.dot(dbv,xv)+y.dot(bv,dxv),
            y.dot(y.cross(dav,bv),xv)+y.dot(y.cross(av,dbv),xv)+y.dot(y.cross(av,bv),dxv))


def energy_control(a,b,x):
    a,b,x = y.unit(a),y.unit(b),y.unit(x)
    _,_,_,u,d1,d2,d3=invariants(a,b,x)
    av,bv,xv=a[1:],b[1:],x[1:]
    aa,bb,xx=y.dot(av,av),y.dot(bv,bv),y.dot(xv,xv)
    triple=sum(y.dot(v,v) for v in (y.cross(bv,xv),y.cross(xv,av),y.cross(av,bv)))-3*d3*d3
    expected=[[aa+xx-2*d1*d1,u-d1*d2,-2*d1*d3],
              [u-d1*d2,bb+xx-2*d2*d2,-2*d2*d3],
              [-2*d1*d3,-2*d2*d3,triple]]
    actual=[[F(0) for _ in range(3)] for _ in range(3)]
    for loop,q in enumerate((a,b,x)):
        for axis in (1,2,3):
            directions=[(F(0),)*4 for _ in range(3)]
            directions[loop]=y.native_direction(q,axis)
            values=relational_derivative(a,b,x,*directions)
            for i in range(3):
                for j in range(3):
                    actual[i][j]+=4*values[i]*values[j]
    assert actual==expected and p.psd(actual)
    assert actual[0][0]<=aa+xx and actual[1][1]<=bb+xx and actual[2][2]<=aa+bb+xx
    k=[list(av),list(bv),list(y.cross(av,bv))]
    delta=y.dot(k[2],k[2])
    gram=p.multiply(k,p.transpose(k))
    assert p.psd(p.add(gram,[[delta*v/2 for v in row] for row in p.identity(3)],-1))
    return dict(full_electric_Gram=actual, target_readout_Gram=gram,
                reference_cross_squared=delta,
                inverse_squared_norm_upper=2/delta if delta else None,
                finite_configuration_not_actual_expectation=True)


def mirror_control(parameter):
    z=p.positive(parameter)
    if z>=1:
        raise ValueError('Use a rational parameter in (0,1)')
    c,t=(1-z*z)/(1+z*z),2*z/(1+z*z)
    a,b=(c,t,0,0),(c,0,t,0)
    plus,minus=invariants(a,b,(c,0,0,t)),invariants(a,b,(c,0,0,-t))
    assert plus[:6]==minus[:6] and plus[6]==-minus[6]==t**3
    return dict(parameter=z, imaginary_radius=t, shared_pair_readings=plus[:6],
                signed_triple=(plus[6],minus[6]), orbit_distance=2*t,
                reading_distance=2*t**3, inverse_Lipschitz_ratio=1/t**2)


def loop_word(root,first,second,side):
    current=root
    word=[]
    for axis,sign in (first,second,(first[0],-first[1]),(second[0],-second[1])):
        end=y.shifted(current,axis,sign,side)
        word.append(((current,axis) if sign==1 else (end,axis),sign))
        current=end
    assert current==root
    return tuple(word)


def bouquet(side):
    graph=y.lattice(side)
    root=(0,0,0)
    words=(loop_word(root,(0,-1),(2,1),side),loop_word(root,(1,-1),(2,-1),side),
           loop_word(root,(0,1),(1,1),side))  # A, B, X
    supports=[set(e for e,_ in word) for word in words]
    assert all(len(s)==4 for s in supports)
    assert all(not a&b for a,b in combinations(supports,2))
    support=set.union(*supports)
    faces=[set(e for e,_ in loop_word(v,(i,1),(j,1),side))
           for v in graph['vertices'] for i,j in combinations(range(3),2)]
    counts={k:sum(len(f&support)==k for f in faces) for k in range(1,5)}
    assert counts=={1:18,2:9,3:0,4:3}
    return dict(graph=graph,words=words,support=support,incidence_counts=counts)


def density_budget(links,faces,theta):
    p.size(links)
    p.size(faces)
    theta=p.nonnegative(theta)
    if not theta:
        return dict(theta=theta,exact_ratio=F(1),actual_marginal=True)
    omega=2*faces*theta
    root=p.sqrt_interval(3*links*omega)
    exponent=4*F(22,7)*root['upper']
    return dict(theta=theta,links=links,incident_faces=faces,omega=omega,
                fixed_time=dict(prefactor=F(21,19)**(2*links),exp_exponent=18*omega),
                optimized=dict(prefactor=64**links,exp_exponent_upper=exponent,
                               sqrt_argument=3*links*omega),
                total_volume_independent_for_fixed_support=True,coupling_uniform=False,
                actual_marginal=True)


def solve(matrix,right):
    n=len(matrix)
    if not n or len(right)!=n or any(len(row)!=n for row in matrix):
        raise ValueError('Invalid square linear system')
    rows=[[p.exact(x) for x in row]+[p.exact(b)] for row,b in zip(matrix,right)]
    for col in range(n):
        pivot=next((i for i in range(col,n) if rows[i][col]),None)
        if pivot is None:
            raise ValueError('Singular linear system')
        rows[col],rows[pivot]=rows[pivot],rows[col]
        value=rows[col][col]
        rows[col]=[x/value for x in rows[col]]
        for i in range(n):
            if i!=col:
                factor=rows[i][col]
                rows[i]=[a-factor*b for a,b in zip(rows[i],rows[col])]
    return [row[-1] for row in rows]


def derivative(values):
    return [i*values[i] for i in range(1,len(values))] or [F(0)]


def weighted_projection_control(eta=F(3,4),cutoff=3):
    eta=p.exact(eta)
    p.size(cutoff)
    if abs(eta)>=1:
        raise ValueError('The comparison Haar density must be strictly positive')
    weight=[1,eta]
    integral=lambda f: poly.haar_integral(poly.poly_multiply(weight,f))
    inner=lambda f,g: integral(poly.poly_multiply(f,g))
    basis=[poly.character(n) for n in range(cutoff+1)]
    f=poly.poly_add(poly.character(1),poly.character(cutoff+1),F(1,20))
    f=poly.poly_add(f,poly.character(cutoff+3),F(1,100))
    f=poly.poly_add(f,[integral(f)],-1)
    gram=[[inner(a,b) for b in basis] for a in basis]
    right=[inner(b,f) for b in basis]
    coefficients=solve(gram,right)
    projected=[F(0)]
    haar_projected=[F(0)]
    for coefficient,b in zip(coefficients,basis):
        projected=poly.poly_add(projected,b,coefficient)
        haar_projected=poly.poly_add(haar_projected,b,poly.haar_integral(poly.poly_multiply(b,f)))
    remainder=poly.poly_add(f,projected,-1)
    haar_remainder=poly.poly_add(f,haar_projected,-1)
    energy=integral(poly.poly_multiply([1,0,-1],poly.poly_multiply(derivative(f),derivative(f))))
    tail=inner(remainder,remainder)
    ratio=(1+abs(eta))/(1-abs(eta))
    threshold=cutoff*(cutoff+2)
    assert tail<=ratio*energy/threshold and tail<=inner(haar_remainder,haar_remainder)
    assert not integral(projected) and all(not inner(remainder,b) for b in basis)
    return dict(kind='weighted_Haar_polynomial_comparison_not_actual_vacuum',eta=eta,
                cutoff_twice_spin=cutoff,cutoff_energy=threshold,density_ratio=ratio,
                function=f,actual_metric_projection=projected,Haar_projection=haar_projected,
                actual_energy=energy,function_norm_squared=inner(f,f),tail_squared=tail,
                tail_upper=ratio*energy/threshold,projection_coefficients=coefficients)


def positive_definite(a):
    if not a or any(len(row)!=len(a) for row in a) or a!=p.transpose(a):
        return False
    a=[[p.exact(x) for x in row] for row in a]
    while a:
        pivot=a[0][0]
        if pivot<=0:
            return False
        a=[[a[i][j]-a[i][0]*a[0][j]/pivot for j in range(1,len(a))] for i in range(1,len(a))]
    return True


def response_budget(covariance,output_gram,finite_upper,density_ratio,energy,cutoff):
    c=[[p.exact(x) for x in row] for row in covariance]
    r=[[p.exact(x) for x in row] for row in output_gram]
    if not positive_definite(c) or len(r)!=len(c) or any(len(row)!=len(c) for row in r):
        raise ValueError('Independent centered basis and matching matrices required')
    finite=p.nonnegative(finite_upper)
    ratio=p.positive(density_ratio)
    energy,cutoff=p.nonnegative(energy),p.positive(cutoff)
    if finite>1 or ratio<1 or not p.psd(r) or not p.psd(p.add(c,r,-1)):
        raise ValueError('Invalid contraction Gram or density ratio')
    if not p.psd(p.add([[finite*finite*x for x in row] for row in c],r,-1)):
        raise ValueError('Proposed finite response upper bound fails')
    tail=p.sqrt_interval(ratio*energy/cutoff)
    return dict(finite_response_upper=finite,energy=energy,cutoff=cutoff,density_ratio=ratio,
                omitted_input_norm_upper=min(F(1),tail['upper']),
                energy_constrained_response_upper=min(F(1),finite+tail['upper']),
                complete_native_cutoff_required=True,full_output_Gram_required=True,
                actual_response_hypothesis_discharged=False,
                unrestricted_operator_norm_upper=False)


def three_loop_rank_upper(cutoff):
    p.size(cutoff,0)
    maximum=isqrt(cutoff+1)-1
    modes=[(n*(n+2),(n+1)**2) for n in range(maximum+1)]
    return sum(a[1]*b[1]*c[1] for a,b,c in product(modes,repeat=3) if a[0]+b[0]+c[0]<=cutoff)


def controls():
    triples=((y.QA,y.QB,y.QX),(y.ONE,y.QA,y.QB),(y.QA,y.QA,y.QX),
             ((0,1,0,0),(0,0,1,0),(0,0,0,1)))
    return dict(trace_controls=[trace_control(*row) for row in triples],
                full_energy_controls=[energy_control(*row) for row in triples],
                orientation_blindness=[mirror_control(F(1,n)) for n in (2,4,16,64)],
                local_geometry=[dict(side=n,links=len(bouquet(n)['support']),
                                     incident_face_multiplicities=bouquet(n)['incidence_counts']) for n in (4,5,6)],
                actual_density_budgets=[density_budget(12,30,theta) for theta in (F(0),F(1,16000),F(1),F(25))],
                weighted_polynomial_controls=[weighted_projection_control(eta,n) for eta in (F(0),F(3,4),F(-2,3)) for n in (2,4)],
                finite_response_controls=[response_budget([[1,0],[0,2]],[[F(1,16),0],[0,F(1,2)]],F(1,2),2,1,cutoff)
                                          for cutoff in (16,128,512)],
                unprojected_three_loop_rank_upper={str(c):three_loop_rank_upper(c) for c in (3,6,9,12,24)},
                actual_energy_upper_coefficients_over_sqrt_theta=(F(9),F(9),F(27,2)))


def source_checks():
    manifest=json.loads(SOURCES.read_text())
    if manifest['stage']!='YM86' or manifest['base_commit']!=BASE_COMMIT or manifest['base_tree']!=BASE_TREE:
        raise ValueError('YM86 source base mismatch')
    if len(manifest['upstream_sha256'])!=375 or len(manifest['local_inputs'])!=4:
        raise ValueError('YM86 source manifest shape mismatch')
    if y.PIN.read_text().strip()!=PARENT_HASH or p.canonical_sha(json.loads(y.RESULT.read_text()))!=PARENT_HASH:
        raise ValueError('YM85 result/pin mismatch')
    checked={}
    for path,expected in manifest['upstream_sha256'].items():
        actual=sha256((ROOT/path).read_bytes()).hexdigest()
        if actual!=expected:
            raise ValueError('Changed upstream input: '+path)
        checked[path]=actual
    for path in manifest['local_inputs']:
        checked[path]=sha256((ROOT/path).read_bytes()).hexdigest()
    checked[str(SOURCES.relative_to(ROOT))]=sha256(SOURCES.read_bytes()).hexdigest()
    return checked


def run():
    if sys.version_info[:2]!=(3,12) or sys.flags.optimize:
        raise RuntimeError('Use Python 3.12 without optimization')
    return p.encode(dict(certificate_type='YM86_FINITE_ENERGY_RESPONSE',verdict='PASS',
                         runtime_policy='Python 3.12 only',model=MODEL,
                         input_sha256=source_checks(),controls=controls(),
                         claim_status='SMOOTH_LOCAL_MEASUREMENTS__ACTUAL_ENERGY_TAIL__CONTINUUM_OPEN'))


def check(cert,result=RESULT,pin=PIN):
    if cert.get('model')!=MODEL or cert.get('certificate_type')!='YM86_FINITE_ENERGY_RESPONSE':
        raise ValueError('YM86 scope mismatch')
    digest=p.canonical_sha(cert)
    if json.loads(result.read_text())!=cert or pin.read_text().strip()!=digest:
        raise ValueError('YM86 fresh result/pin mismatch')
    return digest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument('--write',action='store_true')
    group.add_argument('--check',action='store_true')
    args=parser.parse_args()
    cert=run()
    digest=p.canonical_sha(cert)
    if args.write:
        RESULT.write_text(json.dumps(cert,sort_keys=True,indent=2)+'\n')
        PIN.write_text(digest+'\n')
    else:
        check(cert)
    print('YM86 PASS',digest)
    print('Smooth local measurements and actual energy tails: proved; full 4D continuum gap: open')


if __name__=='__main__':
    main()
