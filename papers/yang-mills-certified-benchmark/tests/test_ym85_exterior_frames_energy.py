import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb, cos, sin, sqrt
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'certificates'))
import ym85_exterior_frames_energy as y


def decimal(x):
    return Decimal(x.numerator)/Decimal(x.denominator)


def complex_add(a, b):
    return (a[0]+b[0], a[1]+b[1])


def complex_mul(a, b):
    return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])


def matrix(q):
    a,b,c,d = q
    return [[(a,b), (c,d)], [(-c,d), (a,-b)]]


def matrix_product(a,b):
    return [[complex_add(complex_mul(a[i][0], b[0][j]), complex_mul(a[i][1], b[1][j]))
             for j in range(2)] for i in range(2)]


def rank(matrix):
    rows = [list(row) for row in matrix]
    r = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(r, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        value = rows[r][col]
        rows[r] = [x/value for x in rows[r]]
        for i in range(r+1, len(rows)):
            factor = rows[i][col]
            rows[i] = [a-factor*b for a,b in zip(rows[i],rows[r])]
        r += 1
        if r == len(rows):
            break
    return r


def float_frame(a,b):
    av,bv = a[1:],b[1:]
    r = sqrt(sum(v*v for v in av))
    n = [v/r for v in av]
    parallel = sum(x*z for x,z in zip(n,bv))
    perp = [v-parallel*x for v,x in zip(bv,n)]
    s = sqrt(sum(v*v for v in perp))
    second = [v/s for v in perp]
    third = [n[1]*second[2]-n[2]*second[1], n[2]*second[0]-n[0]*second[2], n[0]*second[1]-n[1]*second[0]]
    return [n,second,third]


class FrameGeometryTests(unittest.TestCase):
    def test_quaternion_product_against_independent_complex_matrices(self):
        for a in (y.ONE,y.QA,y.QB,y.QX):
            for b in (y.ONE,y.QA,y.QB,y.QX):
                self.assertEqual(matrix(y.qmul(a,b)), matrix_product(matrix(a),matrix(b)))
                self.assertEqual(y.qmul(a,y.qdagger(a)), y.ONE)
                self.assertEqual(y.dot(y.qmul(a,b),y.qmul(a,b)), 1)

    def test_frame_is_oriented_orthonormal_and_gauge_covariant(self):
        original = y.body_frame(y.QA,y.QB)
        expected = y.canonical(y.QX,original)
        for g in (y.ONE,y.QA,y.QB,y.QX):
            rotate = lambda q: y.qmul(y.qmul(g,q),y.qdagger(g))
            frame = y.body_frame(rotate(y.QA),rotate(y.QB))
            for i in range(3):
                for j in range(3):
                    self.assertEqual(y.dot(frame['axes'][i],frame['axes'][j]), F(i==j))
            self.assertEqual(y.cross(frame['axes'][0],frame['axes'][1]),frame['axes'][2])
            self.assertEqual(y.canonical(rotate(y.QX),frame),expected)

    def test_two_nonparallel_reference_vectors_leave_only_central_stabilizer(self):
        def cross_matrix(v):
            a,b,c = v
            return [[0,-c,b],[c,0,-a],[-b,a,0]]
        self.assertEqual(rank(cross_matrix(y.QA[1:])+cross_matrix(y.QB[1:])),3)
        self.assertEqual(rank(cross_matrix(y.QA[1:])+cross_matrix(y.QA[1:])),2)

    def test_singular_frames_and_inexact_fixture_lengths_are_rejected(self):
        for a,b in ((y.ONE,y.QB),(y.QA,y.QA),(y.QA,y.ONE)):
            with self.assertRaises(ValueError):
                y.body_frame(a,b)
        with self.assertRaises(ValueError):
            y.rational_sqrt(2)
        with self.assertRaises(ValueError):
            y.unit((1,1,0,0))
        with self.assertRaises(TypeError):
            y.unit((1.0,0,0,0))

    def test_every_cubic_edge_has_four_disjoint_correct_detours(self):
        for side in (4,5):
            graph = y.lattice(side)
            for edge in graph['edges']:
                paths = y.detours(edge,side)
                self.assertEqual(len(paths),4)
                used = set()
                for path in paths:
                    current = graph['ends'][edge][0]
                    for alternate,sign in path:
                        left,right = graph['ends'][alternate]
                        if sign < 0:
                            left,right = right,left
                        self.assertEqual(current,left)
                        current=right
                        self.assertNotIn(alternate,used)
                        self.assertNotEqual(alternate,edge)
                        used.add(alternate)
                    self.assertEqual(current,graph['ends'][edge][1])

    def test_four_deleted_links_leave_a_complete_exterior_in_adversarial_cases(self):
        graph = y.lattice(4)
        edge = graph['edges'][0]
        alternate_edges = [e for path in y.detours(edge,4) for e,_ in path]
        for obstruct in combinations(alternate_edges,3):
            removed = {edge,*obstruct}
            tree = y.exterior_tree(graph,removed)
            self.assertEqual(len(tree['tree']),63)
            self.assertEqual(len(tree['chords']),125)
            self.assertFalse(tree['tree'] & removed)

    def test_large_deletion_can_disconnect_exterior(self):
        graph = y.lattice(4)
        root = graph['vertices'][0]
        star = {e for e,ends in graph['ends'].items() if root in ends}
        self.assertEqual(len(star),6)
        with self.assertRaisesRegex(ValueError,'disconnected'):
            y.exterior_tree(graph,star)

    def test_local_vertex_gauges_preserve_all_canonical_link_variables(self):
        for side in (4,5):
            out = y.lattice_frame_control(side)
            self.assertTrue(out['all_canonical_links_gauge_invariant'])
            self.assertEqual(out['exterior_cycle_rank'],2*side**3-3)

    def test_tree_transports_obey_endpoint_gauge_law(self):
        graph = y.lattice(4)
        tree = y.exterior_tree(graph,graph['edges'][:2])
        palette = (y.QA,y.QB,y.QX,y.ONE)
        field = {e:palette[i%4] for i,e in enumerate(graph['edges'])}
        gauges = {v:palette[(2*v[0]+v[1]+3*v[2]+1)%4] for v in graph['vertices']}
        transformed = y.gauge_transform(graph,field,gauges)
        before,after = y.tree_transports(tree,field),y.tree_transports(tree,transformed)
        for v in graph['vertices']:
            self.assertEqual(after[v], y.qmul(y.qmul(gauges[tree['root']],before[v]),y.qdagger(gauges[v])))


class PhysicalResponseTests(unittest.TestCase):
    def test_finite_reference_is_positive_normalized_and_gauge_invariant(self):
        data = y.finite_reference_model()
        states,weights = data['states'],data['weights']
        lookup = dict(zip(states,weights))
        self.assertEqual(sum(weights),1)
        self.assertGreater(min(weights),0)
        for x,z,r,t in states:
            self.assertEqual(lookup[x,z,r,t],lookup[-x,-z,-r,t])
        g=data['gauge']
        self.assertEqual(y.p.multiply(g,g),g)
        self.assertEqual(rank(g),8)
        self.assertEqual(y.p.multiply(data['response'],g),y.p.multiply(g,data['response']))

    def test_conditional_binary_tables_give_the_claimed_maximal_correlation(self):
        data=y.finite_reference_model()
        for r in (-1,1):
            for z,t in data['ts'].items():
                table={(x,b):sum(w for state,w in zip(data['states'],data['weights'])
                                  if state==(x,b,r,z)) for x in (-1,1) for b in (-1,1)}
                total=sum(table.values())
                table={k:v/total for k,v in table.items()}
                self.assertEqual(sum(x*b*w for (x,b),w in table.items()),t)
                self.assertEqual(sum(x*w for (x,b),w in table.items()),0)
                self.assertEqual(sum(b*w for (x,b),w in table.items()),0)

    def test_physical_witness_retains_rare_conditional_operator_norm(self):
        for probability in (F(1,2),F(1,100),F(1,1000000)):
            data=y.finite_reference_model(rare_probability=probability)
            witness=[F(x*r) if z==1 else F(0) for x,b,r,z in data['states']]
            image=y.previous.apply(data['response'],witness)
            norm=lambda v:sum(w*x*x for w,x in zip(data['weights'],v))
            self.assertEqual(norm(witness),probability)
            self.assertEqual(norm(image)/norm(witness),F(49,64))
            self.assertEqual(y.previous.apply(data['gauge'],witness),witness)
            self.assertEqual(data['physical_norm'],data['full_norm'])

    def test_finite_physical_partial_floor_matches_conditional_flip_floor(self):
        data=y.finite_reference_model()
        states,weights=data['states'],data['weights']
        witness=[F(x*r) if z==-1 else F(0) for x,b,r,z in states]
        index={s:i for i,s in enumerate(states)}
        energy=sum(w*(witness[index[-x,b,r,z]]-value)**2/2
                   for (x,b,r,z),w,value in zip(states,weights,witness))
        px=y.p.conditional_projection(states,weights,[0])
        residual=[x-z for x,z in zip(witness,y.previous.apply(px,witness))]
        variance=sum(w*x*x for w,x in zip(weights,residual))
        conditional_floor=min(2/(1-t*t) for t in data['ts'].values())
        self.assertEqual(energy/variance,conditional_floor)
        self.assertEqual(y.previous.apply(data['gauge'],witness),witness)

    def test_finite_reference_rejects_missing_full_support(self):
        for kwargs in (dict(t_plus=1),dict(t_minus=-1),dict(rare_probability=0),dict(rare_probability=1)):
            with self.assertRaises(ValueError):
                y.finite_reference_model(**kwargs)


class EnergyDomainTests(unittest.TestCase):
    def test_four_link_product_derivatives_give_full_plaquette_gram(self):
        quaternions=(y.QA,y.QB,y.QX,y.unit((F(5,13),0,F(12,13),0)))
        frame=y.body_frame(y.QA,y.QB)
        signs=(1,1,-1,-1)
        factors=[q if sign==1 else y.qdagger(q) for q,sign in zip(quaternions,signs)]
        def product(values):
            out=y.ONE
            for q in values:
                out=y.qmul(out,q)
            return out
        holonomy=product(factors)
        marker=y.canonical(holonomy,frame)[1:]
        gram=[[F(0) for _ in range(3)] for _ in range(3)]
        for e in range(4):
            for axis in (1,2,3):
                dq=y.native_direction(quaternions[e],axis)
                if signs[e]<0:
                    dq=y.scale(y.qmul(y.qmul(y.qdagger(quaternions[e]),dq),y.qdagger(quaternions[e])),-1)
                modified=list(factors)
                modified[e]=dq
                derivative=product(modified)
                values=[y.dot(a,derivative[1:]) for a in frame['axes']]
                for i in range(3):
                    for j in range(3):
                        gram[i][j]+=values[i]*values[j]
        self.assertEqual(gram,[[F(i==j)-marker[i]*marker[j] for j in range(3)] for i in range(3)])

    def test_frame_differentials_match_independent_geodesic_differences(self):
        step=1e-6
        for changing in ('a','b'):
            for axis in (1,2,3):
                da=y.native_direction(y.QA,axis) if changing=='a' else (F(0),)*4
                db=y.native_direction(y.QB,axis) if changing=='b' else (F(0),)*4
                expected=y.frame_differential(y.QA,y.QB,da,db)
                def displaced(q,dq,amount):
                    if not any(dq):
                        return list(map(float,q))
                    return [float(x)*cos(amount/2)+2*float(v)*sin(amount/2) for x,v in zip(q,dq)]
                plus=float_frame(displaced(y.QA,da,step),displaced(y.QB,db,step))
                minus=float_frame(displaced(y.QA,da,-step),displaced(y.QB,db,-step))
                for i in range(3):
                    for j in range(3):
                        self.assertAlmostEqual((plus[i][j]-minus[i][j])/(2*step),float(expected[i][j]),places=8)

    def test_exterior_chord_cost_matches_inverse_perpendicular_radius(self):
        previous_cost=F(0)
        for denominator in (2,4,8,16,32,64):
            t=F(1,denominator)
            b=(0,(1-t*t)/(1+t*t),2*t/(1+t*t),0)
            out=y.pointwise_energy(y.QA,b,y.QX)
            expected=1/(16*out['frame_s']**2)
            self.assertEqual(out['exterior_chord_energy'],[0,expected,expected])
            self.assertGreater(expected,previous_cost)
            previous_cost=expected

    def test_haar_perpendicular_squared_radius_has_uniform_moments(self):
        for n in range(20):
            moment=sum((comb(n,k)*y.sphere_moment((0,0,k,n-k)) for k in range(n+1)),F(0))
            self.assertEqual(moment,F(1,n+1))

    def test_first_axis_inverse_radius_integral_has_exact_convergent_series(self):
        for n in (1,2,4,16,64):
            partial=sum((y.sphere_moment((k,0,0,0)) for k in range(n)),F(0))
            remainder=F(2*comb(2*n,n),4**n)
            self.assertEqual(partial+remainder,2)
            self.assertGreater(remainder,0)

    def test_second_axis_positive_inverse_radius_series_diverges(self):
        for exponent in range(1,9):
            n=2**exponent
            harmonic=sum((F(1,k) for k in range(1,n+1)),F(0))
            self.assertGreaterEqual(harmonic,1+F(exponent,2))
            # Every term integrates to 1/(k+1) under the uniform u law.
            for k in range(min(n,20)):
                integral=sum((F((-1)**j*comb(k,j),j+1) for j in range(k+1)),F(0))
                self.assertEqual(integral,F(1,k+1))

    def test_logarithmic_cutoff_energy_interval_contains_independent_log(self):
        with localcontext() as context:
            context.prec=80
            for level in (1,4,16,64):
                out=y.singularity_control(level)
                exact=Decimal(level)*Decimal(2).ln()/8
                bounds=out['Haar_exterior_energy_interval']
                self.assertLess(decimal(bounds['lower']),exact)
                self.assertGreater(decimal(bounds['upper']),exact)
                self.assertEqual(out['epsilon'],F(1,2**level))
                for key in ('lower','upper'):
                    self.assertEqual(out['actual_cutoff_energy_coefficient_of_minimum_density'][key],bounds[key]/2)

    def test_actual_stiffness_thresholds_and_physical_clock(self):
        for s in (F(5),F(11,2),F(9),F(100)):
            out=y.stiffness_control(s)
            self.assertEqual(out['specified_marker_partial_quotient_lower'],2*s/9-1)
            if s>=9:
                self.assertGreaterEqual(out['specified_marker_partial_quotient_lower'],s/9)
                for a in (F(1,10),F(1,1000)):
                    kappa=1/(a*s)
                    self.assertGreaterEqual(kappa*out['specified_marker_partial_quotient_lower'],1/(9*a))
        with self.assertRaises(ValueError):
            y.stiffness_control(4)

    def test_covariance_and_partial_energy_matrix_bound_for_a_separate_fixture(self):
        # This finite law only checks matrix algebra, not actual vacuum moments.
        small=F(1,10)
        scalar=(1-small*small)/(1+small*small)
        radius=2*small/(1+small*small)
        markers=[tuple(sign*radius*F(i==axis) for i in range(3)) for axis in range(3) for sign in (-1,1)]
        m=[[sum(f[i]*f[j] for f in markers)/6 for j in range(3)] for i in range(3)]
        g=[[F(i==j)-m[i][j] for j in range(3)] for i in range(3)]
        delta=1-scalar
        lower=(1-2*delta)/(2*delta)
        self.assertTrue(y.p.psd(y.p.add(g,m,-lower)))
        self.assertEqual(sum(m[i][i] for i in range(3)),1-scalar*scalar)
        self.assertEqual(sum(m[i][i] for i in range(3)),radius*radius)


class EvidenceTests(unittest.TestCase):
    def test_reproduction_and_scope_do_not_promote_partial_results_to_gap(self):
        cert=y.run()
        self.assertEqual(len(cert['input_sha256']),373)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        for name in ('partial_marker_stiffness_is_full_physical_gap',
                     'second_and_third_normalized_frame_markers_have_finite_full_energy',
                     'full_lattice_weak_coupling_gap_proved',
                     'continuum_state_survival_proved',
                     'positive_four_dimensional_continuum_mass_gap'):
            self.assertFalse(cert['model'][name])

    def test_changed_result_pin_and_scope_are_rejected(self):
        cert=y.run()
        with tempfile.TemporaryDirectory() as directory:
            result,pin=Path(directory)/'result',Path(directory)/'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            pin.write_text(y.p.canonical_sha(cert))
            changed=copy.deepcopy(cert)
            changed['controls']['first_axis_Haar_inverse_squared_radius']='3'
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed=copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap']=True
            with self.assertRaises(ValueError):
                y.check(changed,result,pin)

    def test_upstream_source_drift_is_rejected(self):
        manifest=json.loads(y.SOURCES.read_text())
        path=next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][path]='0'*64
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'manifest'
            source.write_text(json.dumps(manifest))
            with patch.object(y,'SOURCES',source):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    y.source_checks()

    def test_cli_is_read_only_and_rejects_optimized_runtime(self):
        inputs=(y.RESULT,y.PIN,y.SOURCES)
        before=[path.read_bytes() for path in inputs]
        checked=subprocess.run([sys.executable,'-B',str(y.HERE/'ym85_exterior_frames_energy.py'),'--check'],capture_output=True,text=True)
        self.assertEqual(checked.returncode,0,checked.stderr)
        self.assertIn('YM85 PASS',checked.stdout)
        self.assertEqual(before,[path.read_bytes() for path in inputs])
        optimized=subprocess.run([sys.executable,'-B','-O',str(y.HERE/'ym85_exterior_frames_energy.py'),'--check'],capture_output=True,text=True)
        self.assertNotEqual(optimized.returncode,0)
        self.assertIn('without optimization',optimized.stderr)


if __name__=='__main__':
    unittest.main()
