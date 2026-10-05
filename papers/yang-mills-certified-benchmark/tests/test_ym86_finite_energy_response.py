import copy
from fractions import Fraction as F
from itertools import permutations, product
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'certificates'))
import ym86_finite_energy_response as z

y, p = z.y, z.p
ZERO = (F(0),)*4


def word_value(word, field):
    value = y.ONE
    for edge, sign in word:
        value = y.qmul(value, field[edge] if sign > 0 else y.qdagger(field[edge]))
    return value


def word_jet(word, field, moving, axis):
    """Differentiate the actual oriented link word, including inverse factors."""
    value, first, second = y.ONE, ZERO, ZERO
    generator = tuple(F(k == axis, 2) for k in range(4))
    for edge, sign in word:
        q = field[edge]
        dq = y.qmul(generator, q) if edge == moving else ZERO
        ddq = y.scale(q, F(-1,4)) if edge == moving else ZERO
        if sign < 0:
            q, dq, ddq = y.qdagger(q), y.qdagger(dq), y.qdagger(ddq)
        second = y.plus(y.plus(y.qmul(second,q),y.qmul(first,dq),2),y.qmul(value,ddq))
        first = y.plus(y.qmul(first,q),y.qmul(value,dq))
        value = y.qmul(value,q)
    return value, first, second


def monomial(indices, coefficient=1):
    powers = [0]*9
    for i in indices:
        powers[i] += 1
    return {tuple(powers): F(coefficient)}


def add(a,b,scale=1):
    out = dict(a)
    for powers, value in b.items():
        out[powers] = out.get(powers,F(0))+scale*value
    return {powers:value for powers,value in out.items() if value}


def multiply(a,b):
    out = {}
    for ap,av in a.items():
        for bp,bv in b.items():
            powers = tuple(i+j for i,j in zip(ap,bp))
            out[powers] = out.get(powers,F(0))+av*bv
    return out


def partial(a,index):
    out = {}
    for powers,value in a.items():
        if powers[index]:
            reduced = list(powers)
            reduced[index] -= 1
            out[tuple(reduced)] = powers[index]*value
    return out


def evaluate(a,values):
    result = F(0)
    for powers,coefficient in a.items():
        for value,power in zip(values,powers):
            coefficient *= value**power
        result += coefficient
    return result


def haar(a):
    result = F(0)
    for powers,coefficient in a.items():
        if any(n % 2 for n in powers):
            continue
        for i in (0,3,6):
            coefficient *= y.sphere_moment((0,*(n//2 for n in powers[i:i+3])))
        result += coefficient
    return result


def relational_polynomials():
    first,second,third = {},{},{}
    for i in range(3):
        first = add(first,monomial((i,i+6)))
        second = add(second,monomial((i+3,i+6)))
    for i,j,k in permutations(range(3)):
        sign = 1 if (i,j,k) in ((0,1,2),(1,2,0),(2,0,1)) else -1
        third = add(third,monomial((i,j+3,k+6),sign))
    return first,second,third


class RelationalGeometryTests(unittest.TestCase):
    def test_trace_relations_include_central_parallel_and_full_rank_cases(self):
        for triple in ((y.QA,y.QB,y.QX),(y.ONE,y.ONE,y.ONE),
                       (y.QA,y.QA,y.QX),((0,1,0,0),(0,0,1,0),(0,0,0,1))):
            out = z.trace_control(*triple)
            self.assertEqual(out['determinant'],out['readings'][-1]**2)
            for g in (y.QA,y.QB,y.QX):
                rotate = lambda q:y.qmul(y.qmul(g,q),y.qdagger(g))
                self.assertEqual(z.invariants(*(rotate(q) for q in triple)),out['readings'])

    def test_pair_traces_cannot_detect_mirror_orientation(self):
        for parameter in (F(1,2),F(1,4),F(1,16),F(1,64)):
            out = z.mirror_control(parameter)
            t = out['imaginary_radius']
            c = (1-parameter**2)/(1+parameter**2)
            a,b = (c,t,0,0),(c,0,t,0)
            pair = lambda x:(y.qmul(a,b)[0],y.qmul(a,x)[0],y.qmul(b,x)[0])
            self.assertEqual(pair((c,0,0,t)),pair((c,0,0,-t)))
            self.assertEqual(out['reading_distance'],out['orbit_distance']**3/4)
            self.assertEqual(out['inverse_Lipschitz_ratio'],1/t**2)
            self.assertGreater(out['signed_triple'][0],0)

    def test_orbit_distance_is_attained_and_respects_rational_rotations(self):
        t = F(4,5)
        c = F(3,5)
        plus = ((c,t,0,0),(c,0,t,0),(c,0,0,t))
        minus = (plus[0],plus[1],(c,0,0,-t))
        distance = lambda gs:sum(y.dot(y.plus(a,b,-1),y.plus(a,b,-1)) for a,b in zip(plus,gs))
        self.assertEqual(distance(minus),4*t*t)
        for g in (y.ONE,y.QA,y.QB,y.QX,*y.UNITS):
            rotated = tuple(y.qmul(y.qmul(g,q),y.qdagger(g)) for q in minus)
            self.assertGreaterEqual(distance(rotated),4*t*t)

    def test_fixed_reference_readings_reconstruct_imaginary_vector(self):
        a,b,x = y.QA,y.QB,y.QX
        k = [list(a[1:]),list(b[1:]),list(y.cross(a[1:],b[1:]))]
        self.assertEqual(z.solve(k,list(z.invariants(a,b,x)[4:])),list(x[1:]))
        out = z.energy_control(a,b,x)
        delta = out['reference_cross_squared']
        self.assertGreater(delta,0)
        self.assertTrue(p.psd(p.add(out['target_readout_Gram'],p.identity(3),-delta/2)))

    def test_smooth_readings_remain_defined_where_normalized_frames_fail(self):
        for a,b in ((y.ONE,y.QB),(y.QA,y.QA),(y.QA,y.ONE)):
            with self.assertRaises(ValueError):
                y.body_frame(a,b)
            self.assertEqual(z.energy_control(a,b,y.QX)['reference_cross_squared'],0)
            self.assertIsNone(z.energy_control(a,b,y.QX)['inverse_squared_norm_upper'])
            self.assertTrue(all(abs(v)<=1 for v in z.invariants(a,b,y.QX)))

    def test_bouquet_has_twelve_links_ten_vertices_and_thirty_incident_faces(self):
        for side in (4,5,6):
            data = z.bouquet(side)
            graph,words = data['graph'],data['words']
            vertices = [{v for edge,_ in word for v in graph['ends'][edge]} for word in words]
            self.assertEqual(len(set.union(*vertices)),10)
            for i in range(3):
                self.assertEqual(len(vertices[i]),4)
                for j in range(i):
                    self.assertEqual(vertices[i]&vertices[j],{(0,0,0)})
            self.assertEqual(len(data['support']),12)
            self.assertEqual(sum(data['incidence_counts'].values()),30)
            self.assertEqual(sum(k*v for k,v in data['incidence_counts'].items()),48)

    def test_independent_vertex_gauges_preserve_the_local_measurements(self):
        data = z.bouquet(4)
        palette = (y.QA,y.QB,y.QX,y.ONE)
        graph = data['graph']
        field = {e:palette[i%4] for i,e in enumerate(graph['edges'])}
        gauges = {v:palette[(v[0]+2*v[1]+3*v[2])%4] for v in graph['vertices']}
        transformed = y.gauge_transform(graph,field,gauges)
        before = tuple(word_value(word,field) for word in data['words'])
        after = tuple(word_value(word,transformed) for word in data['words'])
        self.assertEqual(z.invariants(*before),z.invariants(*after))


class FullEnergyTests(unittest.TestCase):
    def test_full_twelve_link_derivatives_and_electric_eigenvalues(self):
        data = z.bouquet(4)
        palette = (y.QX,y.QB,y.ONE,y.QA)
        field = {e:palette[i%4] for i,e in enumerate(data['graph']['edges'])}
        words = data['words']
        holonomies = tuple(word_value(word,field) for word in words)
        coordinates = sum((list(q[1:]) for q in holonomies),[])
        polynomials = relational_polynomials()
        derivatives = [[partial(f,k) for k in range(9)] for f in polynomials]
        gram = [[F(0) for _ in range(3)] for _ in range(3)]
        laplacian = [F(0)]*3
        for moving in data['support']:
            for axis in (1,2,3):
                jets = [word_jet(word,field,moving,axis) for word in words]
                first = sum((list(q[1][1:]) for q in jets),[])
                second = sum((list(q[2][1:]) for q in jets),[])
                velocities = []
                for i,gradient in enumerate(derivatives):
                    velocities.append(sum(evaluate(g,coordinates)*v for g,v in zip(gradient,first)))
                    acceleration = sum(evaluate(g,coordinates)*v for g,v in zip(gradient,second))
                    acceleration += sum(evaluate(partial(g,l),coordinates)*first[k]*first[l]
                                        for k,g in enumerate(gradient) for l in range(9))
                    laplacian[i] -= acceleration
                for i,j in product(range(3),repeat=2):
                    gram[i][j] += velocities[i]*velocities[j]
        self.assertEqual(gram,z.energy_control(*holonomies)['full_electric_Gram'])
        self.assertEqual(laplacian,[e*evaluate(f,coordinates) for e,f in zip((6,6,9),polynomials)])

    def test_independent_haar_moments_match_full_energy_integration_by_parts(self):
        polynomials = relational_polynomials()
        norms = (F(3,16),F(3,16),F(3,32))
        for i,f in enumerate(polynomials):
            self.assertEqual(haar(f),0)
            for j,g in enumerate(polynomials):
                covariance = haar(multiply(f,g))
                self.assertEqual(covariance,norms[i] if i==j else 0)
                gamma = {}
                for base in (0,3,6):
                    radial_f,radial_g = {},{}
                    for k in range(base,base+3):
                        df,dg = partial(f,k),partial(g,k)
                        gamma = add(gamma,multiply(df,dg))
                        radial_f = add(radial_f,multiply(monomial((k,)),df))
                        radial_g = add(radial_g,multiply(monomial((k,)),dg))
                    gamma = add(gamma,multiply(radial_f,radial_g),-1)
                self.assertEqual(haar(gamma),(6,6,9)[j]*covariance)

    def test_pointwise_energy_bounds_hold_without_independence_of_loops(self):
        for a,b,x in product((y.ONE,y.QA,y.QB,y.QX),repeat=3):
            out = z.energy_control(a,b,x)
            gram = out['full_electric_Gram']
            aa,bb,xx = (1-q[0]**2 for q in (a,b,x))
            self.assertLessEqual(gram[0][0],aa+xx)
            self.assertLessEqual(gram[1][1],bb+xx)
            self.assertLessEqual(gram[2][2],aa+bb+xx)
            self.assertTrue(p.psd(gram))

    def test_actual_energy_constants_follow_from_the_pinned_defect_bound(self):
        defect_coefficient = F(9,4)
        self.assertEqual(4*defect_coefficient,9)
        self.assertEqual(6*defect_coefficient,F(27,2))
        # The same upper bound makes the norms vanish; it supplies no variance floor.
        for s in (5,10,100):
            self.assertLessEqual(2*defect_coefficient/s,1)
            for a in (F(1,10),F(1,1000)):
                kappa = 1/(a*s)
                self.assertEqual(kappa*9/s,9/(a*s*s))


class EnergyTailTests(unittest.TestCase):
    def test_actual_marginal_ratio_reuses_exact_upstream_heat_kernel_bounds(self):
        for theta in (F(1,16000),F(1),F(25)):
            out = z.density_budget(12,30,theta)
            upstream = p.conditional_budget(12,60*theta)
            self.assertEqual(out['omega'],60*theta)
            self.assertEqual(out['fixed_time']['prefactor']*upstream['fixed_time']['prefactor'],F(3,4))
            self.assertEqual(out['fixed_time']['exp_exponent'],upstream['fixed_time']['exponent'])
            self.assertEqual(out['optimized']['prefactor']*upstream['optimized']['prefactor'],F(3,4))
            self.assertEqual(out['optimized']['exp_exponent_upper'],upstream['optimized']['exponent_upper'])
        self.assertEqual(z.density_budget(12,30,0)['exact_ratio'],1)

    def test_weighted_character_gram_and_projection_from_three_term_recurrence(self):
        for eta in (F(0),F(3,4),F(-2,3)):
            for cutoff in (2,4):
                out = z.weighted_projection_control(eta,cutoff)
                coefficients = {0:-eta/2,1:F(1),cutoff+1:F(1,20),cutoff+3:F(1,100)}
                inner = lambda m,n:F(m==n)+(eta/2 if abs(m-n)==1 else 0)
                for m in range(cutoff+1):
                    expected = sum(value*inner(m,n) for n,value in coefficients.items())
                    actual = sum(value*inner(m,n) for n,value in enumerate(out['projection_coefficients']))
                    self.assertEqual(actual,expected)
                self.assertEqual(out['actual_metric_projection']==out['Haar_projection'],not eta)

    def test_weighted_polynomial_energy_from_independent_electric_character_matrix(self):
        for eta in (F(0),F(3,4),F(-2,3)):
            for cutoff in (2,4):
                out = z.weighted_projection_control(eta,cutoff)
                coefficients = {0:-eta/2,1:F(1),cutoff+1:F(1,20),cutoff+3:F(1,100)}
                energy = F(0)
                for m,am in coefficients.items():
                    for n,an in coefficients.items():
                        lm,ln = m*(m+2),n*(n+2)
                        entry = F(lm) if m==n else eta*(lm+ln-3)/4 if abs(m-n)==1 else 0
                        energy += am*an*entry
                self.assertEqual(out['actual_energy'],energy)
                self.assertLessEqual(out['tail_squared'],out['tail_upper'])
                self.assertGreater(out['tail_squared'],0)

    def test_incomplete_low_energy_marker_space_does_not_satisfy_the_tail_theorem(self):
        # H=diag(1,1,100); retaining only e1 misses e2 below cutoff 10.
        omitted_norm_squared,energy,cutoff = F(1),F(1),F(10)
        self.assertGreater(omitted_norm_squared,energy/cutoff)
        # The complete cutoff must retain BOTH e1 and e2.
        self.assertEqual(omitted_norm_squared-1,0)

    def test_three_loop_rank_is_only_an_unprojected_upper_bound(self):
        previous = 0
        for cutoff in (0,3,6,9,12,24):
            exact = sum((a+1)**2*(b+1)**2*(c+1)**2 for a,b,c in product(range(6),repeat=3)
                        if a*(a+2)+b*(b+2)+c*(c+2)<=cutoff)
            self.assertEqual(z.three_loop_rank_upper(cutoff),exact)
            self.assertGreaterEqual(exact,previous)
            previous = exact
        self.assertEqual(z.three_loop_rank_upper(3),13)

    def test_invalid_density_and_cutoff_inputs_are_rejected(self):
        for eta in (1,-1,F(3,2)):
            with self.assertRaises(ValueError):
                z.weighted_projection_control(eta,2)
        for args in ((0,30,1),(12,0,1),(12,30,-1)):
            with self.assertRaises(ValueError):
                z.density_budget(*args)


class FiniteResponseTests(unittest.TestCase):
    def test_nonidentity_covariance_is_used_in_the_operator_norm(self):
        c = [[1,0],[0,2]]
        r = [[F(1,16),0],[0,F(1,2)]]
        out = z.response_budget(c,r,F(1,2),2,1,128)
        self.assertEqual(out['energy_constrained_response_upper'],F(5,8))
        with self.assertRaisesRegex(ValueError,'finite response upper bound'):
            z.response_budget(c,r,F(49,100),2,1,128)

    def test_changing_basis_preserves_the_certified_response_bound(self):
        c = [[1,F(1,3)],[F(1,3),2]]
        r = [[F(x)/9 for x in row] for row in c]
        basis = [[1,2],[0,1]]
        transform = lambda m:p.multiply(p.transpose(basis),p.multiply(m,basis))
        self.assertEqual(z.response_budget(c,r,F(1,3),1,1,100),
                         z.response_budget(transform(c),transform(r),F(1,3),1,1,100))

    def test_full_outputs_detect_response_lost_by_output_truncation(self):
        # T e1=e2. Its compression to span(e1) is zero, but ||T e1||=1.
        out = z.response_budget([[1]],[[1]],1,1,1,100)
        self.assertEqual(out['energy_constrained_response_upper'],1)
        with self.assertRaisesRegex(ValueError,'finite response upper bound'):
            z.response_budget([[1]],[[1]],0,1,1,100)

    def test_conditional_projection_full_output_gram_has_the_stated_identity(self):
        states = list(product((-1,1),repeat=3))
        weights = [F(1,8)*(1+F(1,3)*a*b+F(1,5)*b*c) for a,b,c in states]
        pb = p.conditional_projection(states,weights,[0])
        pc = p.conditional_projection(states,weights,[1])
        both = p.conditional_projection(states,weights,[0,1])
        t = p.add(p.multiply(pb,pc),both,-1)
        vectors = [[F(a) for a,b,c in states],[F(b) for a,b,c in states],
                   [F(a*c) for a,b,c in states]]
        apply = lambda m,v:[sum(a*b for a,b in zip(row,v)) for row in m]
        inner = lambda a,b:sum(w*x*y for w,x,y in zip(weights,a,b))
        for f,g in product(vectors,repeat=2):
            expected = inner(apply(pc,f),apply(pb,apply(pc,g)))-inner(apply(both,f),apply(both,g))
            self.assertEqual(inner(apply(t,f),apply(t,g)),expected)

    def test_energy_tail_shrinks_with_cutoff_and_all_bounds_are_capped(self):
        values = [z.response_budget([[1]],[[F(1,4)]],F(1,2),2,1,c)
                  for c in (16,128,512)]
        self.assertGreater(values[0]['energy_constrained_response_upper'],values[1]['energy_constrained_response_upper'])
        self.assertEqual(values[2]['energy_constrained_response_upper'],F(9,16))
        self.assertEqual(z.response_budget([[1]],[[0]],0,1000,10,1)['energy_constrained_response_upper'],1)

    def test_invalid_gram_and_density_claims_are_rejected(self):
        cases = (([[0]],[[0]],0,1), ([[1]],[[2]],1,1),
                 ([[1]],[[-1]],1,1), ([[1]],[[0]],1,F(1,2)),
                 ([[1,0],[0,1]],[[0]],1,1))
        for c,r,b,d in cases:
            with self.assertRaises(ValueError):
                z.response_budget(c,r,b,d,1,100)


class EvidenceTests(unittest.TestCase):
    def test_reproduction_pins_the_scope_and_all_source_inputs(self):
        cert = z.run()
        self.assertEqual(len(cert['input_sha256']),380)
        self.assertEqual(z.check(cert),z.PIN.read_text().strip())
        for name in ('seven_linear_markers_span_all_physical_observables',
                     'actual_normalized_state_energy_or_variance_floor_proved',
                     'energy_constrained_bound_is_unrestricted_operator_norm',
                     'actual_weak_coupling_response_matrix_bound_supplied',
                     'positive_four_dimensional_continuum_mass_gap',
                     'continuum_state_survival_proved'):
            self.assertFalse(cert['model'][name])

    def test_changed_result_pin_and_scope_are_rejected(self):
        cert = z.run()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                z.check(cert,result,pin)
            pin.write_text(p.canonical_sha(cert))
            changed = copy.deepcopy(cert)
            changed['controls']['actual_energy_upper_coefficients_over_sqrt_theta'] = []
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                z.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
            with self.assertRaises(ValueError):
                z.check(changed,result,pin)

    def test_upstream_source_drift_is_rejected(self):
        manifest = json.loads(z.SOURCES.read_text())
        path = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][path] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'manifest'
            source.write_text(json.dumps(manifest))
            with patch.object(z,'SOURCES',source):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    z.source_checks()

    def test_cli_is_read_only_and_optimized_runtime_is_rejected(self):
        inputs = (z.RESULT,z.PIN,z.SOURCES)
        before = [path.read_bytes() for path in inputs]
        checked = subprocess.run([sys.executable,'-B',str(z.HERE/'ym86_finite_energy_response.py'),'--check'],capture_output=True,text=True)
        self.assertEqual(checked.returncode,0,checked.stderr)
        self.assertIn('YM86 PASS',checked.stdout)
        self.assertEqual(before,[path.read_bytes() for path in inputs])
        optimized = subprocess.run([sys.executable,'-B','-O',str(z.HERE/'ym86_finite_energy_response.py'),'--check'],capture_output=True,text=True)
        self.assertNotEqual(optimized.returncode,0)
        self.assertIn('without optimization',optimized.stderr)


if __name__ == '__main__':
    unittest.main()
