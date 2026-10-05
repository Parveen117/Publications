import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym88_spatial_score_contraction as y

p,poly,old = y.p,y.poly,y.previous


def dec(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


def weighted_norm(matrix):
    # The path graph on these matrix indices; independent of the lattice builder.
    return max(sum(2**abs(i-j)*abs(v) for j,v in enumerate(row))
               for i,row in enumerate(matrix))


def walk_partial(source,theta,degree):
    out = old.zeros(len(source))
    power = p.identity(len(source))
    for coefficient in y.catalan_coefficients(theta,degree):
        power = p.multiply(power,source)
        out = p.add(out,old.scale(power,coefficient))
    return out


class LocalBoundsTests(unittest.TestCase):
    def test_enlarged_endpoint_and_free_limit(self):
        self.assertEqual(y.THETA_MAX/F(1,8000),20)
        endpoint = y.endpoint_control()['endpoint']
        self.assertEqual(endpoint['conditional_gap_lower'],F(13,40))
        self.assertEqual(endpoint['global_bare_gap_lower'],F(3,10))
        self.assertEqual(endpoint['weighted_transfer_row_upper'],F(2,13))
        self.assertEqual(endpoint['unweighted_transfer_row_upper'],F(1,13))
        free = y.parameters(0)
        for key in ('wave_log_gradient_upper','diagonal_log_Hessian_upper',
                    'weighted_mixed_Hessian_row_upper','weighted_transfer_row_upper'):
            self.assertEqual(free[key],0)
        self.assertEqual(free['global_bare_gap_lower'],F(1,2))

    def test_exponential_enclosures_against_independent_high_precision_values(self):
        with localcontext() as context:
            context.prec = 90
            for theta in (F(0),F(1,100000),F(1,8000),F(1,2000),F(1,1000),F(1,500),y.THETA_MAX):
                cell = y.parameters(theta)
                M = dec(4*theta)
                contraction = (9*M).exp()*(dec(F(1,20))+9*M)
                B = dec(F(9,2))*M*dec(F(21,19)**2)*(45*M).exp()/(1-contraction)
                self.assertLessEqual((9*M).exp(),dec(cell['exp9M_upper']))
                self.assertLessEqual((45*M).exp(),dec(cell['exp45M_upper']))
                self.assertLessEqual(B,dec(cell['wave_log_gradient_upper']))
                self.assertLessEqual(B,dec(42*theta))
                self.assertLessEqual(B/2+3*B*B,dec(35*theta))

    def test_inherited_derivative_formula_agrees_inside_its_original_window(self):
        theta = F(1,16000)
        inherited = y.local.parameters(theta)['wave_log_gradient']
        with localcontext() as context:
            context.prec = 180
            M = dec(4*theta)
            actual = dec(F(9,2))*M*dec(F(21,19)**2)*(45*M).exp()/(1-(9*M).exp()*(dec(F(1,20))+9*M))
            self.assertLessEqual(dec(inherited.lo),actual)
            self.assertLessEqual(actual,dec(inherited.hi))
            self.assertLessEqual(actual,dec(y.parameters(theta)['wave_log_gradient_upper']))
        with self.assertRaises(ValueError):
            y.local.parameters(y.THETA_MAX)

    def test_window_is_not_extrapolated_to_weak_coupling(self):
        for theta in (-1,F(1,399),1,25):
            with self.assertRaises(ValueError):
                y.parameters(theta)
        for theta in (True,0.001,None):
            with self.assertRaises(TypeError):
                y.parameters(theta)

    def test_quadratic_large_branch_is_real_and_continuation_is_essential(self):
        for u in (0,F(1,40),F(3,10),1,100):
            self.assertGreaterEqual(y.barrier_polynomial(u),0)
        for u in (F(1,30),F(1,10),F(1,4)):
            self.assertLess(y.barrier_polynomial(u),0)
        self.assertGreater(F(3,10),10*y.THETA_MAX)

    def test_small_branch_implies_linear_bound_and_invariant_ball(self):
        for theta in (F(0),F(1,8000),F(1,1000),y.THETA_MAX):
            cell = y.parameters(theta)
            b,c = 10*theta,cell['walk_resolvent_factor']
            self.assertLessEqual(140*theta+F(2,40),F(2,5))
            self.assertLessEqual(c*(6*theta+2*b*b),b)
            self.assertLessEqual(4*c*b,F(2,13))
            self.assertLessEqual(48*c*c*theta,F(48,169))


class NativeGeometryTests(unittest.TestCase):
    def test_native_bracket_sign_on_polynomial_coordinates(self):
        fixtures = [y.variable(i,2) for i in range(8)]
        fixtures.append(poly.mul(fixtures[0],poly.mul(fixtures[2],fixtures[5])))
        for value in fixtures:
            for i,j in permutations(range(3),2):
                bracket = poly.add(y.D(y.D(value,j),i),poly.scale(y.D(y.D(value,i),j),-1))
                expected = poly.add(*(poly.scale(y.D(value,k),y.geometry.epsilon(i,j,k)) for k in range(3)))
                self.assertEqual(bracket,expected)
            self.assertEqual(y.D(y.D(value,3),0),y.D(y.D(value,0),3))

    def test_covariant_Hessian_is_symmetric_while_ordered_derivatives_are_not(self):
        chi = y.variable(0,2)
        H = y.hessian(chi,2)
        self.assertEqual(H,p.transpose(H))
        self.assertNotEqual(y.D(y.D(chi,1),0),y.D(y.D(chi,0),1))
        self.assertEqual(H[0][1],poly.add(y.D(y.D(chi,1),0),poly.scale(y.D(chi,2),F(-1,2))))

    def test_Ricci_and_Casimir_have_different_native_constants(self):
        out = y.geometry.lie_ricci_control()
        self.assertEqual(out['Ricci'],old.scale(p.identity(3),F(1,2)))
        self.assertEqual(out['fundamental_Casimir'],old.scale(p.identity(4),F(3,4)))
        for i in range(8):
            x = y.variable(i,2)
            lap = poly.add(*(y.D(y.D(x,k),k) for k in range(6)))
            self.assertEqual(lap,poly.scale(x,F(-3,4)))

    def test_exact_two_and_three_factor_tensor_equations(self):
        for links,count in ((2,9),(3,27)):
            out = y.tensor_control(links)
            self.assertEqual(out['exact_mixed_components_checked'],count)
            self.assertEqual(out['residuals_nonzero_if_Ricci_is_dropped'],count)
            self.assertEqual(out['components_changed_by_ordinary_instead_of_rough_Laplacian'],count)
            self.assertEqual(out['mixed_tensor_damping'],1)
            self.assertIn('not_actual_Wilson',out['kind'])

    def test_native_Wilson_linearization_and_quaternion_double_insertions(self):
        W = y.local.model.loop_polynomial()
        self.assertEqual(y.local.model.electric(W,4),poly.scale(W,3))
        qs = tuple(tuple(F(x) for x in v) for v in ((F(3,5),F(4,5),0,0),
                    (0,F(5,13),F(12,13),0),(F(1,2),F(1,2),F(1,2),F(1,2)),(0,0,0,1)))
        self.assertTrue(all(sum(x*x for x in q)==1 for q in qs))
        u,v = (F(3,5),F(4,5),F(0)),(F(0),F(5,13),F(12,13))
        for e,f in permutations(range(4),2):
            derivative = y.local.derivative(y.local.derivative(W,f,v),e,u)
            expected = y.local.word_jet(qs,(e,u),(f,v))
            self.assertEqual(y.local.model.evaluate(derivative,qs),expected)
            self.assertLessEqual(abs(expected),F(1,4))

    def test_control_carrier_validation(self):
        for n in (1,4,True):
            with self.assertRaises(ValueError):
                y.tensor_control(n)
        with self.assertRaises(ValueError):
            y.variable(8,2)


class SpatialWalkTests(unittest.TestCase):
    def test_periodic_graph_incidence_has_volume_independent_row_sum(self):
        for side in (4,5,6):
            rows,distances = y.lattice_graph(side)
            self.assertEqual(len(rows),3*side**3)
            for i,row in enumerate(rows):
                self.assertNotIn(i,row)
                self.assertEqual(sum(row.values()),12)
                for j,multiplicity in row.items():
                    self.assertEqual(rows[j][i],multiplicity)
                    self.assertLessEqual(abs(distances[i]-distances[j]),1)
            self.assertEqual(sum(2**distances[j]*v for j,v in rows[0].items()),24)

    def test_walk_support_and_return_diagonal_are_retained(self):
        control = y.graph_control(4,6)
        steps = control['walk_controls']
        self.assertEqual(steps[0]['diagonal_return_multiplicity'],0)
        self.assertGreater(steps[1]['diagonal_return_multiplicity'],0)
        for n,row in enumerate(steps,1):
            self.assertEqual(row['total_walk_multiplicity'],12**n)
            self.assertLessEqual(row['weighted_walk_sum'],24**n)

    def test_weighted_norm_is_submultiplicative_on_independent_path_matrices(self):
        X = [[F(0),F(1,2),F(-1,7)],[F(2,3),F(1),F(0)],[F(0),F(1,9),F(2)]]
        Y = [[F(1),F(0),F(1,4)],[F(-3),F(0),F(1,3)],[F(1,2),F(1),F(0)]]
        self.assertLessEqual(weighted_norm(p.multiply(X,Y)),weighted_norm(X)*weighted_norm(Y))
        self.assertGreater(weighted_norm(X),max(sum(abs(v) for v in row) for row in X))

    def test_Catalan_coefficients_satisfy_independent_quadratic_convolution(self):
        for theta in (0,F(1,8000),y.THETA_MAX):
            c = y.parameters(theta)['walk_resolvent_factor']
            coefficients = y.catalan_coefficients(theta,14)
            self.assertEqual(coefficients[0],c)
            for n in range(2,15):
                expected = 2*c*sum(coefficients[j-1]*coefficients[n-j-1] for j in range(1,n))
                self.assertEqual(coefficients[n-1],expected)

    def test_walk_series_and_full_tail_against_exact_two_vertex_matrix_function(self):
        theta = y.THETA_MAX
        c = y.parameters(theta)['walk_resolvent_factor']
        a = 3*theta
        C = [[0,a],[a,0]]
        with localcontext() as context:
            context.prec = 100
            fun = lambda z:(1-(1-8*dec(c*c)*z).sqrt())/(4*dec(c))
            plus,minus = fun(dec(a)),fun(-dec(a))
            diagonal,off = (plus+minus)/2,(plus-minus)/2
            self.assertLessEqual(diagonal+2*off,dec(10*theta))
            for degree in (0,1,2,4,8,16):
                partial = walk_partial(C,theta,degree)
                diagonal_error,off_error = diagonal-dec(partial[0][0]),off-dec(partial[0][1])
                self.assertGreater(diagonal_error,0)
                self.assertGreater(off_error,0)
                self.assertLessEqual(diagonal_error+2*off_error,dec(y.walk_tail(theta,degree)))
                if degree>=2:
                    self.assertGreater(partial[0][0],0)

    def test_finite_walk_truncation_is_not_used_as_an_upper_bound(self):
        theta = y.THETA_MAX
        C = [[0,3*theta],[3*theta,0]]
        lower = walk_partial(C,theta,1)
        later = walk_partial(C,theta,5)
        self.assertLess(lower[0][0],later[0][0])
        self.assertLess(lower[0][1],later[0][1])
        self.assertLessEqual(weighted_norm(p.add(later,lower,-1)),y.walk_tail(theta,1))

    def test_Picard_iteration_is_enclosed_by_the_certified_contraction(self):
        theta = y.THETA_MAX
        c = y.parameters(theta)['walk_resolvent_factor']
        C = [[0,3*theta],[3*theta,0]]
        X = old.zeros(2)
        b,q = 10*theta,4*c*10*theta
        previous_increment = None
        for _ in range(7):
            next_X = old.scale(p.add(C,old.scale(p.multiply(X,X),2)),c)
            increment = weighted_norm(p.add(next_X,X,-1))
            self.assertLessEqual(weighted_norm(next_X),b)
            self.assertTrue(all(next_X[i][j]>=X[i][j] for i in range(2) for j in range(2)))
            if previous_increment is not None:
                self.assertLessEqual(increment,q*previous_increment)
            X,previous_increment = next_X,increment

    def test_endpoint_score_tail_includes_the_distance_one_factor(self):
        self.assertEqual(y.score_tail(y.THETA_MAX,1),F(1,13))
        for r in (1,2,4,10,30):
            self.assertEqual(y.score_tail(y.THETA_MAX,r),F(2,13*2**r))
            self.assertEqual(y.score_tail(0,r),0)
        self.assertEqual(y.walk_tail(y.THETA_MAX,0),F(39,1210))
        self.assertEqual(y.walk_tail(0,0),0)

    def test_invalid_tail_and_periodic_inputs_are_rejected(self):
        for distance in (0,-1,True,F(1,2)):
            with self.assertRaises(ValueError):
                y.score_tail(y.THETA_MAX,distance)
        for degree in (-1,True,F(1,2)):
            with self.assertRaises(ValueError):
                y.walk_tail(y.THETA_MAX,degree)
        with self.assertRaises(ValueError):
            y.lattice_graph(3)


class SweepClosureTests(unittest.TestCase):
    def setUp(self):
        self.K = [[0,F(1,26),F(1,26)],[F(1,26),0,F(1,52)],[F(1,26),F(1,52),0]]
        self.rho = F(13,40)

    def test_metric_matches_independently_solved_full_Lyapunov_system(self):
        for order in permutations(range(3)):
            cert = y.sweep_metric(self.K,self.rho,order)
            solved = old.lyapunov_solution(cert['transfer'],cert['loss'])
            self.assertEqual(cert['metric'],solved)
            self.assertEqual(cert['conditional_bare_gap_lower'],F(3,10))

    def test_each_update_loses_exactly_its_conditional_variance_budget(self):
        cert = y.sweep_metric(self.K,self.rho)
        data = old.sweep_data(self.K,[self.rho]*3)
        for i,update in enumerate(data['updates']):
            difference = p.add(cert['metric'],p.multiply(p.transpose(update),p.multiply(cert['metric'],update)),-1)
            expected = old.zeros(3)
            expected[i][i] = 1/self.rho
            self.assertEqual(difference,expected)

    def test_global_metric_bound_agrees_with_direct_block_Ricci_bound(self):
        for theta in (0,F(1,8000),F(1,1000),y.THETA_MAX):
            cell = y.parameters(theta)
            gap = cell['conditional_gap_lower']*(1-cell['unweighted_transfer_row_upper'])
            self.assertEqual(gap,cell['global_bare_gap_lower'])
            self.assertEqual(gap,F(1,2)-2*(35*theta+10*theta/2))

    def test_free_metric_bound_is_dimension_independent(self):
        for n in (1,3,7):
            cert = y.sweep_metric(old.zeros(n),F(1,2))
            self.assertEqual(cert['metric'],old.scale(p.identity(n),2))
            self.assertEqual(cert['conditional_bare_gap_lower'],F(1,2))

    def test_invalid_or_noncontracting_metrics_are_rejected(self):
        for K in ([[0,1],[1,0]],[[F(1,10),0],[0,0]],[[0,F(1,10)],[F(1,20),0]],[[0,-1],[-1,0]]):
            with self.assertRaises(ValueError):
                y.sweep_metric(K,1)
        with self.assertRaises(ValueError):
            y.sweep_metric(self.K,0)


class EvidenceTests(unittest.TestCase):
    def test_reproduction_pins_sources_and_limits_claims_to_the_proved_window(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']),394)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        for name in ('actual_Wilson_vacuum_spatial_Hessian_bound',
                     'actual_full_exterior_conditional_score_summability',
                     'actual_small_coupling_sweep_gate_discharged'):
            self.assertTrue(cert['model'][name])
        for name in ('actual_weak_coupling_sweep_gate_discharged',
                     'continuum_state_survival_proved','positive_four_dimensional_continuum_mass_gap',
                     'analytic_proofs_formally_machine_verified',
                     'finite_walk_truncation_alone_is_an_upper_bound'):
            self.assertFalse(cert['model'][name])

    def test_changed_result_pin_and_scope_are_rejected(self):
        cert = y.run()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            pin.write_text(p.canonical_sha(cert))
            changed = copy.deepcopy(cert)
            changed['controls']['endpoint']['window_endpoint_ratio'] = 200
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
            with self.assertRaisesRegex(ValueError,'scope mismatch'):
                y.check(changed,result,pin)

    def test_upstream_source_drift_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        key = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][key] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'manifest'
            source.write_text(json.dumps(manifest))
            with patch.object(y,'SOURCES',source):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    y.source_checks()

    def test_default_and_check_are_read_only_and_optimization_is_rejected(self):
        inputs = (y.RESULT,y.PIN,y.SOURCES)
        before = [path.read_bytes() for path in inputs]
        script = str(y.HERE/'ym88_spatial_score_contraction.py')
        for args in ([],['--check']):
            checked = subprocess.run([sys.executable,'-B',script,*args],capture_output=True,text=True)
            self.assertEqual(checked.returncode,0,checked.stderr)
            self.assertIn('YM88 PASS',checked.stdout)
            self.assertEqual(before,[path.read_bytes() for path in inputs])
        optimized = subprocess.run([sys.executable,'-B','-O',script,'--check'],capture_output=True,text=True)
        self.assertNotEqual(optimized.returncode,0)
        self.assertIn('without optimization',optimized.stderr)


if __name__=='__main__':
    unittest.main()
