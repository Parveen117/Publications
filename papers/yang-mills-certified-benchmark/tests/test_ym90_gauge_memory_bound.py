import copy
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym90_gauge_memory_bound as y

p,a = y.p,y.algebra


def diag(values):
    return [[v if i==j else F(0) for j in range(len(values))] for i,v in enumerate(values)]


def full_matrix(C,B,D):
    Bt = p.transpose(B)
    return [row+Bt[i] for i,row in enumerate(C)]+[row+D[i] for i,row in enumerate(B)]


def moment_fixture():
    out = y.matrix_control()
    return out['C'],out['B_hidden_by_coarse'],out['D']


class ActualBudgetTests(unittest.TestCase):
    def test_endpoint_constants_and_free_limit(self):
        out = y.actual_budget(F(1,400),2)
        self.assertEqual(out['conditional_floor'],F(3,10))
        self.assertEqual(out['cross_Hessian_norm_squared'],F(1,6400))
        self.assertEqual(out['source_energy_factor'],F(1,480))
        self.assertEqual(out['relative_memory'],F(1,144))
        self.assertEqual(out['resolvent_relative_form_error'],F(1,143))
        free = y.actual_budget(0,2)
        for name in ('source_energy_factor','relative_memory','resolvent_relative_form_error'):
            self.assertEqual(free[name],0)

    def test_block_factors_do_not_amplify_relative_memory(self):
        for b in (2,3,4,8,16,64):
            out = y.actual_budget(F(1,400),b)
            self.assertEqual(out['relative_memory'],F(1,144))
            self.assertGreaterEqual(out['coarse_diagonal_gap_lower'],b*F(227,480))
            self.assertLess(out['physical_two_block_gap_lower'],out['existing_YM88_direct_gap_lower'])
            self.assertEqual(out['unrestricted_two_block_gap_lower'],F(1,16))

    def test_frequency_error_decreases_and_has_the_correct_high_frequency_scale(self):
        previous = F(1)
        for s in (0,F(1,10),1,10,100):
            out = y.actual_budget(F(1,400),3,s)
            alpha = out['memory_frequency_factor']
            self.assertEqual(alpha,F(1,480)/(F(3,10)+s))
            self.assertLess(alpha,previous)
            self.assertLessEqual(out['resolvent_relative_form_error'],F(1,143))
            previous = alpha

    def test_conditional_floor_is_derived_from_the_pointwise_Hessian(self):
        last = F(-1)
        for theta in (0,F(1,10000),F(1,2000),F(1,800),F(1,400)):
            fine = y.previous.previous.parameters(theta)
            out = y.actual_budget(theta,2)
            self.assertEqual(out['conditional_floor'],F(1,2)-2*(fine['diagonal_log_Hessian_upper']+5*theta))
            self.assertGreater(out['relative_memory'],last)
            last = out['relative_memory']

    def test_source_inequality_does_not_by_itself_assert_a_subunit_gate(self):
        out = y.source_budget(1,1)
        self.assertEqual(out['relative_memory'],4)
        with self.assertRaises(ValueError):
            y.previous.two_block_floor(1,1,out['relative_memory'])

    def test_invalid_windows_sizes_and_inexact_inputs_are_rejected(self):
        for theta in (-1,F(1,399)):
            with self.assertRaises(ValueError):
                y.actual_budget(theta,2)
        for b in (0,1,True,F(2)):
            with self.assertRaises(ValueError):
                y.actual_budget(0,b)
        for theta in (True,0.001):
            with self.assertRaises(TypeError):
                y.actual_budget(theta,2)
        with self.assertRaises(ValueError):
            y.actual_budget(0,2,-1)
        with self.assertRaises(ValueError):
            y.source_budget(1,0)


class GaugeSectorTests(unittest.TestCase):
    def test_actual_network_incidence_and_native_charge_floor(self):
        for b in (2,3,4):
            out = y.network_control(b)
            self.assertEqual(out['max_degree'],6)
            self.assertLessEqual(out['max_internal_endpoints'],2)
            self.assertEqual(out['universal_charged_gap_lower'],F(1,16))
            self.assertEqual(out['internal_vertices'],3*4**3*(b-1))

    def test_transported_path_copy_has_norm_squared_b_not_b_squared(self):
        turns = y.local.y.turns()
        for b in (2,3,7):
            T = [[F(0)]*6 for _ in range(6*b)]
            for path in range(2):
                prefix = y.local.y.ONE
                for e in range(b):
                    for k in range(3):
                        unit = tuple(F(i==k) for i in range(3))
                        rotated = y.previous.adjoint(y.local.y.qdagger(prefix),unit)
                        for axis,value in enumerate(rotated):
                            T[3*(path*b+e)+axis][3*path+k] = value
                    prefix = y.local.y.qmul(prefix,turns[(e+path)%len(turns)])
            self.assertEqual(p.multiply(p.transpose(T),T),a.scale(p.identity(6),b))
            # Cross-coordinate derivative operator with norm exactly 1/5.
            J = [[F(i==j,5) for j in range(6*b)] for i in range(6)]
            N = p.multiply(J,T)
            self.assertTrue(p.psd(p.add(a.scale(p.identity(6),F(b,25)),p.multiply(p.transpose(N),N),-1)))
            self.assertEqual(F(4*b,25)/F(3,10)/b,F(8,15))

    def test_overlapping_paths_break_the_disjoint_copy_identity(self):
        T = [[F(1),F(1)],[F(1),F(0)],[F(0),F(1)]]
        gram = p.multiply(p.transpose(T),T)
        self.assertFalse(p.psd(p.add(a.scale(p.identity(2),2),gram,-1)))

    def test_source_selection_allows_the_invariant_floor(self):
        out = y.sector_control()
        self.assertEqual(out['invariant_source_memory'],F(1,64))
        self.assertEqual(out['charged_source_memory'],F(1,4))
        A = diag([F(1),F(1),F(0)])
        for key,commutes in (('B_invariant',True),('B_charged',False)):
            H = full_matrix([[F(1)]],out[key],out['D'])
            self.assertEqual(p.multiply(A,H)==p.multiply(H,A),commutes)
            self.assertTrue(p.psd(H))


class DynamicMemoryTests(unittest.TestCase):
    def test_projected_full_inverse_matches_the_memory_resolvent(self):
        C,B,D = moment_fixture()
        H = full_matrix(C,B,D)
        for case in y.matrix_control()['frequencies']:
            inverse = a.inverse(p.add(H,a.scale(p.identity(5),case['damping'])))
            self.assertEqual(y.adapted.sub(inverse,[0,1]),case['compressed_resolvent'])

    def test_resolvent_error_order_is_checked_as_a_matrix_not_entrywise(self):
        for case in y.matrix_control()['frequencies']:
            error = p.add(case['compressed_resolvent'],case['Markov_resolvent'],-1)
            upper = a.scale(case['Markov_resolvent'],case['relative_form_error'])
            self.assertTrue(p.psd(error))
            self.assertTrue(p.psd(p.add(upper,error,-1)))

    def test_moment_bounds_and_noncommutation(self):
        out = y.matrix_control()
        for n,M in enumerate(out['moments']):
            self.assertTrue(p.psd(M))
            self.assertTrue(p.psd(p.add(a.scale(out['C'],F(1,12)*F(1,2)**n),M,-1)))
        self.assertNotEqual(out['first_moment_commutator'],a.zeros(2))
        self.assertEqual(p.transpose(out['first_moment_commutator']),a.scale(out['first_moment_commutator'],-1))

    def test_stieltjes_Hankel_forms_remain_positive_without_commuting_moments(self):
        moments = y.matrix_control()['moments']
        for shift in (0,1):
            self.assertTrue(p.psd(y.adapted.hankel(moments,3,shift)))

    def test_finite_signed_remainders_hold_beyond_the_power_series_radius(self):
        C,B,D = moment_fixture()
        Bt = p.transpose(B)
        for s in (F(0),F(1,3),F(2),F(7)):
            Sigma = y.adapted.self_energy(Bt,D,s)
            partial = a.zeros(2)
            for n in range(6):
                partial = p.add(partial,a.scale(y.adapted.moment(Bt,D,n),(-s)**n))
                signed = a.scale(p.add(Sigma,partial,-1),(-1)**(n+1))
                bound = y.remainder_budget(F(1,12),2,s,n)
                self.assertTrue(p.psd(signed))
                self.assertTrue(p.psd(p.add(a.scale(C,bound),signed,-1)))

    def test_analytic_remainder_controls_both_signs_inside_its_radius(self):
        C,B,D = moment_fixture()
        for s in (F(-3,2),F(-1,3),F(1,3),F(3,2)):
            partial = a.zeros(2)
            Sigma = y.adapted.self_energy(p.transpose(B),D,s)
            for n in range(5):
                partial = p.add(partial,a.scale(y.adapted.moment(p.transpose(B),D,n),(-s)**n))
                error = p.add(Sigma,partial,-1)
                bound = a.scale(C,y.remainder_budget(F(1,12),2,s,n,analytic=True))
                for sign in (-1,1):
                    self.assertTrue(p.psd(p.add(bound,error,sign)))
        for s in (-2,2,3):
            with self.assertRaises(ValueError):
                y.remainder_budget(F(1,12),2,s,1,analytic=True)
        with self.assertRaises(ValueError):
            y.remainder_budget(F(1,12),2,-1,1)

    def test_time_kernel_and_integrated_tail_with_certified_exponential_intervals(self):
        C,B,D = moment_fixture()
        for t in (F(0),F(1,3),F(2)):
            exps = [p.exp_negative_interval(D[i][i]*t) for i in range(3)]
            floor_exp = p.exp_negative_interval(2*t)['lower']
            kernel_upper = p.multiply(p.transpose(B),p.multiply(diag([e['upper'] for e in exps]),B))
            tail_upper = p.multiply(p.transpose(B),p.multiply(diag([e['upper']/D[i][i] for i,e in enumerate(exps)]),B))
            self.assertTrue(p.psd(p.add(a.scale(C,F(1,6)*floor_exp),kernel_upper,-1)))
            self.assertTrue(p.psd(p.add(a.scale(C,F(1,12)*floor_exp),tail_upper,-1)))

    def test_clock_scaling_preserves_eta_and_rescales_frequency(self):
        C,B,D = moment_fixture()
        for kappa in (F(1,7),F(2),F(13)):
            for s in (F(0),F(1,3),F(4)):
                scaled = y.matrix_case(a.scale(C,kappa),a.scale(B,kappa),a.scale(D,kappa),2*kappa,kappa/6,s)
                original = y.matrix_case(C,B,D,2,F(1,6),s/kappa)
                self.assertEqual(scaled['memory'],a.scale(original['memory'],kappa))
                self.assertEqual(scaled['frequency_factor'],original['frequency_factor'])
                self.assertEqual(scaled['compressed_resolvent'],a.scale(original['compressed_resolvent'],1/kappa))

    def test_two_block_floor_against_the_full_operator(self):
        C,B,D = moment_fixture()
        c,d,eta = F(5,2),F(2),F(1,12)
        self.assertTrue(p.psd(p.add(C,a.scale(p.identity(2),c),-1)))
        floor = y.previous.two_block_floor(c,d,eta)
        self.assertTrue(p.psd(p.add(full_matrix(C,B,D),a.scale(p.identity(5),floor),-1)))

    def test_zero_frequency_memory_alone_does_not_give_the_stronger_frequency_bound(self):
        out = y.insufficient_zero_frequency_control()
        self.assertEqual(out['zero_frequency_memory'],F(1,4))
        self.assertGreater(out['shifted_memory'],out['unsupported_frequency_upper'])
        with self.assertRaisesRegex(ValueError,'Source energy'):
            y.matrix_case([[1]],[[5]],[[100]],1,F(1,4),1)

    def test_invalid_hidden_and_source_certificates_are_rejected(self):
        C,B,D = moment_fixture()
        for d,sigma in ((3,F(1,6)),(2,F(1,1000)),(2,2)):
            with self.assertRaises(ValueError):
                y.matrix_case(C,B,D,d,sigma,1)
        with self.assertRaises(ValueError):
            y.matrix_case([[1,1],[0,1]],B,D,2,F(1,6),1)


class NativeParentTests(unittest.TestCase):
    def test_exact_Gauss_identities_and_nonzero_physical_source(self):
        out = y.parent_native_control()
        self.assertEqual(out['Gauss_polynomial_identities'],24)
        self.assertEqual(out['reference']['source'],-F(8,125))
        self.assertEqual(out['gauge_vertices_checked'],[0,1,2,3])

    def test_independent_word_jets_check_the_source_on_noncommuting_configurations(self):
        turns = y.local.y.turns()
        for offset in range(6):
            qs = tuple(turns[(offset+2*e)%len(turns)] for e in range(5))
            values = dict(enumerate(qs))
            gamma = sum(y.adapted.word_derivative(values,y.PARENT_W_WORD,e,k)*
                        y.adapted.word_derivative(values,y.PARENT_F_WORD,e,k)
                        for e in range(5) for k in (1,2,3))
            out = y.parent_point(qs,F(1,10))
            self.assertEqual(out['source'],-F(1,5)*gamma)

    def test_engineered_parent_ground_state_equation_uses_native_Casimir(self):
        W,_ = y.parent_polynomials()
        turns = y.local.y.turns()
        for offset in range(4):
            qs = tuple(turns[(offset+e)%len(turns)] for e in range(5))
            tau = F(1,10)
            out = y.parent_point(qs,tau)
            minus_delta_chi = tau*y.local.model.evaluate(y.local.model.electric(W,5),qs)
            gradient_square = sum((tau*y.local.model.evaluate(y.native.derivative(W,e,k),qs))**2
                                  for e in range(5) for k in (1,2,3))
            self.assertEqual(minus_delta_chi-gradient_square,out['parent_potential'])

    def test_native_conditional_Hessian_and_cross_operator_norm(self):
        W,_ = y.parent_polynomials()
        tau = F(1,10)
        points = [(y.local.y.ONE,)*5,tuple(y.local.y.turns()[:5])]
        for qs in points:
            evaluate = lambda f:tau*y.local.model.evaluate(f,qs)
            J = [[evaluate(y.native.derivative(y.native.derivative(W,e,k),4,l))
                  for e in range(4) for k in (1,2,3)] for l in (1,2,3)]
            bound = a.scale(p.identity(3),tau*tau/8)
            product = p.multiply(J,p.transpose(J))
            self.assertTrue(p.psd(p.add(bound,product,-1)))
            if qs==(y.local.y.ONE,)*5:
                self.assertEqual(product,bound)
            raw = [[evaluate(y.native.derivative(y.native.derivative(W,4,k),4,l))
                    for k in (1,2,3)] for l in (1,2,3)]
            hessian = a.scale(p.add(raw,p.transpose(raw)),F(1,2))
            self.assertEqual(hessian,a.scale(p.identity(3),-evaluate(W)/4))

    def test_Haar_moments_independently_derive_both_positive_series(self):
        for n in range(10):
            moment = y.local.y.native_moment((2*n,0,0,0))
            next_moment = y.local.y.native_moment((2*n+2,0,0,0))
            self.assertEqual(F(2**(2*n),factorial(2*n))*moment,F(1,factorial(n)*factorial(n+1)))
            self.assertEqual(F(2**(2*n),factorial(2*n))*(moment-next_moment),F(3,2*factorial(n)*factorial(n+2)))
        for i in range(1,4):
            for j in range(1,4):
                self.assertEqual(y.local.y.phi(y.poly.mul(y.poly.COORD[i],y.poly.COORD[j])),F(i==j,4))

    def test_fixed_coarse_rotations_leave_Haar_integral_coefficients_unchanged(self):
        for q in y.local.y.turns()[:4]:
            W = y.poly.add(*(y.poly.scale(y.poly.COORD[i],q[i]) for i in range(4)))
            power = y.poly.ONE
            for n in range(9):
                self.assertEqual(y.local.y.phi(power),y.local.y.native_moment((n,0,0,0)))
                power = y.poly.mul(power,W)

    def test_positive_series_tail_enclosures_are_nested_and_strict(self):
        for tau in (F(1,10),F(1,4)):
            for weighted in (False,True):
                last = y.positive_parent_series(tau,0,weighted)
                for order in range(1,11):
                    now = y.positive_parent_series(tau,order,weighted)
                    self.assertLess(last['lower'],now['lower'])
                    self.assertLess(now['upper'],last['upper'])
                    self.assertLess(now['lower'],now['upper'])
                    last = now
                self.assertLess(last['upper']-last['lower'],F(1,10**28))

    def test_Haar_coarse_law_has_positive_memory_source_with_full_integral_bounds(self):
        out = y.parent_integral_control(F(1,10))
        self.assertTrue(out['coarse_marginal_is_exact_product_Haar'])
        self.assertGreater(out['source_norm_squared']['lower'],0)
        self.assertEqual(out['conditional_floor'],F(9,20))
        self.assertEqual(out['source_energy_factor'],F(1,90))
        self.assertEqual(out['relative_memory'],F(2,81))
        self.assertLess(out['source_norm_squared']['upper'],out['source_energy_factor']*out['probe_coarse_energy'])
        free = y.parent_integral_control(0)
        self.assertEqual(free['source_norm_squared'],dict(lower=0,upper=0))
        self.assertFalse(free['source_nonzero'])

    def test_actual_Wilson_leakage_witness_is_kept_distinct_from_parent(self):
        out = y.adapted.native_controls()
        self.assertEqual(out['first_order_hidden_leakage'],-F(32,75))
        self.assertEqual(out['fine_boundary_links'],12)
        self.assertTrue(out['full_vacuum_not_replaced_by_this_first_order_coefficient'])

    def test_parent_controls_reject_invalid_parameters_and_nonunit_quaternions(self):
        for tau in (-1,F(1,3)):
            with self.assertRaises(ValueError):
                y.parent_integral_control(tau)
        with self.assertRaises(ValueError):
            y.parent_point([(1,1,0,0)]*5,F(1,10))


class EvidenceTests(unittest.TestCase):
    def test_source_chain_and_actual_scope_reproduce(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']),408)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        self.assertTrue(cert['model']['actual_small_theta_dynamic_memory_gate_discharged'])
        for name in ('actual_weak_coupling_dynamic_memory_gate_discharged',
                     'positive_four_dimensional_continuum_mass_gap',
                     'coarse_static_Haar_law_forces_zero_memory',
                     'coupling_is_bounded_on_coarse_L2','resolvent_order_implies_semigroup_order',
                     'engineered_parent_is_original_Wilson_model','analytic_proofs_formally_machine_verified'):
            self.assertFalse(cert['model'][name])

    def test_result_pin_and_scope_drift_are_rejected(self):
        cert = y.run()
        with tempfile.TemporaryDirectory() as directory:
            result,pin = Path(directory)/'result',Path(directory)/'pin'
            result.write_text(json.dumps(cert))
            pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            pin.write_text(p.canonical_sha(cert))
            changed = copy.deepcopy(cert)
            changed['controls']['actual_endpoint']['relative_memory'] = '0'
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
            with self.assertRaisesRegex(ValueError,'scope mismatch'):
                y.check(changed,result,pin)

    def test_upstream_hash_drift_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        key = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][key] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'sources'
            path.write_text(json.dumps(manifest))
            with patch.object(y,'SOURCES',path):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    y.source_checks()

    def test_cli_is_read_only_and_rejects_optimization(self):
        inputs = (y.RESULT,y.PIN,y.SOURCES)
        before = [path.read_bytes() for path in inputs]
        script = str(y.HERE/'ym90_gauge_memory_bound.py')
        for args in ([],['--check']):
            run = subprocess.run([sys.executable,'-B',script,*args],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            self.assertIn('YM90 PASS',run.stdout)
            self.assertEqual(before,[path.read_bytes() for path in inputs])
        run = subprocess.run([sys.executable,'-B','-O',script,'--check'],capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('without optimization',run.stderr)


if __name__=='__main__':
    unittest.main()
