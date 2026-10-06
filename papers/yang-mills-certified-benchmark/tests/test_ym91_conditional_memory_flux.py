import copy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym91_conditional_memory_flux as y

p,a = y.p,y.a


def diagonal(values):
    return [[v if i==j else F(0) for j in range(len(values))] for i,v in enumerate(values)]


class ConditionalTensorTests(unittest.TestCase):
    def test_inverse_order_in_noncommuting_conditional_and_hidden_operators(self):
        out = y.matrix_control()
        A,D,R = out['conditional_operator'],out['hidden_operator'],out['scores']
        self.assertNotEqual(p.multiply(A,D),p.multiply(D,A))
        for s in (F(0),F(1,7),F(3),F(100)):
            case = y.memory_comparison(A,D,R,s)
            self.assertTrue(p.psd(case['dynamic_memory']))
            self.assertTrue(p.psd(p.add(case['conditional_tensor'],case['dynamic_memory'],-1)))
            direct = p.multiply(p.transpose(R),p.multiply(a.inverse(p.add(D,a.scale(p.identity(3),s))),R))
            self.assertEqual(case['dynamic_memory'],direct)

    def test_conditional_tensors_need_not_commute_at_different_frequencies(self):
        out = y.matrix_control()
        T0,T1 = [case['conditional_tensor'] for case in out['frequencies'][:2]]
        self.assertNotEqual(p.multiply(T0,T1),p.multiply(T1,T0))
        self.assertTrue(p.psd(p.add(T0,T1,-1)))

    def test_shifted_full_residual_enclosure_transfers_only_its_upper_bound(self):
        out = y.matrix_control()
        for case,residual in zip(out['frequencies'],out['conditional_residual_enclosures']):
            low,high = residual['enclosure']['lower'],residual['enclosure']['upper']
            self.assertEqual(residual['exact_energy_dual_score'],case['conditional_tensor'])
            self.assertTrue(p.psd(p.add(case['conditional_tensor'],low,-1)))
            self.assertTrue(p.psd(p.add(high,case['conditional_tensor'],-1)))
            self.assertTrue(p.psd(p.add(high,case['dynamic_memory'],-1)))

    def test_exact_conditional_solution_still_is_not_a_dynamic_lower_bound(self):
        out = y.conditional_lower_is_not_dynamic_lower()
        self.assertEqual(out['exact_conditional_trial_lower'],[[F(1)]])
        self.assertEqual(out['dynamic_memory'],[[F(1,100)]])
        self.assertFalse(p.psd(p.add(out['dynamic_memory'],out['exact_conditional_trial_lower'],-1)))

    def test_source_basis_congruence_preserves_all_cross_terms(self):
        out = y.matrix_control()
        rotation = a.matrix([[F(4,5),F(-3,5)],[F(3,5),F(4,5)]])
        changed = y.memory_comparison(out['conditional_operator'],out['hidden_operator'],p.multiply(out['scores'],rotation),1)
        original = y.memory_comparison(out['conditional_operator'],out['hidden_operator'],out['scores'],1)
        for key in ('conditional_tensor','dynamic_memory'):
            self.assertEqual(changed[key],p.multiply(p.transpose(rotation),p.multiply(original[key],rotation)))

    def test_actual_small_theta_bound_is_recovered_without_extending_its_window(self):
        for theta in (0,F(1,800),F(1,400)):
            old = y.previous.actual_budget(theta,3)
            d = old['conditional_floor']
            a2 = old['cross_Hessian_norm_squared']
            self.assertEqual(4*a2/(d*d),old['relative_memory'])
        self.assertEqual(y.previous.actual_budget(F(1,400),3)['relative_memory'],F(1,144))
        with self.assertRaises(ValueError):
            y.previous.actual_budget(F(1,399),3)

    def test_operator_order_symmetry_and_frequency_hypotheses_are_enforced(self):
        for A,D in (([[1]],[[F(1,2)]]),([[0]],[[1]]),([[1,1],[0,1]],[[2,1],[1,2]])):
            with self.assertRaises(ValueError):
                y.memory_comparison(A,D,[[1]] if len(A)==1 else [[1],[1]])
        with self.assertRaises(ValueError):
            y.memory_comparison([[1]],[[2]],[[1]],-1)


class FluxTests(unittest.TestCase):
    def test_every_flux_has_exact_full_divergence_and_an_upper_Gram(self):
        out = y.matrix_control()
        N,R = out['gradient'],out['scores']
        for case in out['flux_hierarchy']:
            self.assertEqual(p.multiply(p.transpose(N),case['corrected_flux']),R)
            self.assertTrue(p.psd(p.add(case['upper_Gram'],case['exact_conditional_tensor'],-1)))

    def test_nested_divergence_free_corrections_monotonically_reach_the_minimum(self):
        cases = y.matrix_control()['flux_hierarchy']
        for old,new in zip(cases,cases[1:]):
            improvement = p.add(old['upper_Gram'],new['upper_Gram'],-1)
            self.assertTrue(p.psd(improvement))
            self.assertNotEqual(improvement,a.zeros(2))
        self.assertEqual(cases[-1]['upper_Gram'],cases[-1]['exact_conditional_tensor'])

    def test_corrected_flux_is_orthogonal_to_all_supplied_cycles(self):
        out = y.matrix_control()
        J = out['flux_hierarchy'][-1]['corrected_flux']
        self.assertEqual(p.multiply(p.transpose(out['cycles']),J),a.zeros(2))
        expected = p.multiply(out['gradient'],p.multiply(a.inverse(out['conditional_operator']),out['scores']))
        self.assertEqual(J,expected)

    def test_cycle_basis_change_preserves_the_corrected_certificate(self):
        out = y.matrix_control()
        changed = p.multiply(out['cycles'],a.matrix([[1,2],[0,3]]))
        result = y.flux_certificate(out['gradient'],out['scores'],out['initial_flux'],changed)
        self.assertEqual(result['upper_Gram'],out['flux_hierarchy'][-1]['upper_Gram'])

    def test_nonzero_divergence_residual_and_invalid_cycles_are_rejected(self):
        out = y.matrix_control()
        bad = copy.deepcopy(out['initial_flux'])
        bad[0][0] += F(1,1000)
        with self.assertRaisesRegex(ValueError,'full divergence residual'):
            y.flux_certificate(out['gradient'],out['scores'],bad)
        with self.assertRaisesRegex(ValueError,'divergence free'):
            y.flux_certificate(out['gradient'],out['scores'],out['initial_flux'],[[1],[0],[0],[0],[0]])
        duplicate = [[row[0],row[0]] for row in out['cycles']]
        with self.assertRaisesRegex(ValueError,'independent'):
            y.flux_certificate(out['gradient'],out['scores'],out['initial_flux'],duplicate)

    def test_exact_flux_can_stay_small_when_an_unused_conditional_mode_is_slow(self):
        for epsilon in (F(1,4),F(1,100),F(1,10000)):
            N = diagonal([epsilon,F(1)])
            R = [[F(0)],[F(1,4)]]
            result = y.flux_certificate(N,R,R)
            self.assertEqual(result['upper_Gram'],[[F(1,16)]])
            crude = F(1,16)/(epsilon*epsilon)
            self.assertGreater(crude,result['upper_Gram'][0][0])


class NativeComparisonTests(unittest.TestCase):
    def test_character_normalization_parity_and_native_eigenvalue(self):
        for n in (1,3,5,9,15):
            coefficients = y.character_coefficients(n)
            polynomial = {(k,0,0,0):v for k,v in enumerate(coefficients) if v}
            self.assertEqual(y.evaluate(coefficients,F(1)),1)
            self.assertEqual(y.evaluate(coefficients,F(-1)),-1)
            self.assertTrue(all(k%2 for k,v in enumerate(coefficients) if v))
            self.assertEqual(y.local.y.phi(polynomial),0)
            self.assertEqual(y.local.y.phi(y.poly.mul(polynomial,polynomial)),F(1,(n+1)**2))
            residual = y.poly.add(y.poly.lap(polynomial),y.poly.scale(polynomial,-F(n*(n+2),4)))
            self.assertEqual(y.poly.sphere_reduce(residual),{})

    def test_native_memory_decreases_while_curvature_and_Hessian_test_deteriorate(self):
        previous_eta,previous_ricci = F(1),F(1)
        for n in (1,3,5,9,31,101):
            out = y.native_budget(n)
            self.assertEqual(out['relative_memory_upper'],F(2,15*(n+1)**2+1))
            self.assertLess(out['relative_memory_upper'],previous_eta)
            self.assertLess(out['pointwise_weighted_Ricci_at_minus_one'],previous_ricci)
            if n>=5:
                self.assertTrue(out['YM90_Hessian_gate_impossible'])
                self.assertLess(out['pointwise_weighted_Ricci_at_minus_one'],0)
            self.assertGreater(out['physical_assembled_gap_lower'],0)
            previous_eta,previous_ricci = out['relative_memory_upper'],out['pointwise_weighted_Ricci_at_minus_one']

    def test_actual_native_cross_Hessian_and_conditional_gap_upper_are_certified(self):
        for n in (1,3,5,9):
            out = y.native_control(n)
            self.assertEqual(out['Rayleigh_energy']/out['Rayleigh_variance'],F(3,4))
            self.assertEqual(out['central_cross_Hessian_Gram'],a.scale(p.identity(3),out['budget']['cross_Hessian_norm_squared_lower']))
            self.assertEqual(out['budget']['conditional_floor'],F(9,20))
        self.assertEqual(y.native_budget(5)['best_possible_YM90_Hessian_gate_ratio_lower'],F(2450,729))

    def test_conditional_covariant_Hessian_at_the_negative_center_has_the_stated_sign(self):
        W,_ = y.previous.parent_polynomials()
        one = y.local.y.ONE
        qs = (one,tuple(-v for v in one),one,one,one)
        for n in (1,5,9):
            coeffs = y.character_coefficients(n)
            chi_prime = F(1,4)*y.evaluate(y.derivative(coeffs),F(-1))/(2*F(3,4))
            raw = [[chi_prime*y.local.model.evaluate(y.native.derivative(y.native.derivative(W,4,k),4,l),qs)
                    for k in (1,2,3)] for l in (1,2,3)]
            H = a.scale(p.add(raw,p.transpose(raw)),F(1,2))
            Ricci = p.add(a.scale(p.identity(3),F(1,2)),a.scale(H,2),-1)
            self.assertEqual(Ricci,a.scale(p.identity(3),y.native_budget(n)['pointwise_weighted_Ricci_at_minus_one']))

    def test_native_word_jets_independently_check_the_physical_source(self):
        turns = y.local.y.turns()
        for n in (1,3,5):
            for offset in (0,2):
                qs = tuple(turns[(offset+2*e)%len(turns)] for e in range(5))
                point = y.native_flux_point(n,F(1,4),qs)
                values = dict(enumerate(qs))
                gamma = sum(y.previous.adapted.word_derivative(values,y.previous.PARENT_W_WORD,e,k)*
                            y.previous.adapted.word_derivative(values,y.previous.PARENT_F_WORD,e,k)
                            for e in range(5) for k in (1,2,3))
                df = y.evaluate(y.derivative(y.character_coefficients(n)),point['W'])
                self.assertEqual(point['source'],-F(1,4)*df/point['density']*gamma)

    def test_engineered_parent_potential_matches_direct_chain_differentiation(self):
        W,_ = y.previous.parent_polynomials()
        qs = tuple(y.local.y.turns()[:5])
        for n in (1,3,9):
            point = y.native_flux_point(n,F(1,4),qs)
            coeffs = y.character_coefficients(n)
            df = y.evaluate(y.derivative(coeffs),point['W'])
            ddf = y.evaluate(y.derivative(y.derivative(coeffs)),point['W'])
            rho = point['density']
            first = F(1,4)*df/(2*rho)
            second = F(1,4)*ddf/(2*rho)-F(1,16)*df*df/(4*rho*rho)
            electric = y.local.model.evaluate(y.local.model.electric(W,5),qs)
            gamma = sum(y.local.model.evaluate(y.native.derivative(W,e,k),qs)**2 for e in range(5) for k in (1,2,3))
            self.assertEqual(first*electric-second*gamma,point['parent_potential'])

    def test_flux_divergence_is_exact_for_noncommuting_points_and_certified_constants(self):
        qs = tuple(y.local.y.turns()[:5])
        for n in (1,5,9):
            constant = y.reciprocal_integral(n)['certified_flux_constant']
            point = y.native_flux_point(n,F(1,4),qs,constant)
            for row in point['flux_components']:
                self.assertEqual(row['score'],row['negative_weighted_divergence'])

    def test_physical_source_is_nonzero_and_gauge_invariant(self):
        out = y.native_control(1)
        self.assertEqual(out['reference']['source'],-F(8,115))
        self.assertEqual(out['gauge_vertices_checked'],4)
        for n in (3,5,9):
            self.assertNotEqual(y.native_control(n)['reference']['source'],0)

    def test_full_reciprocal_series_enclosures_nest_and_improve_the_closed_bound(self):
        for n in (1,3,9):
            last = y.reciprocal_integral(n,degree=0)
            for degree in (1,2,4,8):
                out = y.reciprocal_integral(n,degree=degree)
                self.assertGreaterEqual(out['reciprocal_integral']['lower'],last['reciprocal_integral']['lower'])
                self.assertLessEqual(out['reciprocal_integral']['upper'],last['reciprocal_integral']['upper'])
                self.assertLess(out['reciprocal_integral']['lower'],out['reciprocal_integral']['upper'])
                self.assertLessEqual(2*out['flux_defect']['upper'],y.native_budget(n)['relative_memory_upper'])
                last = out
            self.assertLess(last['full_tail_upper'],F(1,10**10))

    def test_finite_reciprocal_upper_endpoint_gives_an_actually_feasible_flux_cost(self):
        for n in (1,5,9):
            out = y.reciprocal_integral(n,degree=3)
            c = out['certified_flux_constant']
            high = out['reciprocal_integral']['upper']
            self.assertEqual(1-2*c+c*c*high,out['flux_defect']['upper'])
            refined = y.reciprocal_integral(n,degree=8)['reciprocal_integral']
            for value in refined.values():
                self.assertLessEqual(1-2*c+c*c*value,out['flux_defect']['upper'])
            self.assertGreater(out['flux_defect']['lower'],0)

    def test_free_density_and_failed_large_block_gate_are_not_overclaimed(self):
        free = y.native_budget(9,0)
        self.assertEqual(free['relative_memory_upper'],0)
        integral = y.reciprocal_integral(9,0)
        self.assertEqual(integral['reciprocal_integral'],dict(lower=1,upper=1))
        failed = y.native_budget(1,F(1,4),factor=100)
        self.assertFalse(failed['strict_memory_gate'])
        self.assertEqual(failed['physical_assembled_gap_lower'],0)

    def test_invalid_native_orders_densities_and_points_are_rejected(self):
        for n in (0,2,4,True):
            with self.assertRaises(ValueError):
                y.native_budget(n)
        for epsilon in (-1,1,2):
            with self.assertRaises(ValueError):
                y.native_budget(3,epsilon)
        with self.assertRaises(TypeError):
            y.native_budget(3,0.25)
        with self.assertRaises(ValueError):
            y.native_flux_point(1,F(1,4),[(1,1,0,0)]*5)


class FailureAndScopeTests(unittest.TestCase):
    def test_small_observed_memory_does_not_control_a_hidden_invisible_mode(self):
        for epsilon in (F(1,4),F(1,16),F(1,256)):
            out = y.invisible_hidden_mode(epsilon)
            H = out['operator']
            self.assertEqual(p.multiply(H,[[0],[0],[1]]),[[0],[0],[epsilon]])
            self.assertEqual(a.inverse(H)[0][0],F(4,3))
            self.assertEqual(out['relative_memory'],F(1,4))

    def test_exterior_average_misses_a_true_periodic_gradient(self):
        for M in (4,16,64):
            out = y.averaged_exterior_counterexample(F(1,M*M),M)
            self.assertEqual(sum(w*s for w,s in zip(out['arc_lengths'],out['slopes'])),0)
            self.assertEqual(out['exterior_average'],F(1,M))
            self.assertEqual(out['energy_relative_response'],M)

    def test_zero_frequency_flux_control_does_not_supply_a_decay_rate(self):
        old = y.previous.insufficient_zero_frequency_control()
        self.assertGreater(old['shifted_memory'],old['unsupported_frequency_upper'])
        # Fixed inverse-source norm and arbitrarily large time-zero source norm.
        for rate in (1,100,10000):
            self.assertEqual(F(rate,4)/rate,F(1,4))
            self.assertGreaterEqual(F(rate,4),F(1,4))


class EvidenceTests(unittest.TestCase):
    def test_sources_scope_and_result_reproduce(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']),415)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        self.assertTrue(cert['model']['actual_finite_Wilson_conditional_memory_domination'])
        for name in ('actual_Wilson_small_theta_window_extended','actual_weak_coupling_memory_gate_discharged',
                     'positive_four_dimensional_continuum_mass_gap','native_comparison_is_original_Wilson_model',
                     'lower_conditional_tensor_bounds_dynamic_memory_from_below','analytic_proofs_formally_machine_verified'):
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
            changed['controls']['native_budgets'][0]['relative_memory_upper'] = '0'
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['native_comparison_is_original_Wilson_model'] = True
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
        script = str(y.HERE/'ym91_conditional_memory_flux.py')
        for args in ([],['--check']):
            run = subprocess.run([sys.executable,'-B',script,*args],capture_output=True,text=True)
            self.assertEqual(run.returncode,0,run.stderr)
            self.assertIn('YM91 PASS',run.stdout)
            self.assertEqual(before,[path.read_bytes() for path in inputs])
        run = subprocess.run([sys.executable,'-B','-O',script,'--check'],capture_output=True,text=True)
        self.assertNotEqual(run.returncode,0)
        self.assertIn('without optimization',run.stderr)


if __name__=='__main__':
    unittest.main()
