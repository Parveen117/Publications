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
import ym87_conditional_energy_transport as y

p,old,poly = y.p,y.previous,y.poly
q = old.y


def apply(a,v):
    return [sum(x*z for x,z in zip(row,v)) for row in a]


def quadratic(a,v):
    return sum(x*z for x,z in zip(v,apply(a,v)))


def evaluate(values,x):
    result = F(0)
    for coefficient in reversed(values):
        result = result*x+coefficient
    return result


def decimal(x):
    return Decimal(x.numerator)/Decimal(x.denominator)


def linear_Gaussian_update(coefficients,precision,index):
    out = list(coefficients)
    value = out[index]
    out[index] = F(0)
    for j in range(len(out)):
        if j!=index:
            out[j] -= precision[index][j]*value/precision[index][index]
    return out


def precision_transfer(precision):
    n = len(precision)
    return [[abs(precision[j][i])/precision[i][i] if i!=j else F(0)
             for i in range(n)] for j in range(n)]


class ActualScoreBudgetTests(unittest.TestCase):
    def test_heat_smoothing_coefficient_has_a_strict_independent_margin(self):
        with localcontext() as context:
            context.prec = 60
            self.assertLess(Decimal(9)/2*(Decimal(1)/2).exp(),8)
            for m in (F(1,4000),F(1),F(100)):
                time = 1/(4*m)
                prefactor = (decimal(m*time).exp() *
                             ((2*decimal(time)).sqrt()**(-1)+decimal(m)*(2*decimal(time)).sqrt()))
                self.assertLess(prefactor*prefactor,8*decimal(m))
        self.assertLess(F(3),F(256,81))

    def test_actual_score_budgets_match_both_upstream_density_and_gap_bounds(self):
        for theta in (F(1,16000),F(1),F(25)):
            out = y.actual_score_budget(theta,4,13,4)
            density = old.density_budget(1,4,theta)
            gap = p.conditional_budget(4,26*theta)
            self.assertEqual(out['one_link_potential_bound']*out['heat_time'],F(1,4))
            for entry in out['wave_log_gradient_squared_upper_options']:
                d = density[entry['method']]
                exponent = d.get('exp_exponent',d.get('exp_exponent_upper'))
                self.assertEqual(entry['prefactor'],32*theta*d['prefactor']**2)
                self.assertEqual(entry['exp_exponent_upper'],2*exponent)
            for entry in out['conditional_energy_transfer_squared_upper_options']:
                d,r = density[entry['density_method']],gap[entry['gap_method']]
                de = d.get('exp_exponent',d.get('exp_exponent_upper'))
                re = r.get('exponent',r.get('exponent_upper'))
                self.assertEqual(entry['prefactor']*r['prefactor'],512*theta*d['prefactor']**2)
                self.assertEqual(entry['exp_exponent_upper'],2*de+re)
            self.assertFalse(out['spatial_decay_or_sweep_contraction_certified'])

    def test_free_score_vanishes_and_invalid_block_budgets_are_rejected(self):
        out = y.actual_score_budget(0,4,13,4)
        self.assertEqual(out['wave_log_gradient_squared_upper'],0)
        self.assertEqual(out['conditional_energy_transfer_squared_upper'],0)
        self.assertEqual(out['exact_free_conditional_gap'],F(3,4))
        for args in ((-1,1,4,1),(1,0,4,1),(1,1,0,1),(1,1,4,0)):
            with self.assertRaises(ValueError):
                y.actual_score_budget(*args)


class ConditionalEnergyTests(unittest.TestCase):
    def test_conditional_derivative_includes_the_normalization_score(self):
        # Positive comparison density 1+eta*b+eps*a*T_n(b), with Haar a.
        eta,eps,b = F(1,5),F(1,4),F(3,5)
        t = y.chebyshev(5)
        value,derivative = evaluate(t,b),evaluate(old.derivative(t),b)
        marginal = 1+eta*b
        density = [1,eps*value/marginal]
        f,df = [b*b,1,b],[2*b,0,1]
        avg = lambda v:poly.haar_integral(poly.poly_multiply(density,v))
        mean = avg(f)
        raw_score_term = poly.haar_integral(poly.poly_multiply(f,[eta,eps*derivative]))/marginal
        centered_score_term = raw_score_term-mean*eta/marginal
        expected = eps*((1+eta*b)*derivative-eta*value)/(4*marginal*marginal)+2*b+F(1,4)
        self.assertEqual(avg(df)+centered_score_term,expected)
        self.assertNotEqual(avg(df)+raw_score_term,expected)

    def test_native_Chebyshev_is_the_trace_of_an_ordered_SU2_power(self):
        carrier = q.unit((F(3,5),F(4,5),0,0))
        value = q.ONE
        for n in range(1,18):
            value = q.qmul(value,carrier)
            self.assertEqual(value[0],evaluate(y.chebyshev(n),carrier[0]))

    def test_four_link_energy_matches_independent_quaternion_power_derivatives(self):
        carrier = q.QX
        for n in (3,5,9):
            total = F(0)
            for axis in (1,2,3):
                direction = q.native_direction(carrier,axis)
                value,first = q.ONE,(F(0),)*4
                for _ in range(n):
                    first = q.plus(q.qmul(first,carrier),q.qmul(value,direction))
                    value = q.qmul(value,carrier)
                total += 4*first[0]*first[0]
            derivative = evaluate(old.derivative(y.chebyshev(n)),carrier[0])
            self.assertEqual(total,(1-carrier[0]**2)*derivative**2)

    def test_character_orthogonality_independently_fixes_native_energy(self):
        for n in (3,5,9,17):
            character_energy = F(n*(n+2)+(n-2)*n,4)
            self.assertEqual(character_energy,F(n*n,2))
            out = y.gauge_energy_counterexample(n,F(2,3))
            self.assertEqual(out['output_full_energy'],F(4,9)*character_energy/16)
            self.assertEqual(out['output_variance'],F(4,9)/32)

    def test_density_ratio_L2_response_and_gap_floor_stay_fixed_as_energy_grows(self):
        rows = [y.gauge_energy_counterexample(n) for n in (3,5,9,17)]
        for key in ('density_ratio','conditional_L2_correlation_squared',
                    'input_full_energy','full_and_conditional_Poincare_floor'):
            self.assertEqual(len({row[key] for row in rows}),1)
        self.assertEqual(rows[0]['density_ratio'],3)
        self.assertEqual(rows[0]['conditional_L2_correlation_squared'],F(1,32))
        self.assertEqual(rows[0]['full_and_conditional_Poincare_floor'],F(1,4))
        self.assertGreater(rows[-1]['energy_amplification'],1)
        self.assertLess(rows[0]['energy_amplification'],rows[-1]['energy_amplification'])

    def test_energy_dual_score_bound_is_saturated_at_the_declared_exterior(self):
        eps = F(1,2)
        for n in (3,5,9,17):
            out = y.gauge_energy_counterexample(n,eps)
            squared_conditional_derivative = eps*eps*n*n/16
            self.assertEqual(out['energy_dual_score_norm_at_declared_exterior']*out['input_full_energy'],
                             squared_conditional_derivative)

    def test_invalid_comparison_densities_and_noncentered_orders_are_rejected(self):
        for n,eps in ((2,F(1,2)),(3,0),(3,1),(3,F(3,2))):
            with self.assertRaises(ValueError):
                y.gauge_energy_counterexample(n,eps)


class PoissonResidualTests(unittest.TestCase):
    def test_Galerkin_lower_matrix_requires_the_full_omitted_residual_for_an_upper(self):
        out = y.poisson_control([[3,0,0],[0,8,0],[0,0,15]],[[1,2],[3,1],[2,-1]],
                               [[F(1,3),F(2,3)],[F(3,8),F(1,8)],[0,0]],3)
        self.assertEqual(out['cross'],y.zeros(2))
        self.assertNotEqual(out['exact_energy_dual_score'],out['enclosure']['lower'])
        exact_tail = y.scale(y.outer([2,-1]),F(1,15))
        self.assertEqual(p.add(out['exact_energy_dual_score'],out['enclosure']['lower'],-1),exact_tail)
        self.assertTrue(p.psd(p.add(out['enclosure']['upper'],out['exact_energy_dual_score'],-1)))

    def test_dropping_cross_terms_can_make_the_claimed_upper_bound_false(self):
        out = y.poisson_control([[3]],[[3]],[[F(1,2)]],3)
        naive_upper = out['energy_gram'][0][0]+out['residual_gram'][0][0]/3
        self.assertEqual(naive_upper,F(3,2))
        self.assertEqual(out['exact_energy_dual_score'],[[3]])
        self.assertEqual(out['enclosure']['upper'],[[3]])
        self.assertLess(naive_upper,3)

    def test_variational_lower_matrix_may_be_negative_without_invalidating_enclosure(self):
        out = y.poisson_control([[1]],[[1]],[[3]],1)
        self.assertEqual(out['enclosure']['lower'],[[-3]])
        self.assertEqual(out['exact_energy_dual_score'],[[1]])
        self.assertEqual(out['enclosure']['upper'],[[1]])

    def test_exact_solution_has_zero_residual_and_collapsed_enclosure(self):
        operator = [[3,1],[1,3]]
        scores = [[1,2],[3,-1]]
        trials = p.multiply(y.inverse(operator),scores)
        out = y.poisson_control(operator,scores,trials,2)
        self.assertEqual(out['residual_gram'],y.zeros(2))
        self.assertEqual(out['enclosure']['lower'],out['exact_energy_dual_score'])
        self.assertEqual(out['enclosure']['upper'],out['exact_energy_dual_score'])

    def test_residual_certificate_is_covariant_under_an_orthogonal_source_basis(self):
        operator,scores,trials = [[3,0],[0,8]],[[1,2],[3,-1]],[[F(1,3),F(2,3)],[0,0]]
        rotation = [[F(3,5),F(-4,5)],[F(4,5),F(3,5)]]
        transformed = p.multiply(rotation,p.multiply(operator,p.transpose(rotation)))
        before = y.poisson_control(operator,scores,trials,3)
        after = y.poisson_control(transformed,p.multiply(rotation,scores),p.multiply(rotation,trials),3)
        for key in ('exact_energy_dual_score','energy_gram','cross','residual_gram','enclosure'):
            self.assertEqual(before[key],after[key])

    def test_false_conditional_floor_and_incompatible_grams_are_rejected(self):
        with self.assertRaises(ValueError):
            y.poisson_control([[3]],[[1]],[[0]],4)
        with self.assertRaises(ValueError):
            y.poisson_control([[3]],[[1,2]],[[0]],3)
        with self.assertRaises(ValueError):
            y.residual_enclosure([[1]],[[0]],[[-1]],1)


class SweepGapTests(unittest.TestCase):
    def test_directed_updates_match_gaussian_conditional_means_in_every_order(self):
        precision = [[F(2),F(-1,5),F(-1,4)],
                     [F(-1,5),F(3),F(-1,6)],[F(-1,4),F(-1,6),F(4)]]
        coefficients = precision_transfer(precision)
        gaps = [precision[i][i] for i in range(3)]
        covariance = y.inverse(precision)
        for order in permutations(range(3)):
            data = y.sweep_data(coefficients,gaps,order)
            initial,current = [F(1),F(2),F(3)],[F(1),F(2),F(3)]
            total_loss = F(0)
            for i in order:
                total_loss += current[i]**2/gaps[i]
                current = linear_Gaussian_update(current,precision,i)
            self.assertEqual(current,apply(data['transfer'],initial))
            self.assertEqual(total_loss,quadratic(data['loss'],initial))
            self.assertEqual(total_loss,quadratic(covariance,initial)-quadratic(covariance,current))
            self.assertEqual(p.add(covariance,p.multiply(p.transpose(data['transfer']),p.multiply(covariance,data['transfer'])),-1),data['loss'])

    def test_signed_gradients_require_the_componentwise_triangle_bound(self):
        precision = [[F(2),F(-1,2)],[F(-1,2),F(1)]]
        coefficients = precision_transfer(precision)
        data = y.sweep_data(coefficients,[2,1])
        current = [F(1),F(-2)]
        for i,update in enumerate(data['updates']):
            bound = apply(update,[abs(x) for x in current])
            current = linear_Gaussian_update(current,precision,i)
            self.assertTrue(all(abs(a)<=b for a,b in zip(current,bound)))

    def test_gaussian_star_has_the_exact_gap_when_the_unweighted_row_test_fails(self):
        out = y.gaussian_star()
        self.assertEqual(out['unweighted_conditional_correlation_row'],F(6,5))
        self.assertEqual(out['exact_Gaussian_gap'],F(2,5))
        self.assertEqual(out['certificate']['beta'],F(5,2))
        self.assertEqual(out['certificate']['transfer'][0],[F(9,25)]+[F(3,10)]*4)
        self.assertTrue(p.psd(p.add(out['precision'],p.identity(5),F(-2,5))))
        self.assertEqual(apply(p.add(out['precision'],p.identity(5),F(-2,5)),out['sharp_gradient']),[0]*5)
        self.assertTrue(out['weighted_row_and_Hessian_criteria_can_also_succeed'])

    def test_independent_closed_star_covariance_and_sweep_series_agree(self):
        a,m = F(3,10),4
        denominator = 1-m*a*a
        expected = [[F(0) for _ in range(5)] for _ in range(5)]
        expected[0][0] = 1/denominator
        for j in range(1,5):
            expected[0][j] = expected[j][0] = a/denominator
            for k in range(1,5):
                expected[j][k] = F(j==k)+a*a/denominator
        out = y.gaussian_star()['certificate']
        self.assertEqual(out['metric'],expected)
        self.assertEqual(y.lyapunov_solution(out['transfer'],out['loss']),expected)
        partial,term = y.zeros(5),p.identity(5)
        for _ in range(6):
            partial = p.add(partial,p.multiply(p.transpose(term),p.multiply(out['loss'],term)))
            term = p.multiply(out['transfer'],term)
        remainder = p.multiply(p.transpose(term),p.multiply(expected,term))
        self.assertEqual(p.add(partial,remainder),expected)

    def test_nonnormal_sweep_can_pass_while_its_Euclidean_norm_exceeds_one(self):
        precision = [[F(32),F(-4)],[F(-4),F(1)]]
        k = precision_transfer(precision)
        covariance = y.inverse(precision)
        out = y.gap_certificate(k,[32,1],covariance,F(9,4))
        self.assertEqual(out['transfer'],[[F(1,2),4],[0,0]])
        self.assertGreater(sum(v*v for v in apply(out['transfer'],[0,1])),1)
        self.assertEqual(out['Lyapunov_remainder'],y.zeros(2))
        self.assertEqual(out['conditional_bare_gap_lower'],F(4,9))

    def test_free_certificate_is_volume_independent_and_uses_native_three_quarters(self):
        for n in (1,3,8):
            g = y.scale(p.identity(n),F(4,3))
            out = y.gap_certificate(y.zeros(n),[F(3,4)]*n,g,F(4,3))
            self.assertEqual(out['conditional_bare_gap_lower'],F(3,4))
            self.assertEqual(out['transfer'],y.zeros(n))

    def test_unsupported_lyapunov_and_metric_norm_claims_are_rejected(self):
        k = [[0,F(3,10)],[F(3,10),0]]
        data = y.sweep_data(k,[1,1])
        g = y.lyapunov_solution(data['transfer'],data['loss'])
        with self.assertRaisesRegex(ValueError,'norm upper bound'):
            y.gap_certificate(k,[1,1],g,F(1))
        with self.assertRaisesRegex(ValueError,'Lyapunov inequality'):
            y.gap_certificate(k,[1,1],y.scale(g,F(1,2)),10)
        with self.assertRaisesRegex(ValueError,'positive definite'):
            y.gap_certificate(k,[1,1],[[-1,0],[0,1]],10)
        with self.assertRaises(ValueError):
            y.gap_certificate([[0,2],[2,0]],[1,1],y.scale(p.identity(2),100),100)

    def test_invalid_orders_coefficients_and_shapes_are_rejected(self):
        for order in ((0,0),(0,),(-1,1),(True,0)):
            with self.assertRaises(ValueError):
                y.sweep_data([[0,1],[1,0]],[1,1],order)
        for k,gaps in (([[1,0],[0,0]],[1,1]),([[0,-1],[0,0]],[1,1]),
                       ([[0,0],[0,0]],[1]),([[0,0],[0,0]],[0,1])):
            with self.assertRaises(ValueError):
                y.sweep_data(k,gaps)


class EvidenceTests(unittest.TestCase):
    def test_reproduction_pins_all_sources_and_keeps_actual_gap_gate_open(self):
        cert = y.run()
        self.assertEqual(len(cert['input_sha256']),387)
        self.assertEqual(y.check(cert),y.PIN.read_text().strip())
        for name in ('actual_weak_coupling_sweep_gate_discharged',
                     'actual_spatial_summability_under_full_conditioning_proved',
                     'comparison_families_are_original_Wilson_vacua',
                     'continuum_state_survival_proved',
                     'positive_four_dimensional_continuum_mass_gap'):
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
            changed['controls']['heat_constant']['final_squared_coefficient'] = '1'
            result.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):
                y.check(cert,result,pin)
            changed = copy.deepcopy(cert)
            changed['model']['actual_weak_coupling_sweep_gate_discharged'] = True
            with self.assertRaises(ValueError):
                y.check(changed,result,pin)

    def test_upstream_source_drift_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        name = next(iter(manifest['upstream_sha256']))
        manifest['upstream_sha256'][name] = '0'*64
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'manifest'
            source.write_text(json.dumps(manifest))
            with patch.object(y,'SOURCES',source):
                with self.assertRaisesRegex(ValueError,'Changed upstream input'):
                    y.source_checks()

    def test_cli_is_read_only_and_rejects_optimized_runtime(self):
        inputs = (y.RESULT,y.PIN,y.SOURCES)
        before = [path.read_bytes() for path in inputs]
        checked = subprocess.run([sys.executable,'-B',str(y.HERE/'ym87_conditional_energy_transport.py'),'--check'],capture_output=True,text=True)
        self.assertEqual(checked.returncode,0,checked.stderr)
        self.assertIn('YM87 PASS',checked.stdout)
        self.assertEqual(before,[path.read_bytes() for path in inputs])
        optimized = subprocess.run([sys.executable,'-B','-O',str(y.HERE/'ym87_conditional_energy_transport.py'),'--check'],capture_output=True,text=True)
        self.assertNotEqual(optimized.returncode,0)
        self.assertIn('without optimization',optimized.stderr)


if __name__=='__main__':
    unittest.main()
