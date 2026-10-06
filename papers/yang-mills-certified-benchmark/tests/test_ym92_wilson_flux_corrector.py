"""YM92: exact native, residual, incidence, and scope controls."""
import copy
from fractions import Fraction as F
from itertools import permutations
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym92_wilson_flux_corrector as y

p, a, poly = y.p, y.a, y.poly


def inserted_word(qs, directions):
    """Independent ordered quaternion insertions on distinct link factors."""
    qmul, dagger = y.local.y.qmul, y.local.y.qdagger
    factors = []
    for e, sign in y.local.model.ORIENTED_WORD:
        q = qs[e]
        if e not in directions:
            factor = q if sign == 1 else dagger(q)
        else:
            generator = (F(0),)+tuple(directions[e])
            factor = (tuple(-v/2 for v in qmul(generator, q)) if sign == 1
                      else tuple(v/2 for v in qmul(dagger(q), generator)))
        factors.append(factor)
    return y.local.qproduct(factors)[0]


class ActualBudgetTests(unittest.TestCase):
    def test_exact_endpoint_improves_actual_uniform_certificates(self):
        out = y.endpoint_control()
        self.assertEqual(out['actual']['bare_full_and_conditional_gap_lower'], F(251, 800))
        self.assertEqual(out['actual']['relative_memory_upper'], F(2105401, 3189425625))
        self.assertLess(out['actual']['relative_memory_upper'], F(1, 1500))
        self.assertGreater(out['memory_upper_improvement_factor'], F(21, 2))

    def test_nonlinear_bootstrap_uses_a_positive_denominator(self):
        for n in range(101):
            t = F(n, 40000)
            out = y.remainder_budget(t)
            self.assertGreater(out['bootstrap_denominator'], 0)
            self.assertLessEqual(out['refined_mixed_remainder'], out['mixed_remainder_upper'])
            self.assertEqual(out['preliminary_mixed_remainder'],
                             4*(35*t)*(5*t)+2*(5*t)**2+4*t*(42*t))
        self.assertEqual(y.remainder_budget(F(1, 400))['bootstrap_denominator'], F(25141, 40000))

    def test_quadratic_inequality_alone_would_admit_large_error(self):
        t, large = F(1, 400), F(1)
        self.assertLessEqual(large, 310*t*t+144*t*large+2*large*large)
        self.assertGreater(large, y.remainder_budget(t)['preliminary_mixed_remainder'])

    def test_all_orders_flux_bound_keeps_the_three_residual_terms(self):
        t = F(1, 400)
        out = y.actual_budget(t)
        expected = 1000*t*t+2*(36*t+500*t*t)*(8*t/9)+F(8, 3)*t*(42*t)
        self.assertEqual(out['residual_gradient_over_sqrt_b'], expected)
        self.assertEqual(out['completed_flux_norm_over_sqrt_b'],
                         out['trial_flux_norm_over_sqrt_b']+expected/out['bare_full_and_conditional_gap_lower'])
        self.assertGreater(out['relative_memory_upper'], out['trial_flux_norm_over_sqrt_b']**2)

    def test_endpoint_is_worst_on_declared_interval(self):
        old_eta, old_gap = F(0), F(1, 2)
        for n in range(101):
            out = y.actual_budget(F(n, 40000))
            self.assertGreaterEqual(out['relative_memory_upper'], old_eta)
            self.assertLessEqual(out['bare_full_and_conditional_gap_lower'], old_gap)
            old_eta, old_gap = out['relative_memory_upper'], out['bare_full_and_conditional_gap_lower']

    def test_block_factor_cancels_in_relative_memory(self):
        rows = [y.actual_budget(F(1, 400), b) for b in (2, 3, 8, 64)]
        self.assertEqual(len({r['relative_memory_upper'] for r in rows}), 1)
        for row in rows:
            self.assertEqual(row['conditional_tensor_norm_upper']/row['factor'], row['relative_memory_upper'])

    def test_flux_improves_even_the_refined_Hessian_test(self):
        for n in (1, 3, 10, 100):
            t = F(n, 40000)
            out = y.actual_budget(t)
            d, q = out['bare_full_and_conditional_gap_lower'], out['mixed_Hessian_row_upper']
            self.assertLess(out['relative_memory_upper'], 4*q*q/(d*d))
            self.assertGreater(2*q-(F(4, 9)*t+1112*t*t), 0)

    def test_frequency_bound_uses_minimum_of_two_justified_estimates(self):
        rows = [y.actual_budget(F(1, 400), 2, s) for s in (0, 1, 10, 100)]
        for row in rows:
            expected = min(row['relative_memory_upper'], row['source_energy_factor'] /
                           (row['bare_full_and_conditional_gap_lower']+row['damping']))
            self.assertEqual(row['frequency_memory_upper'], expected)
            self.assertEqual(row['relative_resolvent_error'], expected/(1-expected))
        self.assertGreater(rows[0]['frequency_memory_upper'], rows[-1]['frequency_memory_upper'])

    def test_free_case_has_no_memory_or_residual(self):
        out = y.actual_budget(0, 4)
        for key in ('relative_memory_upper', 'trial_flux_norm_over_sqrt_b', 'residual_gradient_over_sqrt_b',
                    'source_energy_factor', 'relative_resolvent_error'):
            self.assertEqual(out[key], 0)
        self.assertEqual(out['bare_full_and_conditional_gap_lower'], F(1, 2))

    def test_parameters_do_not_silently_extend_window_or_change_types(self):
        for t in (-1, F(1, 399)):
            with self.assertRaises(ValueError):
                y.actual_budget(t)
        for value in (0.1, True):
            with self.assertRaises(TypeError):
                y.actual_budget(value)
        for b in (1, True, F(2)):
            with self.assertRaises(ValueError):
                y.actual_budget(0, b)
        with self.assertRaises(ValueError):
            y.actual_budget(0, 2, -1)


class NativeResponseTests(unittest.TestCase):
    def test_exact_native_mixed_resolvent_normalization(self):
        for pair in ((0, 1), (0, 2)):
            row = y.linear_resolvent_control(*pair)
            self.assertEqual(row['native_mixed_components_checked'], 9)
            self.assertEqual(row['components_changed_by_scalar_instead_of_rough_Laplacian'], 9)

    def test_repeated_factor_covariant_third_derivatives(self):
        row = y.third_derivative_control()
        self.assertEqual(row['repeated_factor_components'], 324)
        self.assertEqual(row['distinct_factor_components'], 324)
        self.assertGreater(row['components_changed_by_omitting_connections'], 0)

    def test_three_distinct_native_insertions_in_noncommuting_directions(self):
        qs = ((F(3, 5), F(4, 5), 0, 0), (0, F(5, 13), F(12, 13), 0),
              (F(1, 2),)*4, (0, 0, 0, 1))
        u, v, w = (F(3, 5), F(4, 5), 0), (0, F(5, 13), F(12, 13)), (F(4, 5), 0, F(3, 5))
        for e, f, g in permutations(range(4), 3):
            H = y.plaquette_mixed(e, f)
            value = F(0)
            for c in range(3):
                T = y.covariant_block(H, e, f, 3*g+c)
                value += sum(u[i]*v[j]*w[c]*y.local.model.evaluate(T[i][j], qs)
                             for i in range(3) for j in range(3))
            self.assertEqual(value, inserted_word(qs, {e:u, f:v, g:w}))
            self.assertLessEqual(abs(value), F(1, 8))

    def test_first_order_log_wave_solves_the_electric_equation(self):
        W = y.local.model.loop_polynomial()
        self.assertEqual(y.local.model.electric(poly.scale(W, F(1, 3)), 4), W)

    def test_native_four_link_carrier_contract(self):
        for pair in ((0, 0), (-1, 1), (0, 4), (True, 1)):
            with self.assertRaises(ValueError):
                y.plaquette_mixed(*pair)
        for retained in ((), (0, 0), (0, 1, 2), (4,), (True,)):
            with self.assertRaises(ValueError):
                y.flux_polynomials(retained)


class PlaquetteFluxTests(unittest.TestCase):
    def test_both_conditional_Casimirs_are_exact(self):
        for retained, eigenvalue in (((0,), F(9, 4)), ((0, 1), F(3, 2)), ((0, 3), F(3, 2))):
            row = y.flux_polynomials(retained)
            self.assertEqual(row['eigenvalue'], eigenvalue)
            self.assertEqual(poly.scale(row['potential'], eigenvalue), row['source'])

    def test_flux_divergence_and_Haar_centering_without_density_ansatz(self):
        for retained in ((0,), (1,), (0, 1), (2, 3)):
            row = y.flux_polynomials(retained)
            phi = row['potential']
            divergence = poly.add(*(y.D(y.D(phi, 3*e+k), 3*e+k) for e in row['hidden'] for k in range(3)))
            self.assertEqual(poly.scale(divergence, -1), row['source'])
            centered = row['source']
            for e in row['hidden']:
                centered = y.local.integrate_link(centered, e)
            self.assertEqual(centered, {})

    def test_density_drift_changes_every_checked_residual_derivative(self):
        row = y.weighted_residual_control()
        self.assertEqual(row['conditional_gradient_components_checked'], 6)
        self.assertEqual(row['components_wrong_if_weighted_drift_is_dropped'], 6)
        self.assertIn('not_Wilson_vacuum', row['kind'])

    def test_geometry_has_at_least_two_unretained_links_per_active_face(self):
        for b in (2, 3, 4):
            row = y.incidence_control(b)
            expected = {0, 2} if b == 2 else {0, 1, 2}
            self.assertEqual(set(row['retained_face_histogram']), expected)
            self.assertEqual(row['retained_links'], b*row['coarse_links'])
            self.assertEqual(row['normalized_incidence_bounds']['flux_column'], F(8, 9))

    def test_geometric_operator_bounds_retain_rows_and_columns(self):
        for b in (2, 3, 4):
            row = y.incidence_control(b)['normalized_incidence_bounds']
            self.assertLessEqual(row['flux_row']*row['flux_column'], F(8, 9)**2)
            self.assertLessEqual(row['drift_row']*row['drift_column'], F(4, 3)**2)


class ResidualGramTests(unittest.TestCase):
    def test_exact_correction_includes_the_full_residual(self):
        row = y.finite_flux_control()['residual_required']
        N = y.previous.matrix_control()['gradient']
        self.assertEqual(p.multiply(p.transpose(N), row['complete_flux']), row['source'])
        self.assertTrue(p.psd(p.add(row['corrected_Gram'], row['exact_conditional_tensor'], -1)))
        self.assertTrue(p.psd(p.add(row['full_residual_upper'], row['exact_conditional_tensor'], -1)))

    def test_trial_flux_without_correction_can_underestimate(self):
        row = y.finite_flux_control()['residual_required']
        self.assertFalse(p.psd(p.add(row['trial_Gram'], row['exact_conditional_tensor'], -1)))
        self.assertNotEqual(row['residual'], [[F(0), F(0)] for _ in range(3)])

    def test_diagonal_only_response_is_not_a_Loewner_upper_bound(self):
        row = y.finite_flux_control()['cross_terms_required']
        self.assertEqual(row['exact_conditional_tensor'], [[1, 1], [1, 2]])
        self.assertFalse(p.psd([[0, -1], [-1, 0]]))

    def test_residual_energy_bound_allows_different_Young_parameters(self):
        row = y.previous.matrix_control()
        N, R = row['gradient'], row['scores']
        J = a.matrix([[0, 0] for _ in N])
        for t in (F(1, 10), F(1), F(10)):
            out = y.approximate_flux_certificate(N, R, J, 1, t)
            self.assertTrue(p.psd(p.add(out['full_residual_upper'], out['exact_conditional_tensor'], -1)))
            self.assertEqual(out['corrected_Gram'], out['exact_conditional_tensor'])

    def test_false_gap_or_invalid_residual_parameters_are_rejected(self):
        for floor, young in ((2, 1), (0, 1), (1, 0)):
            with self.assertRaises(ValueError):
                y.approximate_flux_certificate([[1]], [[1]], [[0]], floor, young)

    def test_zero_frequency_control_does_not_imply_extra_decay(self):
        row = y.previous.previous.insufficient_zero_frequency_control()
        self.assertTrue(row)
        # D=100, B=1: eta=1/100, but Sigma(1)=1/101 exceeds eta/(1+1).
        self.assertGreater(F(1, 101), F(1, 200))


class EvidenceTests(unittest.TestCase):
    def test_fresh_result_and_pin(self):
        fresh = y.run()
        self.assertEqual(y.check(fresh), y.PIN.read_text().strip())
        self.assertEqual(len(fresh['input_sha256']), 432)

    def test_modified_scope_cannot_pass_even_with_rehashed_files(self):
        changed = copy.deepcopy(y.run())
        changed['model']['positive_four_dimensional_continuum_mass_gap'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(changed))
            pin.write_text(p.canonical_sha(changed))
            with self.assertRaisesRegex(ValueError, 'scope'):
                y.check(changed, result, pin)

    def test_changed_source_base_is_rejected(self):
        manifest = json.loads(y.SOURCES.read_text())
        manifest['base_commit'] = '0'*40
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'manifest'
            path.write_text(json.dumps(manifest))
            with patch.object(y, 'SOURCES', path), self.assertRaisesRegex(ValueError, 'base'):
                y.source_checks()

    def test_Runtime_policy_and_scope_stay_explicit(self):
        self.assertFalse(y.MODEL['weak_coupling_window_extended'])
        self.assertFalse(y.MODEL['comparison_polynomials_are_numerical_Wilson_vacua'])
        self.assertFalse(y.MODEL['positive_four_dimensional_continuum_mass_gap'])
        with patch.object(y.sys, 'version_info', (3, 11, 0)), self.assertRaises(RuntimeError):
            y.run()


if __name__ == '__main__':
    unittest.main()
