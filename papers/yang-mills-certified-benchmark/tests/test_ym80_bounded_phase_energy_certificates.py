import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from math import factorial
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym80_bounded_phase_energy_certificates as z


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


def gaussian_moment(degree):
    if degree % 2:
        return F(0)
    n = degree//2
    return F(factorial(2*n), 2**n*factorial(n))


def multiply(left, right):
    out = [F(0)]*(len(left)+len(right)-1)
    for i, first in enumerate(left):
        for j, second in enumerate(right):
            out[i+j] += first*second
    return out


def gaussian_expectation(poly):
    return sum((coefficient*gaussian_moment(degree) for degree, coefficient in enumerate(poly)), F(0))


class BoundedPhaseMathematicsTests(unittest.TestCase):
    def test_regularizer_has_certified_retention_and_third_moment_constants(self):
        row = z.regularizer_budget()
        self.assertEqual(row['eta'], F(7, 30))
        self.assertEqual(row['inverse_sine_budget'], F(183, 5)*F(7, 30)**2)
        self.assertEqual(row['beta'], F(1500, 4489))
        self.assertEqual(row['commutator_K'], 952)
        self.assertEqual(row['commutator_Z_upper'], 419280)
        self.assertEqual(row['squared_second_moment_constant_upper'], F(339994914968, 2401))
        self.assertEqual(row['second_moment_constant'], 12000)
        self.assertLess(row['squared_second_moment_constant_upper'], 12000**2)
        self.assertLessEqual(row['beta'], row['beta_sqrt_upper']**2)
        self.assertLessEqual(row['beta_sqrt_interval']['lower']**2, row['beta'])
        self.assertGreaterEqual(row['beta_sqrt_interval']['upper']**2, row['beta'])

    def test_fivefold_gain_uses_actual_directed_constants_not_rounded_decimals(self):
        row = z.regularizer_budget()
        legacy = F(1280, 1463)**2/(2**18*F(15, 16))
        self.assertGreater(row['asymptotic_unit_coefficient_lower'], 5*legacy)
        conservative = F(1500, 4489)**2/(12000*F(29, 50))
        self.assertEqual(conservative/legacy, F(12039575625, 2337530036))
        self.assertGreaterEqual(row['asymptotic_unit_coefficient_lower'], conservative)

    def test_retention_Jensen_bound_is_valid_for_a_nonconstant_two_point_variable(self):
        weights, values = (F(2, 5), F(3, 5)), (F(1, 9), F(7, 3))
        mean = sum(w*x for w, x in zip(weights, values))
        retained = sum(w/(1+x) for w, x in zip(weights, values))
        self.assertGreater(retained, 1/(1+mean))
        eta = z.regularizer_budget()['eta']
        beta = 1/(1+F(183, 5)*eta**2)
        self.assertEqual(beta, z.regularizer_budget()['beta'])
        # The old mixed-energy bound E[rho'] is not E[(rho')^2].
        with localcontext() as context:
            context.prec = 90
            mixed = sum(decimal(w)/(1+decimal(x)).sqrt() for w, x in zip(weights, values))
            self.assertGreater(mixed, decimal(retained))

    def test_phase_moment_bound_tracks_correlated_energy_not_independent_extrema(self):
        row = z.phase_moment_budget(3, 7, F(4, 9), time=2, beta_sqrt_upper=F(2, 3))
        self.assertEqual(row['phase_energy_lower'], F(4, 3))
        self.assertEqual(row['phase_energy_upper'], 3)
        self.assertEqual(row['phase_second_moment_upper'], 30)
        self.assertEqual(row['correlated_denominator'], F(227, 3))
        self.assertEqual(row['reflected_norm_squared_lower'], F(64, 681))
        self.assertEqual(row['centered_norm_squared_lower'], F(8, 81))

    def test_retention_quotient_monotonicity_against_independent_square_roots(self):
        with localcontext() as context:
            context.prec = 90
            e, r, time, beta = F(3, 2), F(7, 3), F(2, 5), F(1, 4)
            row = z.phase_moment_budget(e, r, beta, time=time, beta_sqrt_upper=F(1, 2))
            bounds = []
            for u in (beta, F(1, 3), F(3, 4), F(1)):
                ud, ed, rd, td = map(decimal, (u, e, r, time))
                bound = td**2*ed**2*ud**2/(1+ud*td*ed+ud*td**2*ed**2+ud.sqrt()*rd*td**2*ed)
                bounds.append(bound)
                self.assertLessEqual(decimal(row['reflected_norm_squared_lower']), bound)
            self.assertTrue(all(left < right for left, right in zip(bounds, bounds[1:])))

    def test_moment_budget_bounds_independent_two_atom_heat_measure(self):
        # m0=3/4,m1=1,m2=3/2, and m3 need not enter the phase budget here.
        energies, masses = (F(1), F(2)), (F(1, 2), F(1, 4))
        row = z.phase_moment_budget(1, F(1, 2), 1, time=F(3, 2))
        with localcontext() as context:
            context.prec = 90
            actual = sum(decimal(m)*(1-(-decimal(F(3, 2)*E)).exp())**2 for E, m in zip(energies, masses))
            self.assertLess(decimal(row['reflected_norm_squared_lower']), actual)
        self.assertLessEqual(row['centered_norm_squared_lower'], sum(masses))

    def test_phase_budget_zero_source_and_invalid_premises(self):
        empty = z.phase_moment_budget(0, 1, F(1, 2))
        self.assertEqual(empty['reflected_norm_squared_lower'], 0)
        self.assertEqual(empty['centered_norm_squared_lower'], 0)
        for args in ((-1, 1, F(1, 2)), (1, -1, F(1, 2)), (1, 1, 0), (1, 1, F(3, 2))):
            with self.assertRaises(ValueError): z.phase_moment_budget(*args)
        with self.assertRaises(ValueError): z.phase_moment_budget(1, 1, F(1, 2), time=0)
        with self.assertRaises(ValueError): z.phase_moment_budget(1, 1, F(1, 2), beta_sqrt_upper=F(1, 2))

    def test_actual_physical_budget_retains_clock_and_regulator_scaling(self):
        row = z.phase_energy_budget(25, F(1, 8), 2, phase_squared=F(3, 2), time=2, clock_scale=3)
        beta = F(1500, 4489)
        self.assertEqual(row['exact_raw_physical_energy'], 6)
        self.assertEqual(row['smooth_energy_lower'], 6*beta)
        self.assertEqual(row['smooth_energy_upper'], 6)
        self.assertEqual(row['smooth_second_moment_upper'], 12000*3*8*6)
        self.assertEqual(row['energy_budget'], 9)
        self.assertEqual(row['phase_energy_lower'], 9*beta)
        self.assertEqual(row['moment_scale'], 12000*3*8)
        self.assertLessEqual(F(7, 30)**2/5, 1)

    def test_fixed_source_certificate_is_order_a_and_not_a_positive_continuum_floor(self):
        rows = [z.phase_energy_budget(25, F(1, n), 1) for n in (8, 64, 4096)]
        coefficient = z.regularizer_budget()['asymptotic_unit_coefficient_lower']
        for row, n in zip(rows, (8, 64, 4096)):
            self.assertGreater(row['reflected_norm_squared_lower'], 0)
            self.assertLess(n*row['reflected_norm_squared_lower'], coefficient)
        self.assertTrue(all(left['reflected_norm_squared_lower'] > right['reflected_norm_squared_lower']
                            for left, right in zip(rows, rows[1:])))
        self.assertTrue(all(n*row['reflected_norm_squared_lower'] < m*later['reflected_norm_squared_lower']
                            for (n, row), (m, later) in zip(zip((8, 64), rows), zip((64, 4096), rows[1:]))))

    def test_actual_regularization_rejects_outside_analytic_domain(self):
        with self.assertRaises(ValueError): z.phase_energy_budget(24, F(1, 8), 1)
        with self.assertRaises(ValueError): z.phase_energy_budget(25, F(1, 8), 1, eta=3)
        with self.assertRaises(ValueError): z.phase_energy_budget(25, 0, 1)
        with self.assertRaises(ValueError): z.phase_energy_budget(25, F(1, 8), -1)
        with self.assertRaises(ValueError): z.phase_energy_budget(25, F(1, 8), 1, phase_squared=-1)

    def test_OU_intervals_enclose_independent_decimal_correlations(self):
        row = z.ou_phase_control(F(1, 3), raw_energy=2, time=F(3, 2), phase_parameter_squared=F(2, 5))
        with localcontext() as context:
            context.prec = 90
            first = (-Decimal('2.4')*(1-Decimal('-0.5').exp())).exp()
            second = (-Decimal('2.4')*(1-Decimal(-1).exp())).exp()
            truth = dict(phase_correlation_at_T=first, phase_correlation_at_2T=second,
                         phase_centered_mass=1-Decimal('-2.4').exp(), Weyl_reflected_norm=1-2*first+second)
            for key, value in truth.items():
                self.assertLessEqual(decimal(row[key]['lower']), value)
                self.assertGreaterEqual(decimal(row[key]['upper']), value)

    def test_OU_Poisson_moments_validate_phase_chain_constant(self):
        row = z.ou_phase_control(F(5, 2), raw_energy=3, phase_parameter_squared=F(1, 4))
        self.assertEqual(row['phase_energy'], F(3, 4))
        self.assertEqual(row['Poisson_mean'], F(3, 10))
        self.assertEqual(row['phase_second_moment'], F(39, 16))

    def test_phase_chain_matches_exact_quadratic_Hermite_moments(self):
        # Standard Gaussian, L=-d²+x d. G=He2=x²-1, ΓG=4x².
        G, LG, Gamma = [F(-1), 0, 1], [F(-2), 0, 2], [0, 0, F(4)]
        raw_energy = gaussian_expectation(Gamma)
        raw_second = gaussian_expectation(multiply(LG, LG))
        Gamma_second = gaussian_expectation(multiply(Gamma, Gamma))
        self.assertEqual((raw_energy, raw_second, Gamma_second), (4, 8, 48))
        phase = z.phase_chain_moments(F(9, 64), raw_energy, raw_second, Gamma_second)
        self.assertEqual(phase['phase_energy'], F(9, 16))
        self.assertEqual(phase['phase_second_moment'], F(531, 256))
        self.assertTrue(phase['identity_fixture_only'])
        self.assertFalse(phase['pointwise_Gamma_bound_inferred'])
        # μΓ² is48, not(μΓ)²=16: the nonlinear-chain term cannot be factorized.
        self.assertNotEqual(phase['phase_second_moment'], F(9, 64)*8+F(9, 64)**2*16)

    def test_phase_chain_matches_cubic_Hermite_without_constant_Gamma_assumption(self):
        # G=He3=x³-3x, LG=3G, ΓG=9(x²-1)².
        LG = [0, F(-9), 0, F(3)]
        Gamma = [F(9), 0, F(-18), 0, F(9)]
        raw_energy = gaussian_expectation(Gamma)
        raw_second = gaussian_expectation(multiply(LG, LG))
        Gamma_second = gaussian_expectation(multiply(Gamma, Gamma))
        self.assertEqual((raw_energy, raw_second, Gamma_second), (18, 54, 4860))
        row = z.phase_chain_moments(F(1, 9), raw_energy, raw_second, Gamma_second)
        self.assertEqual(row['phase_energy'], 2)
        self.assertEqual(row['phase_second_moment'], 66)

    def test_OU_sharpness_saturates_raw_third_moment_and_phase_energy_retention(self):
        row = z.ou_sharpness_control(8, raw_energy=3, energy_budget=2, beta=F(1, 4))
        self.assertLessEqual(row['rate_interval']['lower'], 16)
        self.assertGreaterEqual(row['rate_interval']['upper'], 16)
        self.assertEqual(row['raw_energy'], F(3, 4))
        self.assertEqual(row['raw_third_moment'], 192)
        self.assertEqual(row['phase_energy'], F(1, 2))
        self.assertLessEqual(row['raw_second_moment_interval']['lower'], 12)
        self.assertGreaterEqual(row['raw_second_moment_interval']['upper'], 12)
        self.assertLessEqual(row['asymptotic_mass_scale_interval']['lower'], F(1, 32))
        self.assertGreaterEqual(row['asymptotic_mass_scale_interval']['upper'], F(1, 32))
        with localcontext() as context:
            context.prec = 90
            exact_mass = 1-(-Decimal(1)/32).exp()
            self.assertLessEqual(decimal(row['phase_centered_mass']['lower']), exact_mass)
            self.assertGreaterEqual(decimal(row['phase_centered_mass']['upper']), exact_mass)

    def test_bounded_OU_phases_can_lose_all_survival_despite_fixed_energy(self):
        for n in (8, 64, 512):
            row = z.ou_sharpness_control(n)
            self.assertEqual(row['phase_energy'], 1)
            self.assertEqual(row['raw_third_moment'], n*n)
            self.assertLessEqual(row['phase_centered_mass']['upper'], F(1, n)+F(1, 10**17))
            self.assertLessEqual(row['Weyl_reflected_norm']['upper'], F(1, n)+F(1, 10**17))
            certified = z.phase_moment_budget(1, n, 1)['reflected_norm_squared_lower']
            self.assertLessEqual(certified, row['Weyl_reflected_norm']['lower'])

    def test_high_amplitude_survival_spends_unbounded_energy_and_loses_delayed_correlation(self):
        rows = [z.high_amplitude_control(n) for n in (8, 64)]
        for n, row in zip((8, 64), rows):
            self.assertEqual(row['phase_energy'], n)
            self.assertTrue(row['violates_fixed_energy_budget'])
            self.assertGreater(row['Weyl_reflected_norm']['lower'], F(3, 5))
        self.assertLess(rows[-1]['phase_centered_correlation_at_T']['upper'],
                        rows[0]['phase_centered_correlation_at_T']['lower'])
        self.assertLess(rows[-1]['phase_centered_correlation_at_T']['upper'], F(1, 10**16))

    def test_fixed_energy_optimization_agrees_with_explicit_phase_choice(self):
        row = z.optimized_phase_budget(25, F(1, 8), 2, physical_energy_budget=3, clock_scale=2)
        direct = z.phase_energy_budget(25, F(1, 8), 2, phase_squared=F(3, 4), clock_scale=2)
        self.assertEqual(row['energy_budget'], 3)
        for key in ('phase_energy_lower', 'phase_energy_upper', 'reflected_norm_squared_lower'):
            self.assertEqual(row[key], direct[key])

    def test_UV_packet_bounds_independent_multilevel_spectral_measures(self):
        fixtures = (((F(1, 2), 1, 2), (F(1, 4), F(1, 4), F(3, 8)), F(1)),
                    ((1, 32), (F(1, 2), F(1, 2048)), F(1, 2)))
        for relative_energies, weights, beta in fixtures:
            for scale in (F(1, 3), 2, 9):
                row = z.ultraviolet_packet_budget(scale, beta)
                energies = tuple(scale*energy for energy in relative_energies)
                mass = sum(weights)
                first = sum(weight*energy for weight, energy in zip(weights, energies))
                second = sum(weight*energy**2 for weight, energy in zip(weights, energies))
                self.assertLessEqual(mass, 1)
                self.assertGreaterEqual(first, beta*scale)
                self.assertLessEqual(second, 2*scale**2)
                lower, upper = row['spectral_window']['lower'], row['spectral_window']['upper']
                above_lower = sum(weight for weight, energy in zip(weights, energies) if energy >= lower)
                above_upper = sum(weight for weight, energy in zip(weights, energies) if energy > upper)
                in_band = sum(weight for weight, energy in zip(weights, energies) if lower <= energy <= upper)
                self.assertGreaterEqual(above_lower, row['mass_above_lower_endpoint_lower'])
                self.assertLessEqual(above_upper, row['mass_above_upper_endpoint_upper'])
                self.assertGreaterEqual(in_band, row['spectral_mass_lower'])
                self.assertEqual(lower/scale, beta/2)
                self.assertEqual(upper/scale, 8/beta)
                self.assertFalse(row['fixed_physical_energy_budget'])
                self.assertFalse(row['packet_is_fixed_physical_energy_window'])

    def test_UV_packet_rejects_impossible_retention_or_nonpositive_scale(self):
        for scale, beta in ((0, F(1, 2)), (-1, F(1, 2)), (1, 0), (1, -1), (1, F(3, 2))):
            with self.assertRaises(ValueError): z.ultraviolet_packet_budget(scale, beta)

    def test_separated_time_budget_agrees_with_conditional_spectral_Jensen(self):
        # A positive-energy measure of mass1/2 atE=2 has energy1;
        # D1=(1/2)(1-e^-2)^2>1/3, allowingδ=1/3.
        row = z.separated_time_budget(1, F(1, 2), F(1, 3), 1)
        with localcontext() as context:
            context.prec = 90
            claimed = Decimal(1)/3*Decimal(-3).exp()
            actual = Decimal(1)/2*Decimal(-2).exp()*(1-Decimal(-2).exp())**2
            value = decimal(row['separated_reflected_norm_squared_lower'])
            self.assertLessEqual(value, claimed)
            self.assertLess(value, actual)
            self.assertGreater(value, 0)
        zero_delay = z.separated_time_budget(1, 0, F(1, 3), 1)
        self.assertEqual(zero_delay['separated_reflected_norm_squared_lower'], F(1, 3))

    def test_separated_time_bound_does_not_manufacture_survival_or_accept_inconsistent_energy(self):
        row = z.separated_time_budget(1, 1, 0, 1)
        self.assertEqual(row['separated_reflected_norm_squared_lower'], 0)
        with self.assertRaises(ValueError): z.separated_time_budget(1, 1, 2, 1)
        with self.assertRaises(ValueError): z.separated_time_budget(1, -1, F(1, 2), 1)
        with self.assertRaises(ValueError): z.separated_time_budget(0, 1, F(1, 2), 1)


class BoundedPhaseCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_manifest_and_parent_bind_the_sources(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 295)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        self.assertEqual(self.cert['parent_result_sha256'], z.PARENT_HASH)
        self.assertEqual(z.previous.PIN.read_text().strip(), z.PARENT_HASH)

    def test_upstream_parent_and_manifest_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'PARENT_HASH', '0'*64):
            with self.assertRaises(ValueError): z.source_checks()
        with patch.object(z, 'BASE_TREE', '0'*40):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_continuum_and_UV_repair_overclaims_are_rejected(self):
        for field in ('uniform_positive_continuum_survival', 'normalization_fixes_UV_escape',
                      'raw_linear_source_order_a_bound_proved', 'positive_continuum_mass_gap'):
            wrong = copy.deepcopy(self.cert)
            wrong['model'][field] = True
            with tempfile.TemporaryDirectory() as directory:
                result, pin = Path(directory)/'result', Path(directory)/'pin'
                result.write_text(json.dumps(wrong))
                pin.write_text(z.p.canonical_sha(wrong))
                with self.assertRaises(ValueError): z.check(wrong, result, pin)
                with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_certificate_survives_JSON_roundtrip_without_binary_floats(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                self.assertTrue(all(isinstance(key, str) for key in value))
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        self.assertEqual(json.loads(json.dumps(self.cert)), self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_model_does_not_turn_a_finite_bound_into_continuum_survival(self):
        self.assertTrue(self.cert['model']['actual_finite_phase_certificate'])
        self.assertTrue(self.cert['model']['actual_new_bounded_smooth_phase_has_order_a_certificate'])
        self.assertFalse(self.cert['model']['raw_linear_source_order_a_bound_proved'])
        self.assertFalse(self.cert['model']['uniform_positive_continuum_survival'])
        self.assertFalse(self.cert['model']['normalization_fixes_UV_escape'])


if __name__ == '__main__':
    unittest.main()
