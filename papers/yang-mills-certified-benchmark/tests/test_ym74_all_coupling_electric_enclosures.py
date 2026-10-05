import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym74_all_coupling_electric_enclosures as z


class ElectricEnclosureMathematicsTests(unittest.TestCase):
    def test_directed_square_roots_include_exact_and_very_small_values(self):
        for value in (F(0), F(1, 4), F(2), F(10**40+1), F(1, 10**100)):
            interval = z.sqrt_interval(value)
            self.assertLessEqual(interval['lower']**2, value)
            self.assertGreaterEqual(interval['upper']**2, value)
            self.assertLessEqual(interval['upper']-interval['lower'], F(1, z.GRID))
        self.assertEqual(z.sqrt_interval(F(1, 4)), dict(lower=F(1, 2), upper=F(1, 2)))

    def test_integer_ceiling_arithmetic_has_no_float_roundoff(self):
        for value in (F(0), F(1, 100), F(4), F(400000000000000000001, 10**20), F(10**40)):
            ceil = z.ceil_sqrt(value)
            self.assertGreaterEqual(ceil*ceil, value)
            if ceil:
                self.assertLess((ceil-1)**2, value)
        self.assertEqual(z.ceil_fraction(F(-5, 3)), -1)
        self.assertEqual(z.ceil_fraction(F(5, 3)), 2)

    def test_exact_domains_reject_binary_floats_booleans_and_negative_budgets(self):
        for bad in (-1, 0.5, True):
            for function in (z.sqrt_interval, z.ceil_sqrt, z.vacuum_budget, z.finite_shell):
                with self.assertRaises(ValueError): function(bad)
        for bad in (-1, F(1, 2), True):
            with self.assertRaises(ValueError): z.product_trial(bad)
        with self.assertRaises(ValueError): z.vacuum_budget(1, 0)
        with self.assertRaises(ValueError): z.representation_rank_bound(3, 0)

    def test_product_trial_matches_independent_native_Haar_integration(self):
        rows = z.trial_controls()['independent_native_integrations']
        self.assertEqual(len(rows), 8)
        self.assertEqual(rows[1]['independent_Haar_norm'], F(5, 4))
        self.assertEqual(rows[1]['independent_Haar_mean'], F(2, 5))
        self.assertEqual(rows[1]['independent_Haar_electric'], F(3, 20))
        for row in rows:
            native = z.product_trial(row['index'])
            self.assertEqual(native['link_scalar_mean'], row['independent_Haar_mean'])
            self.assertEqual(native['link_electric_energy'], row['independent_Haar_electric'])

    def test_beta_density_has_correct_mean_and_electric_inverse_moment(self):
        for n in (1, 2, 7, 100):
            row = z.product_trial(n)
            alpha, beta = row['beta_alpha'], row['beta_beta']
            self.assertEqual(row['link_scalar_mean'], 2*alpha/(alpha+beta)-1)
            self.assertEqual(row['link_electric_energy'], F(n*n, 4)*beta/(alpha-1))
            self.assertEqual(row['plaquette_mean'], row['link_scalar_mean']**4)

    def test_free_constant_trial_and_actual_zero_coupling_bounds(self):
        trial = z.product_trial(0)
        self.assertEqual(trial['link_electric_energy'], 0)
        self.assertEqual(trial['plaquette_mean'], 0)
        row = z.vacuum_budget(0, 192)
        self.assertEqual(row['ground_energy_lower'], 0)
        self.assertEqual(row['ground_energy_upper'], 0)
        self.assertEqual(row['actual_total_electric_energy_upper'], 0)
        self.assertEqual(row['actual_plaquette_deficiency_upper'], 1)

    def test_actual_vacuum_deficit_splits_into_two_nonnegative_forms(self):
        for theta in (F(1, 8000), F(1), F(64), F(10**6)):
            row = z.vacuum_budget(theta, 192)
            deficit = row['actual_shifted_vacuum_deficit_upper']
            self.assertEqual(row['ground_energy_upper']+192*theta, deficit)
            self.assertEqual(row['actual_total_electric_energy_upper'], deficit)
            self.assertEqual(row['actual_plaquette_deficiency_upper'], deficit/(192*theta))
            self.assertGreaterEqual(deficit, 0)
            self.assertLessEqual(deficit, 192*theta)

    def test_optimized_trial_density_obeys_sqrt_theta_bound(self):
        for theta in (F(1, 10**20), F(1, 32), F(1), F(16), F(10**6)):
            row = z.vacuum_budget(theta)
            n = row['trial']['index']
            self.assertGreaterEqual(n*n, 32*theta)
            self.assertLess((n-1)**2, 32*theta)
            density = row['trial']['link_electric_energy']+theta*row['trial']['plaquette_deficiency']
            self.assertLessEqual(density, row['optimized_density_asymptotic_upper'])
            self.assertLessEqual(row['trial']['link_electric_energy'], F(3*n, 16))
            self.assertLessEqual(row['trial']['plaquette_deficiency'], F(6, n))

    def test_weak_coupling_flatness_is_local_and_volume_uniform(self):
        rows = z.trial_controls()['actual_vacuum_weak_coupling_budgets']
        self.assertEqual([row['theta'] for row in rows], [1, 16, 256, 4096, 65536, 10**6])
        self.assertLess(rows[-1]['actual_plaquette_deficiency_upper'], F(22, 10000))
        for row in rows:
            big = z.vacuum_budget(row['theta'], 3*100**3)
            self.assertEqual(big['actual_plaquette_deficiency_upper'], row['actual_plaquette_deficiency_upper'])
            self.assertEqual(row['actual_holonomy_Frobenius_squared_upper'],
                             4*row['actual_plaquette_deficiency_upper'])

    def test_hidden_cutoff_uses_improved_actual_ground_energy(self):
        row = z.electric_cutoff_budget(256, 192, 10**5, physical_band=3, clock=F(1, 2))
        self.assertGreater(row['hidden_bare_floor'], F(10**5)-256*192)
        self.assertEqual(row['hidden_physical_floor'], row['hidden_bare_floor']/2)
        self.assertTrue(row['positive_hidden_floor_certified'])
        self.assertEqual(row['entire_low_energy_subspace_leakage_squared_upper'],
                         (row['vacuum_deficit_upper']+6)/10**5)

    def test_coverage_saturates_safely_and_does_not_assert_an_excitation_exists(self):
        row = z.electric_cutoff_budget(256, 192, 1, 100)
        self.assertEqual(row['entire_low_energy_subspace_leakage_squared_upper'], 1)
        self.assertFalse(row['positive_hidden_floor_certified'])
        self.assertFalse(row['nonvacuum_low_energy_states_asserted'])
        vacuum = z.electric_cutoff_budget(0, 192, 3, 0)
        self.assertEqual(vacuum['vacuum_leakage_squared_upper'], 0)

    def test_whole_band_coverage_against_an_independent_finite_spectral_example(self):
        # Diagonal H has known physical band [0,3]; a direct projector calculation
        # checks the Markov form argument, including a hidden component in that band.
        electric = (F(0), F(2), F(5), F(9))
        shifted = (F(0), F(1), F(3), F(7))
        deficit, band, cutoff = F(2), F(3), F(4)
        for e, h in zip(electric, shifted):
            self.assertLessEqual(e, h+deficit)
        next_electric = F(5)
        exact_leakage_norm_squared = max(F(e > cutoff and h <= band) for e, h in zip(electric, shifted))
        self.assertEqual(exact_leakage_norm_squared, 1)
        self.assertLessEqual(exact_leakage_norm_squared, (deficit+band)/next_electric)
        # Removing the M/clock term would make a false bound of 2/5.
        self.assertGreater(exact_leakage_norm_squared, deficit/next_electric)

    def test_combined_cutoff_gives_hidden_divergence_and_uniform_band_capture(self):
        for growth in (F(1), F(4), F(16), F(64)):
            row = z.combined_cutoff(256, 192, F(1, 100), growth,
                                    reference_energy=2, physical_band=3)
            self.assertGreaterEqual(row['actual_budget']['hidden_physical_floor'], 2*growth)
            self.assertLessEqual(row['actual_budget']['entire_low_energy_subspace_leakage_squared_upper'],
                                 min(F(1), F(3, 2)/growth))
            self.assertIsNone(row['physical_gap_total_width_upper'])

    def test_precision_cutoff_controls_total_gap_width_not_one_endpoint(self):
        theta, count, clock, epsilon = F(16), 192, F(3, 5), F(1, 1000)
        row = z.combined_cutoff(theta, count, clock, 4, precision=epsilon)
        independent = z.ritz_width_budget(theta*count, row['cutoff'], clock)
        self.assertEqual(row['physical_gap_total_width_upper'], independent['physical_total_width_upper'])
        self.assertLessEqual(independent['physical_total_width_upper'], epsilon)
        self.assertEqual(row['precision_cutoff'], 4+z.ceil_fraction(2*theta*count+2*clock*(theta*count)**2/epsilon))

    def test_cutoff_domains_and_separation_premises_are_explicit(self):
        for kwargs in (dict(cutoff=0), dict(cutoff=4, next_electric_lower=3),
                       dict(cutoff=4, clock=0), dict(cutoff=4, physical_band=-1)):
            with self.assertRaises(ValueError): z.electric_cutoff_budget(1, 192, **kwargs)
        with self.assertRaises(ValueError): z.combined_cutoff(1, 192, 1, F(1, 2))
        with self.assertRaises(ValueError): z.combined_cutoff(1, 192, 1, 2, precision=0)
        with self.assertRaises(ValueError): z.ritz_width_budget(2, 7)

    def test_finite_plaquette_coupling_shell_covers_all_four_spin_shifts(self):
        spins = tuple(F(index, 2) for index in range(6))
        for cutoff in (F(3), F(8), F(17, 2)):
            upper = z.finite_shell(cutoff)['exact_coupling_shell_upper']
            count = 0
            for four in product(spins, repeat=4):
                initial = sum(spin*(spin+1) for spin in four)
                if initial <= cutoff:
                    maximal = sum((spin+F(1, 2))*(spin+F(3, 2)) for spin in four)
                    self.assertLessEqual(maximal, upper)
                    count += 1
            self.assertGreater(count, 1)

    def test_finite_rank_bound_uses_Peter_Weyl_multiplicities(self):
        for cutoff in (F(0), F(3, 4), F(2), F(3), F(100, 7)):
            row = z.representation_rank_bound(cutoff, 4)
            n = row['maximum_twice_edge_spin']
            self.assertLessEqual(F(n*(n+2), 4), cutoff)
            self.assertGreater(F((n+1)*(n+3), 4), cutoff)
            dimension = sum((k+1)**2 for k in range(n+1))
            self.assertEqual(row['unrestricted_single_edge_dimension'], dimension)
            self.assertEqual(row['physical_retained_rank_upper'], dimension**4)

    def test_Ritz_root_perfect_square_and_zero_coupling_limits(self):
        self.assertEqual(z.ritz_root_interval(0, 3, 4), dict(lower=F(-1), upper=F(-1)))
        self.assertEqual(z.ritz_root_interval(3, 5, 0), dict(lower=F(3), upper=F(3)))
        self.assertEqual(z.ritz_root_interval(5, 2, 0), dict(lower=F(2), upper=F(2)))
        with self.assertRaises(ValueError): z.ritz_root_interval(1, 2, -1)

    def test_Ritz_root_outward_rounding_matches_independent_decimal_arithmetic(self):
        with localcontext() as context:
            context.prec = 100
            for c, d, q in ((F(-2, 3), F(7, 2), F(5, 13)),
                            (F(10), F(-3), F(7)), (F(1, 10**20), F(2), F(1, 10**30))):
                cv, dv, qv = (Decimal(value.numerator)/Decimal(value.denominator) for value in (c, d, q))
                expected = (cv+dv-((cv-dv)**2+4*qv).sqrt())/2
                row = z.ritz_root_interval(c, d, q)
                low = Decimal(row['lower'].numerator)/Decimal(row['lower'].denominator)
                high = Decimal(row['upper'].numerator)/Decimal(row['upper'].denominator)
                self.assertLessEqual(low, expected)
                self.assertGreaterEqual(high, expected)

    def test_uncertain_Ritz_inputs_use_conservative_endpoint_directions(self):
        base = z.ritz_gap_interval((F(-1, 10), F(-1, 10)), (3, 3), 5, F(1, 4))
        wider = z.ritz_gap_interval((F(-11, 100), F(-9, 100)),
                                    (F(299, 100), F(301, 100)), F(49, 10), F(26, 100))
        self.assertLessEqual(wider['lower'], base['lower'])
        self.assertGreaterEqual(wider['upper'], base['upper'])
        for ground, excited in (((1, 0), (3, 3)), ((0, 0), (3, 2)), ((5, 5), (3, 3))):
            with self.assertRaises(ValueError): z.ritz_gap_interval(ground, excited, 4, 1)
        with self.assertRaises(ValueError): z.ritz_gap_interval((0, 0), (3, 3), 4, 1, 0)

    def test_exact_full_spectra_validate_noncommuting_spectral_enclosures(self):
        rows = z.matrix_controls()['rational_fixtures_not_Yang_Mills_vacua']
        self.assertEqual(len(rows), 3)
        for row in rows:
            self.assertLessEqual(row['certified_gap']['lower'], row['exact_full_gap'])
            self.assertGreaterEqual(row['certified_gap']['upper'], row['exact_full_gap'])
            self.assertNotEqual(row['nonzero_retained_coupling_commutator'][0][1], 0)
            for index in range(2):
                exact_eigen = row['exact_full_spectrum'][index]
                endpoints = row['certified_gap'][('ground_eigenvalue', 'first_excited_eigenvalue')[index]]
                self.assertLessEqual(endpoints['lower'], exact_eigen)
                self.assertGreaterEqual(endpoints['upper'], exact_eigen)

    def test_uncoupled_hidden_level_cannot_be_dropped(self):
        row = z.matrix_controls()['invisible_hidden_level']
        # diag(3,5) retained, hidden energy1: actual full gap2, retained gap2
        # happens to match, but ground eigenvalue1 is missing without the d cap.
        self.assertEqual(row['ground_eigenvalue']['lower'], 1)
        self.assertEqual(row['first_excited_eigenvalue']['lower'], 1)
        self.assertEqual(row['lower'], 0)

    def test_universal_width_bound_decreases_and_counts_arithmetic_uncertainty(self):
        rows = z.cutoff_controls()['fixed_regulator_width_convergence']
        self.assertGreater(rows[0]['physical_total_width_upper'], rows[-1]['physical_total_width_upper'])
        exact = z.ritz_width_budget(2, 100, 3)
        rounded = z.ritz_width_budget(2, 100, 3, F(1, 100))
        self.assertEqual(rounded['physical_total_width_upper']-exact['physical_total_width_upper'], F(3, 100))

    def test_parameter_error_includes_both_eigenvalues_and_both_interval_ends(self):
        row = z.coupling_parameter_budget(192, F(1, 10**6), F(3, 2))
        self.assertEqual(row['physical_energy_error_radius'], F(288, 10**6))
        self.assertEqual(row['physical_gap_error_radius'], F(576, 10**6))
        self.assertEqual(row['total_gap_interval_width_increase_upper'], F(1152, 10**6))
        self.assertTrue(row['physical_clock_held_fixed'])
        for count, radius, clock in ((0, 1, 1), (192, -1, 1), (192, 0.1, 1), (192, 1, 0)):
            with self.assertRaises(ValueError): z.coupling_parameter_budget(count, radius, clock)

    def test_parameter_lipschitz_factor_two_is_sharp_for_independent_matrix_family(self):
        potential_norm, theta, radius, clock = 7, F(2), F(1, 100), F(5, 3)
        # A(theta)=diag(-7theta,7theta), ||V||=7: its gap is exactly14theta.
        exact_change = clock*(2*potential_norm*(theta+radius)-2*potential_norm*theta)
        budget = z.coupling_parameter_budget(potential_norm, radius, clock)
        self.assertEqual(exact_change, budget['physical_gap_error_radius'])
        self.assertGreater(exact_change, budget['physical_energy_error_radius'])

    def test_four_cube_enumeration_includes_all_winding_observations(self):
        geometry = z.geometry_controls()
        self.assertEqual(geometry['simple_four_cycles'], 240)
        self.assertEqual(geometry['plaquettes'], 192)
        self.assertEqual(geometry['winding_four_cycles'], 48)
        self.assertEqual(geometry['physical_R3_dimension'], 1+192+48)
        self.assertEqual(geometry['distinct_plaquette_pair_xors'], 192*191//2)
        self.assertTrue(geometry['periodic_four_cube_graph_bipartite'])

    def test_native_Haar_fourth_moment_fixes_the_sharp_coupling_constant(self):
        row = z.geometry_controls()
        self.assertEqual(row['Haar_W_squared'], F(1, 4))
        self.assertEqual(row['Haar_W_fourth'], F(1, 8))
        self.assertEqual(row['BBstar_plaquette_identity_coefficient'], F(191, 4))
        self.assertEqual(row['BBstar_plaquette_all_ones_coefficient'], F(1, 4))
        self.assertEqual(row['BBstar_winding_identity_coefficient'], 48)
        self.assertEqual(row['BBstar_largest_eigenvalue_over_theta_squared'], F(383, 4))
        self.assertLess(F(383, 4), 192**2)

    def test_actual_cube_certifies_gap_outside_the_old_small_theta_window(self):
        row = z.four_cube_gap(F(1, 48))
        self.assertGreater(row['gap']['lower'], F(4904, 10000))
        self.assertLess(row['gap']['upper'], F(307874, 100000))
        self.assertEqual(row['exact_coupling_squared'], F(383, 4*48**2))
        self.assertTrue(row['actual_finite_cube'])
        self.assertFalse(row['evaluated_interacting_vacuum'])
        self.assertFalse(row['continuum_claim'])

    def test_cube_gap_bound_recovers_free_gap_and_can_be_inconclusive(self):
        free = z.four_cube_gap(0)['gap']
        self.assertEqual((free['lower'], free['upper']), (3, 3))
        inconclusive = z.four_cube_gap(F(1, 40))['gap']
        self.assertEqual(inconclusive['lower'], 0)
        self.assertFalse(inconclusive['positive_gap_certified'])
        self.assertGreater(inconclusive['upper'], 3)

    def test_cube_bounds_respect_sign_symmetry_and_physical_clock(self):
        positive = z.four_cube_gap(F(1, 48))
        negative = z.four_cube_gap(F(-1, 48))
        self.assertEqual(positive['gap'], negative['gap'])
        changed = z.four_cube_gap(F(1, 48), F(7, 3))
        for endpoint in ('lower', 'upper', 'width'):
            self.assertEqual(changed['gap'][endpoint], F(7, 3)*positive['gap'][endpoint])

    def test_required_control_groups_execute_without_a_source_manifest(self):
        self.assertEqual(len(z.actual_cube_controls()['actual_periodic_cube_gap_intervals']), 5)
        self.assertEqual(len(z.cutoff_controls()['all_coupling_cutoff_examples']), 4)
        self.assertTrue(z.trial_controls()['nonphysical_trial_not_used_for_excited_Ritz_eigenvalues'])


class ElectricEnclosureCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError): z.run()

    def test_253_upstream_pins_and_four_local_inputs(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 253)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], z.BASE_TREE)
        for name in ('YM73_PHYSICAL_MASS_MEMORY_METRIC.md', 'YM62_GAUGE_LATTICE_GAP.md'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+name, manifest['upstream_sha256'])

    def test_upstream_tampering_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError): z.source_checks()

    def test_rehashed_continuum_overclaim_is_rejected_even_as_fresh_input(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['continuum_R4_QFT'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError): z.check(wrong, result, pin)
            with self.assertRaises(ValueError): z.check(self.cert, result, pin)

    def test_rehashed_scope_inflation_cannot_create_a_uniform_gap_claim(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['uniform_weak_coupling_positive_gap'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError): z.check(wrong, result, pin)

    def test_certificate_has_no_binary_floats_or_unintended_control_characters(self):
        def visit(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): visit(child)
            elif isinstance(value, (list, tuple)):
                for child in value: visit(child)
        visit(self.cert)
        for path in json.loads(z.SOURCES.read_text())['local_inputs']:
            self.assertFalse(any(byte < 32 and byte not in (9, 10) for byte in (z.ROOT/path).read_bytes()))

    def test_scope_separates_actual_finite_model_advances_and_continuum_inputs(self):
        claims = self.cert['model']
        for field in ('all_coupling_hidden_floor', 'whole_low_energy_capture',
                      'fixed_regulator_arbitrary_precision',
                      'actual_all_coupling_vacuum_energy_and_flatness_bound',
                      'actual_periodic_four_cube_241_dimensional_cut', 'winding_four_cycles_retained'):
            self.assertTrue(claims[field])
        for field in ('uniform_weak_coupling_positive_gap', 'continuum_R4_QFT',
                      'physical_mass_prediction', 'finite_fixture_is_Yang_Mills_vacuum',
                      'coverage_implies_nonvacuum_low_energy_states_exist',
                      'plaquette_flatness_controls_continuum_curvature'):
            self.assertFalse(claims[field])


if __name__ == '__main__':
    unittest.main()
