import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym65_spatial_response_locality as z


class SpatialResponseLocalityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_certificate_and_hash_match(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_is_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_190_upstream_sources_are_pinned(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 190)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn('papers/yang-mills-certified-benchmark/YM64_VACUUM_ADAPTED_TOWER.md',
                      manifest['upstream_sha256'])

    def test_changed_upstream_and_rehashed_claim_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            changed = copy.deepcopy(self.cert)
            changed['model']['Schur_memory_kernel_locality'] = True
            result.write_text(json.dumps(changed))
            pin.write_text(z.p.canonical_sha(changed))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_native_incidence_counts_the_four_diagonal_multiplicities(self):
        _, rows, _ = z.graph(4)
        self.assertEqual({row[i] for i, row in enumerate(rows)}, {4})
        self.assertEqual({sum(row.values()) for row in rows}, {16})
        self.assertEqual({sum(row.values())-row[i] for i, row in enumerate(rows)}, {12})

    def test_open_boundaries_only_reduce_the_incidence_budget(self):
        _, rows, _ = z.graph(4, False)
        sums = [sum(row.values()) for row in rows]
        self.assertEqual(max(sums), 16)
        self.assertLess(min(sums), 16)

    def test_independent_hyperedge_chains_are_overcounted_by_link_walks(self):
        data = self.cert['incidence']
        distance = data['independent_chain_distance']
        self.assertGreater(distance, 1)
        for row in data['independent_chain_vs_walk_counts']:
            self.assertLessEqual(row['plaquette_chains'], row['link_walks'])
            if row['order'] < distance:
                self.assertEqual(row['link_walks'], 0)
            if row['order'] == distance:
                self.assertGreater(row['link_walks'], 0)

    def test_every_link_walk_step_moves_midpoints_by_at_most_one(self):
        self.assertTrue(self.cert['incidence']['midpoint_distance_changes_by_at_most_one'])
        _, rows, _ = z.graph(5)
        initial = {0: 1}
        for n in range(6):
            self.assertEqual(sum(initial.values()), 16**n)
            initial = z.advance(rows, initial)

    def test_propagation_bound_exists_outside_the_gap_window(self):
        bound = z.lr_upper(1, 20, 1, F(1, 100))
        self.assertGreater(bound, 0)
        self.assertLess(bound, F(1, 1000))
        with self.assertRaises(ValueError):
            z.parameters(1)

    def test_propagation_respects_sign_time_reversal_and_the_zero_coupling_limit(self):
        value = z.lr_upper(F(1, 16000), 12, 4, 30)
        self.assertEqual(value, z.lr_upper(F(-1, 16000), 12, 4, -30))
        self.assertEqual(z.lr_upper(0, 1, 100, 100), 0)
        self.assertEqual(z.lr_upper(F(1, 16000), 1, 100, 0), 0)
        self.assertEqual(z.lr_upper(1, 1, 100, 100), 2)

    def test_spatial_rate_uses_both_the_gap_and_propagation_speed(self):
        cell = z.parameters(F(1, 16000))
        self.assertEqual(cell['velocity'], F(1, 250))
        self.assertGreater(cell['rho'].lo, F(1928159848815, 10**13))
        self.assertLess(cell['rho'].hi, F(1928159848816, 10**13))
        rates = [z.parameters(theta)['rho'] for theta in (F(1, 32000), F(1, 16000), F(1, 10000))]
        self.assertTrue(all(x.lo > y.hi for x, y in zip(rates, rates[1:])))

    def test_balance_relations_match_all_three_filter_exponents(self):
        cell = z.parameters(F(1, 16000))
        distance = 32
        r = cell['rho']*z.Iv(distance)
        T = z.Iv(2)*r/cell['gamma']
        a = cell['gamma']/(z.Iv(2)*T)
        first = cell['gamma']*T/z.Iv(2)
        second = cell['log_two']*z.Iv(distance)-z.Iv(cell['velocity'])*T
        third = a*T*T
        for value in (first, second, third):
            self.assertFalse(value.separated_from(r))

    def test_machin_circle_identity_and_gaussian_reflection_are_enclosed(self):
        pi = z.pi_interval()
        self.assertGreater(pi.lo, F(314159265358979323846, 10**20))
        self.assertLess(pi.hi, F(314159265358979323847, 10**20))
        self.assertEqual((z.erfc_interval(0).lo, z.erfc_interval(0).hi), (1, 1))
        for x in (F(1, 4), F(1), F(3)):
            positive, negative = z.erfc_interval(x), z.erfc_interval(-x)
            self.assertFalse((positive+negative).separated_from(z.Iv(2)))

    def test_gaussian_filter_has_the_correct_frequency_and_heat_signs(self):
        for row in self.cert['filter']['scalar_spectral_filters']:
            self.assertGreater(F(row['error_upper']), 0)
            self.assertLess(F(row['error_upper']), F(row['proved_error_upper']))
        correct, target = z.scalar_filter(3, 1), z.refine.ex(-3)
        wrong = correct*z.refine.ex(F(1, 4))
        self.assertGreater((wrong-target).lo, z.exp_negative(F(37, 4)).hi)

    def test_gaussian_control_refuses_unproved_ranges_and_inexact_inputs(self):
        for args in ((0, 0), (1, -1), (1, 3), (4, 1)):
            with self.assertRaises(ValueError):
                z.scalar_filter(*args)
        for x in (True, 0.1, 5):
            with self.assertRaises(ValueError):
                z.erfc_interval(x)

    def test_full_inverse_response_tower_matches_an_independent_energy_hessian(self):
        data = self.cert['sources']
        moments = [F(value) for value in data['inverse_moments']]
        self.assertEqual(moments, [F(576, 625)*F(1, 2)**n for n in range(7)])
        self.assertEqual(2*moments[1], F(data['independent_energy_hessian']))
        self.assertEqual(F(data['independent_energy_hessian']), F(576, 625))

    def test_higher_inverse_orders_have_nontrivial_spatial_bounds_at_large_distance(self):
        rows = {row['inverse_order']: row for row in self.cert['budgets'] if row['distance'] == 128}
        self.assertLess(F(rows[0]['normalized_upper']), F(7036, 10**14))
        self.assertLess(F(rows[1]['normalized_upper']), F(3473, 10**12))
        self.assertLess(F(rows[2]['normalized_upper']), F(8571, 10**11))
        self.assertLess(F(rows[4]['normalized_upper']), F(1741, 10**8))

    def test_gap_normalization_and_clock_scaling_are_not_schur_moment_scaling(self):
        kappa = F(7, 3)
        for n in (0, 1, 2, 4):
            base = z.response_budget(F(1, 16000), 64, 4, n)
            changed = z.response_budget(F(1, 16000), 64, 4, n, kappa)
            self.assertEqual(changed['normalized_upper'], base['normalized_upper'])
            self.assertEqual(changed['absolute_upper'], base['absolute_upper']/kappa**n)
            self.assertEqual(changed['rate'].lo, base['rate'].lo)
            self.assertFalse((changed['cutoff_time']*z.Iv(kappa)).separated_from(base['cutoff_time']))

    def test_free_disjoint_readouts_have_zero_connected_inverse_response(self):
        for n in (0, 1, 2, 8):
            result = z.response_budget(0, 4, 4, n)
            self.assertEqual(result['absolute_upper'], 0)

    def test_all_envelopes_are_clipped_by_the_gap_norm_bound(self):
        for row in self.cert['budgets']:
            self.assertGreaterEqual(F(row['normalized_upper']), 0)
            self.assertLessEqual(F(row['normalized_upper']), 1)

    def test_zero_distance_bad_support_and_closed_gap_endpoint_are_refused(self):
        for args in ((F(1, 8000), 4, 4, 0), (F(1, 16000), 0, 4, 0),
                     (F(1, 16000), 4, False, 0), (F(1, 16000), 4, 4, True)):
            with self.assertRaises(ValueError):
                z.response_budget(*args)
        with self.assertRaises(ArithmeticError):
            z.parameters(F(1, 8000)-F(1, 10**40))

    def test_geometric_shell_sums_are_exact_for_all_required_powers(self):
        self.assertEqual(z.geometric_moments(F(1, 2)), (1, 2, 6, 26, 150))
        self.assertEqual(z.geometric_moments(0), (0, 0, 0, 0, 0))
        for q in (-1, 1, 0.5):
            with self.assertRaises(ValueError):
                z.geometric_moments(q)

    def test_susceptibility_is_summable_with_a_finite_conservative_budget(self):
        bound = F(self.cert['susceptibility_row_sum_upper'])
        self.assertGreater(bound, 3*10**9)
        self.assertLess(bound, 3005479*10**3)
        self.assertEqual(bound, z.susceptibility_sum_upper(F(-1, 16000)))

    def test_simplified_susceptibility_envelope_dominates_the_inverse_bound(self):
        cell = z.parameters(F(1, 16000))
        r = cell['rho']*z.Iv(32)
        upper = z.Iv(2)*(z.Iv(22)*r+z.Iv(5))*z.exp_negative(r)/(z.Iv(3)*cell['gamma'])
        actual_bound = z.response_budget(F(1, 16000), 32, 4, 1)['absolute_upper']*2
        self.assertGreater(upper.lo, actual_bound)

    def test_tiling_removes_only_boundary_order_plaquettes(self):
        for row in self.cert['tiling']['tilings']:
            side = row['block_side']*row['blocks_per_axis']
            self.assertLessEqual(row['crossing_faces'], F(6*side**3,row['block_side']))
        self.assertEqual(self.cert['tiling']['independent_vertex_partition_internal_faces'], 288)

    def test_open_to_periodic_energy_comparison_counts_every_wrapping_face(self):
        for row in self.cert['tiling']['periodic_wrap_faces']:
            L = row['side']
            self.assertEqual(row['wrap_faces'], 3*L*(2*L-1))
            self.assertLessEqual(row['wrap_faces'], 6*L*L)

    def test_bulk_mean_budget_uses_both_boundary_error_and_uniform_susceptibility(self):
        theta_max, side, step = F(1, 16000), 10**16, F(1, 10**14)
        expected = 4*theta_max/side/step+z.susceptibility_sum_upper(theta_max)*step
        self.assertEqual(z.bulk_mean_error(theta_max, side, step), expected)
        self.assertLess(expected, F(1, 10000))
        with self.assertRaises(ValueError):
            z.bulk_mean_error(theta_max, 4, theta_max, theta_max)

    def test_rare_conditioning_amplifies_uniformly_small_unconditional_covariance(self):
        rows = self.cert['conditioning']['examples']
        self.assertEqual({F(row['conditioned_covariance']) for row in rows}, {F(1, 2)})
        self.assertEqual([F(row['all_bounded_left_right_covariance_bound']) for row in rows],
                         [F(1, 2**(d+1)) for d in (4, 8, 16)])

    def test_spatial_decay_dominates_a_polynomial_readout_normalization(self):
        theta = F(1, 16000)
        sequence = [d**2*z.response_budget(theta, d, 4)['absolute_upper'] for d in (32, 64, 128)]
        self.assertTrue(all(x > y for x, y in zip(sequence, sequence[1:])))
        self.assertLess(sequence[-1], F(1, 10**5))

    def test_scope_separates_resolvent_locality_bulk_regularity_and_remaining_gates(self):
        model = self.cert['model']
        for key in ('YM62_Hamiltonian_retained', 'unbounded_on_site_terms_retained',
                    'volume_uniform_propagation', 'volume_uniform_connected_clustering',
                    'every_fixed_inverse_response_order_local', 'blocked_resolvent_readout_bounds',
                    'bulk_plaquette_mean_unique', 'bulk_energy_C1_1'):
            self.assertTrue(model[key])
        for key in ('uniform_conditional_density_locality', 'Schur_memory_kernel_locality',
                    'bulk_second_derivative_limit', 'continuum_R4_QFT',
                    'physical_light_speed_derived', 'physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
