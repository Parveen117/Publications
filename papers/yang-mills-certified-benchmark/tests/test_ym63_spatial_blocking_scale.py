import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym63_spatial_blocking_scale as z


class SpatialBlockingScaleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_evidence_and_hash_are_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_energy_clock_and_magnetic_coefficients_match_exactly(self):
        data = z.scale_map(F(1, 7), F(5, 3), F(3, 2), F(7))
        self.assertEqual(data['kappa'], F(35, 2))
        self.assertEqual(data['theta'], F(42, 25))
        self.assertEqual(data['kappa']*data['theta'], F(147, 5))
        self.assertEqual(data['magnetic_coefficient'], F(147, 5))

    def test_gap_scales_linearly_under_positive_clock_calibration(self):
        base = z.scaled_gap(F(1, 16000), 1)
        changed = z.scaled_gap(F(1, 16000), F(7, 3))
        self.assertEqual(changed['lower'], F(7, 3)*base['lower'])
        self.assertEqual(changed['upper'], F(700, 57))

    def test_fixed_bare_coupling_forces_inverse_spacing_gap_bounds(self):
        rows = self.cert['scales']['fixed_coupling']
        lower = [F(r['spacing'])*F(r['gap_lower']) for r in rows]
        upper = [F(r['spacing'])*F(r['gap_upper']) for r in rows]
        self.assertEqual(len(set(lower)), 1)
        self.assertEqual(set(upper), {F(25600, 19)})
        self.assertTrue(all(r['in_YM62_window'] for r in rows))

    def test_holding_the_clock_fixed_exits_the_certified_coupling_window(self):
        rows = self.cert['scales']['fixed_clock_leaves_window']
        self.assertEqual({F(r['kappa']) for r in rows}, {F(256)})
        self.assertEqual([r['in_YM62_window'] for r in rows], [True, False, False, False])
        self.assertEqual({F(r['theta'])*F(r['spacing'])**2 for r in rows}, {F(1, 16384)})

    def test_endpoint_and_weak_bare_coupling_do_not_receive_a_gap_certificate(self):
        self.assertIsNone(z.scaled_gap(F(1, 8000), 1))
        for coupling in (F(1), F(1, 2), F(1, 1000)):
            cell = z.scale_map(F(1, 100), coupling)
            self.assertFalse(cell['in_YM62_window'])
            self.assertIsNone(z.scaled_gap(cell['theta'], cell['kappa']))

    def test_convention_constants_remain_explicit_and_change_the_ratio(self):
        first = z.scale_map(1, 256, 1, 4)
        second = z.scale_map(1, 256, F(1, 4), 4)
        self.assertEqual(second['theta'], 4*first['theta'])
        self.assertEqual(4*second['kappa'], first['kappa'])
        self.assertTrue(first['in_YM62_window'])
        self.assertFalse(second['in_YM62_window'])

    def test_product_link_preserves_native_reference_coordinate_moments(self):
        coords = z.product_coordinates(3)
        self.assertEqual(z.native.phi(coords[0], 3), 0)
        self.assertEqual(z.native.phi(z.p.mul(coords[0], coords[0]), 3), F(1, 4))
        self.assertEqual(z.native.phi(z.p.mul(coords[0], coords[1]), 3), 0)

    def test_free_blocking_clock_is_b_not_b_squared_or_inverse_b(self):
        f = z.p.mul(z.p.COORD[0], z.p.COORD[1])
        image = z.pullback(f, 3)
        base = z.pullback(z.y.lap(f), 3)
        actual = z.previous.electric(image, 3)
        self.assertEqual(actual, z.p.scale(base, 3))
        self.assertNotEqual(actual, z.p.scale(base, 9))
        self.assertNotEqual(actual, z.p.scale(base, F(1, 3)))

    def test_intermediate_gauge_factors_cancel_in_the_ordered_product(self):
        q = z.y.turns()[:2]
        g = (z.y.turns()[3], z.y.turns()[5], z.y.turns()[1])
        changed = z.previous.gauge_transform(q, g, ((0, 1), (1, 2)))
        product = z.y.qmul(q[0], q[1])
        expected = z.y.qmul(z.y.qmul(g[0], product), z.y.qdagger(g[2]))
        self.assertEqual(z.y.qmul(changed[0], changed[1]), expected)

    def test_two_successive_blocks_have_the_same_fine_paths_as_one_combined_block(self):
        direct = z.blocked_geometry(4, 4)
        first, second = z.blocked_geometry(8, 2), z.blocked_geometry(4, 2)
        composed = tuple(tuple(e for middle in path for e in first['paths'][middle])
                         for path in second['paths'])
        self.assertEqual(composed, direct['paths'])

    def test_each_fine_plaquette_has_a_link_outside_all_retained_paths(self):
        data = z.blocked_geometry(4, 2)
        paths = [e for path in data['paths'] for e in path]
        self.assertEqual(len(paths), len(set(paths)))
        self.assertTrue(all(data['unretained_per_face']))
        self.assertEqual(len(data['fine']['faces']), 1536)

    def test_unretained_link_averaging_annihilates_the_first_order_interaction(self):
        w = z.previous.loop_polynomial()
        self.assertFalse(z.source.site_average(w, 1, 4))
        self.assertGreater(z.native.phi(z.p.mul(w, w), 4), 0)

    def test_second_order_return_is_positive_even_when_first_order_vanishes(self):
        self.assertEqual(z.memory_norm(1536, F(1, 16000)), F(3, 2000000))
        self.assertEqual(z.memory_norm(1536, F(-1, 16000)), z.memory_norm(1536, F(1, 16000)))
        self.assertEqual(z.memory_norm(1536, 0), 0)

    def test_schur_correction_has_the_subtraction_sign(self):
        denominator = F(2)*3-F(1, 3)**2
        direct_resolvent = F(3)/denominator
        self.assertEqual(direct_resolvent, 1/(F(2)-F(1, 3)**2/3))
        self.assertNotEqual(direct_resolvent, 1/(F(2)+F(1, 3)**2/3))

    def test_shared_link_contraction_matches_a_full_native_boundary_polynomial(self):
        data = self.cert['contractions']
        self.assertTrue(data['shared_link_contraction'])
        self.assertEqual(data['outer_terms'], 1024)
        self.assertEqual(F(data['outer_rate']), F(9, 2))
        self.assertEqual(F(data['two_face_full_polynomial_response']), F(7, 216))

    def test_averaging_before_the_second_incident_insertion_erases_the_source(self):
        patch = z.rectangle(2, 1)
        n = len(patch['edges'])
        left = z.previous.loop_polynomial(patch['faces'][0], n)
        right = z.previous.loop_polynomial(patch['faces'][1], n)
        shared = set(e for e, sign in patch['faces'][0]) & set(e for e, sign in patch['faces'][1])
        self.assertFalse(z.source.site_average(left, next(iter(shared)), n))
        self.assertTrue(z.multiply_average(left, right, shared, n))

    def test_bra_ket_and_differentiated_energy_recursions_agree(self):
        rows = self.cert['response']['rectangles']
        self.assertEqual(len(rows), 7)
        for row in rows:
            self.assertEqual(F(row['leading']), F(row['energy_route']))
            self.assertEqual(row['proper_inverse_witnesses'], 2**(row['area']+1)-2)
            self.assertEqual(F(row['free_surface_moment']), F(1, 4**row['area']))

    def test_equal_area_different_shape_has_different_response(self):
        strip, square = z.loop_response(1, 4), z.loop_response(2, 2)
        self.assertEqual(strip['area'], square['area'])
        self.assertEqual(strip['leading'], F(73, 58320))
        self.assertEqual(square['leading'], F(35, 23328))
        self.assertNotEqual(strip['leading'], square['leading'])

    def test_independent_plaquette_factorization_misses_the_exact_coefficient(self):
        ratio = z.loop_response(2, 2)['leading']/z.loop_response(1, 1)['leading']**4
        self.assertEqual(ratio, F(35, 18))
        self.assertNotEqual(ratio, 1)

    def test_no_proper_label_subset_of_the_closed_sheet_has_zero_frontier(self):
        patch = z.rectangle(2, 2)
        boundaries = z.boundary_masks(patch['faces']+(patch['outer'],))
        self.assertEqual([i for i, value in enumerate(boundaries) if not value], [0, 31])
        self.assertGreater(len({v.bit_count() for v in boundaries[1:-1]}), 1)

    def test_one_fixed_inverse_rate_cannot_replace_changing_frontiers(self):
        # Five distinct labels: using rate 3 at every proper insertion gives
        # 5! / (3^4 * 4^4), a different coefficient.
        wrong = F(120, 3**4*4**4)
        self.assertNotEqual(wrong, z.loop_response(2, 2)['leading'])

    def test_center_involution_fixes_response_parity_and_the_next_error_order(self):
        data = self.cert['parity']
        self.assertEqual(data['open_box_face_flips'], 150)
        self.assertEqual(data['rectangle_parity_checks'], 16)
        self.assertTrue(data['odd_periodic_extension_refused'])
        for m, n in ((1, 1), (1, 2), (2, 2), (3, 3)):
            self.assertEqual(z.loop_response(m, n)['fixed_volume_remainder_order'], m*n+2)

    def test_rotating_the_rectangle_preserves_its_coefficient(self):
        self.assertEqual(z.loop_response(2, 3)['leading'], z.loop_response(3, 2)['leading'])
        self.assertEqual(z.loop_response(1, 4)['leading'], z.loop_response(4, 1)['leading'])

    def test_coarse_potential_matching_includes_the_free_clock_factor(self):
        normalized = z.loop_response(2, 2)['matching_normalized_coarse_coefficient']
        self.assertEqual(normalized, F(35, 3888))
        self.assertEqual(2*normalized, F(35, 1944))
        self.assertTrue(self.cert['response']['matching_is_only_one_observable'])

    def test_elementary_loop_response_reproduces_the_previous_chapter(self):
        previous = z.previous.plaquette_controls()['diagonal_source_susceptibility']
        self.assertEqual(z.loop_response(1, 1)['leading'], previous)

    def test_closed_proper_subsets_and_malformed_boundaries_are_refused(self):
        word = z.rectangle(1, 1)['faces'][0]
        with self.assertRaises(ValueError):
            z.inverse_weights((word, word))
        with self.assertRaises(ValueError):
            z.boundary_masks((((0, 1), (0, -1)),))
        with self.assertRaises(ValueError):
            z.boundary_masks(())

    def test_invalid_parameters_and_excessive_exact_enumerations_are_refused(self):
        for value in (0, -1, True, 0.1):
            with self.assertRaises(ValueError):
                z.scale_map(value, 1)
            with self.assertRaises(ValueError):
                z.scale_map(1, value)
        for dims in ((0, 1), (True, 2), (2.0, 2), (1, 13)):
            with self.assertRaises(ValueError):
                z.loop_response(*dims)
        with self.assertRaises(ValueError):
            z.blocked_geometry(4, 1)

    def test_python_312_is_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_upstream_and_full_rg_claim_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            changed = copy.deepcopy(self.cert)
            changed['model']['exact_one_parameter_RG_closure'] = True
            result.write_text(json.dumps(changed))
            pin.write_text(z.native.previous.native.chain.old.canonical_sha(changed))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_all_175_upstream_files_and_four_new_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 175)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn('papers/generalized-euler-evolution/GE3_RETURNING_MEMORY.md', manifest['upstream_sha256'])

    def test_scope_distinguishes_retained_model_and_finite_response_from_continuum(self):
        model = self.cert['model']
        for key in ('YM62_operator_retained', 'free_reference_spatial_cut', 'gauge_compatible_product_links',
                    'exact_free_heat_intertwining', 'interacting_return_memory',
                    'leading_rectangle_response', 'two_sided_scaled_gap_constraints'):
            self.assertTrue(model[key])
        for key in ('exact_one_parameter_RG_closure', 'uniform_area_law',
                    'interacting_reference_conditional_expectation', 'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
