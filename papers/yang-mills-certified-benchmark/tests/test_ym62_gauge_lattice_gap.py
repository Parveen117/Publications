import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym62_gauge_lattice_gap as z


class GaugeLatticeGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.w = z.loop_polynomial()

    def test_fresh_evidence_and_hash_are_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_each_vertex_and_colour_obeys_the_gauss_and_casimir_identities(self):
        self.assertEqual(self.cert['gauge']['gauss_ward_identities'], 12)
        self.assertEqual(self.cert['gauge']['casimir_gauss_commutators'], 12)

    def test_finite_gauge_word_transforms_by_base_vertex_conjugation(self):
        qs = (z.y.ONE, z.y.turns()[2], z.y.turns()[5], (F(5, 13), F(12, 13), F(0), F(0)))
        gs = (z.y.turns()[1], z.y.ONE, z.y.turns()[4], z.y.turns()[3])
        original, transformed = z.loop_value(qs), z.loop_value(z.gauge_transform(qs, gs))
        self.assertEqual(transformed, z.y.qmul(z.y.qmul(gs[0], original), z.y.qdagger(gs[0])))
        self.assertEqual(z.evaluate(self.w, qs), transformed[0])

    def test_wrong_order_has_an_explicit_noninvariance_witness(self):
        data = self.cert['gauge']
        self.assertEqual(F(data['wrong_order_before']), F(-37, 125))
        self.assertEqual(F(data['wrong_order_after']), F(-29, 3125))

    def test_old_open_overlap_is_not_a_local_gauge_scalar(self):
        data = self.cert['gauge']
        self.assertEqual(F(data['open_overlap_before']), F(9, 25))
        self.assertEqual(F(data['open_overlap_after']), F(-351, 625))

    def test_incoming_gauss_sign_cannot_be_replaced_by_plus(self):
        self.assertFalse(z.gauss(self.w, 1, 2))
        wrong = z.p.add(z.native.derivative(self.w, 1, 2), z.right_derivative(self.w, 0, 2))
        self.assertTrue(wrong)

    def test_rank_one_colour_kinetic_operator_fails_local_rotation(self):
        self.assertTrue(self.cert['gauge']['rank_one_colour_operator_rejected'])
        f = z.p.mul(z.p.COORD[0], z.p.COORD[0])
        D = z.y.native_generator
        self.assertNotEqual(D(D(D(f, 2), 1), 1), D(D(D(f, 1), 1), 2))

    def test_four_links_have_electric_rate_three_and_free_variance_quarter(self):
        self.assertEqual(z.electric(self.w, 4), z.p.scale(self.w, 3))
        self.assertEqual(z.native.phi(self.w, 4), 0)
        self.assertEqual(z.native.phi(z.p.mul(self.w, self.w), 4), F(1, 4))

    def test_one_link_gradient_has_the_exact_quaternion_radius_identity(self):
        gradient = z.p.add(*(z.p.mul(df, df) for axis in (1, 2, 3)
                             for df in (z.native.derivative(self.w, 0, axis),)))
        radius = z.reduce(z.p.mul, [z.native.lift(z.y.RADIUS, e, 4) for e in range(4)])
        target = z.p.scale(z.p.add(radius, z.p.scale(z.p.mul(self.w, self.w), -1)), F(1, 4))
        self.assertEqual(gradient, target)

    def test_half_casimir_normalization_is_detectably_different(self):
        self.assertNotEqual(z.p.scale(z.electric(self.w, 4), F(1, 2)), z.p.scale(self.w, 3))
        self.assertEqual(z.y.lap(z.p.COORD[0]), z.p.scale(z.p.COORD[0], F(3, 4)))

    def test_cubic_geometry_has_three_spatial_directions(self):
        cell = z.lattice(4)
        self.assertEqual(len(cell['vertices']), 64)
        self.assertEqual(len(cell['edges']), 192)
        self.assertEqual(len(cell['faces']), 192)
        degrees, pairs = z.incidence(cell['faces'], 192)
        self.assertEqual(set(degrees), {4})
        self.assertEqual(sum(v for (i, j), v in pairs.items() if i == 0), 12)
        open_cell = z.lattice(2, False)
        self.assertEqual((len(open_cell['edges']), len(open_cell['faces'])), (12, 6))

    def test_reversed_links_close_each_plaquette_at_its_base_vertex(self):
        cell = z.lattice(4)
        for word in cell['faces']:
            edges = [cell['ends'][e] if sign == 1 else cell['ends'][e][::-1] for e, sign in word]
            self.assertEqual([x[1] for x in edges], [x[0] for x in edges[1:]+edges[:1]])

    def test_adjacent_plaquettes_have_zero_free_quadratic_cross_response(self):
        left = z.loop_polynomial(z.ORIENTED_WORD, 7)
        right = z.loop_polynomial(((0, 1), (4, 1), (5, -1), (6, -1)), 7)
        self.assertEqual(z.native.phi(z.p.mul(left, right), 7), 0)
        self.assertEqual(z.native.phi(z.p.mul(left, left), 7), F(1, 4))

    def test_response_sign_and_factor_follow_the_energy_potential(self):
        data = self.cert['plaquette']
        e2 = F(data['quadratic_energy_per_plaquette'])
        self.assertEqual(e2, F(-1, 12))
        self.assertEqual(-2*e2, F(data['diagonal_source_susceptibility']))
        self.assertEqual(3*e2, F(data['bulk_quadratic_energy_per_cell']))

    def test_gate_is_sign_symmetric_and_has_positive_outward_rate(self):
        for theta in (F(0), F(1, 16000), F(999, 8000000)):
            positive, negative = z.block_gate(theta), z.block_gate(-theta)
            self.assertEqual(positive['alpha'], negative['alpha'])
            self.assertGreater(positive['gamma'].lo, 0)
            self.assertGreater(positive['margin'], 0)

    def test_strict_endpoint_and_every_larger_budget_are_refused(self):
        for theta in (F(1, 8000), F(1, 4000), F(1)):
            for J in range(1, 31):
                self.assertIsNone(z.block_gate(theta, J))
        self.assertIsNone(z.block_gate(F(1, 16000), R=F(4, 3)))

    def test_two_body_chain_budget_is_recovered_without_changing_its_engine(self):
        new = z.block_gate(F(1, 4096), 6, F(5, 4), 2, 2)
        old = z.blocks.block_gate(1, F(1, 4096), 6, F(5, 4))
        self.assertEqual(new['alpha'], old['alpha'])
        self.assertEqual(new['margin'], old['margin'])

    def test_shared_link_pair_multiplicity_is_not_collapsed(self):
        faces = (((0, 1), (1, 1), (2, 1)), ((0, 1), (1, 1), (3, 1)))
        degree, pairs = z.incidence(faces, 4)
        self.assertEqual(pairs[0, 1], 2)
        self.assertEqual(sum(v for (i, j), v in pairs.items() if i == 0), 2*degree[0])
        with self.assertRaises(ValueError):
            z.incidence((((0, 1), (0, -1)),), 2)

    def test_clipped_block_duplicates_keep_the_full_removal_budget(self):
        labels = z.blocks.old.blocks(13, 1, 5)
        self.assertEqual(len(labels), 65)
        self.assertEqual(len(set(labels)), 13)
        self.assertEqual(self.cert['incidence']['weighted_coordinate_checks'], 364)
        self.assertEqual(self.cert['incidence']['total_fine_coefficient'], 80)

    def test_coarse_skeleton_and_gate_cover_arbitrarily_fine_steps(self):
        h = F(1, 10**12)
        r = z.blocks.skeleton(1, h)
        ell = 5*r
        gate = z.block_gate(F(1, 16000))
        self.assertLessEqual(2*r/(z.EPS*ell)+80*h*F(1, 16000)*ell, gate['alpha'])
        self.assertEqual(gate['R'], F(7, 6))
        self.assertEqual(gate['margin'], F(1, 8))

    def test_time_refinement_is_fixed_volume_and_zero_at_zero_potential(self):
        self.assertEqual(z.refinement_budget(192, 0, 1, F(1, 64)), 0)
        self.assertGreater(z.refinement_budget(192, F(1, 16000), 1, F(1, 64)),
                           z.refinement_budget(6, F(1, 16000), 1, F(1, 64)))
        self.assertFalse(self.cert['refinement']['uniform_in_volume_refinement_error'])

    def test_explicit_example_rate_is_enclosed_and_not_rounded_up(self):
        row = self.cert['gap']['cells'][2]
        self.assertEqual(F(row['theta']), F(1, 16000))
        self.assertGreater(F(row['gamma_lower']), F('0.00308301'))
        self.assertLess(F(row['gamma_lower']), F(row['gamma_upper']))

    def test_physical_loop_variance_is_uniformly_nonzero(self):
        data = self.cert['histories']
        self.assertEqual(F(data['local_projection_error_upper']), F(1, 1500))
        self.assertGreater(F(data['variance_lower']), F(19, 100))
        self.assertEqual(F(data['stated_physical_plaquette_variance_floor']), F(19, 100))

    def test_history_gap_passes_without_mixing_different_vacuum_means(self):
        data = self.cert['histories']
        self.assertTrue(data['periodic_states_are_not_averaged'])
        self.assertTrue(data['uniform_gap_passes_on_all_positive_time_cylinders'])
        self.assertFalse(data['uniqueness_of_subsequential_limits'])
        self.assertFalse(data['time_zero_markov_limit'])

    def test_invalid_geometry_and_inexact_couplings_are_refused(self):
        for side in (0, 3, True, 4.0):
            with self.assertRaises(ValueError):
                z.lattice(side)
        for theta in (0.1, True):
            with self.assertRaises(ValueError):
                z.block_gate(theta)
        for J in (0, -1, True, 5.0):
            with self.assertRaises(ValueError):
                z.block_gate(0, J)
        with self.assertRaises(ValueError):
            z.refinement_budget(6, 1, 1, 2)

    def test_python_312_is_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_upstream_and_continuum_claim_tampering_are_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            bad = copy.deepcopy(self.cert)
            bad['model']['continuum_R4_QFT'] = True
            result.write_text(json.dumps(bad))
            pin.write_text(z.native.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_all_167_upstream_files_and_four_new_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 167)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn(str(z.previous.RESULT.relative_to(z.ROOT)), manifest['upstream_sha256'])

    def test_fixed_spatial_lattice_and_new_model_scope_remain_explicit(self):
        model = self.cert['model']
        self.assertEqual((model['spatial_dimensions'], model['continuous_heat_time_dimensions']), (3, 1))
        for key in ('native_links', 'independent_casimir_electric_operator', 'ordered_plaquette_potential',
                    'local_gauss_sector', 'volume_uniform_interacting_gap',
                    'subsequential_infinite_volume_gapped_histories', 'nontrivial_physical_loop_variance'):
            self.assertTrue(model[key])
        for key in ('unchanged_YM61_correlated_chain', 'continuum_R4_QFT',
                    'all_compact_simple_gauge_groups', 'physical_mass_prediction'):
            self.assertFalse(model[key])


if __name__ == '__main__':
    unittest.main()
