import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import selection_model as sm
from gravity_certificate import in_span, warped_chart
from selection_certificate import controls, metric_curvature


class SelectionTests(unittest.TestCase):
    def test_minimum_coupling_tie_has_three_distinct_signatures(self):
        result = controls()['W17_minimum_edge_tie']
        self.assertEqual(result['edge_count_each'], 3)
        self.assertEqual(set(result['signatures'].values()),
                         {(1, 3, 0), (2, 2, 0), (4, 0, 0)})

    def test_odd_K_cycle_eliminates_all_symmetric_forms(self):
        edges = ((0, 1, 1), (0, 2, 1), (1, 2, 1))
        self.assertIsNone(sm.coupling_form(3, edges))
        maps = [gm.cayley(sm.pair_generator(3, *edge), Q(1, 2)) for edge in edges]
        self.assertEqual(gm.invariant_forms(3, maps), ())

    def test_balanced_cycle_has_exactly_one_form_line(self):
        edges = ((0, 1, 1), (0, 2, 1), (1, 2, -1))
        form = sm.coupling_form(3, edges)
        maps = [gm.cayley(sm.pair_generator(3, *edge), Q(-2, 3)) for edge in edges]
        basis = gm.invariant_forms(3, maps)
        self.assertEqual(len(basis), 1)
        self.assertTrue(in_span(form, basis))

    def test_edge_order_does_not_change_normalization(self):
        edges = ((0, 3, -1), (2, 3, 1), (1, 2, -1))
        self.assertEqual(sm.coupling_form(4, edges), sm.coupling_form(4, edges[::-1]))

    def test_disconnected_graph_is_outside_connected_theorem(self):
        with self.assertRaisesRegex(ValueError, 'connected'):
            sm.coupling_form(4, ((0, 1, 1),))
        # Do not misclassify inconsistent component plus isolated vertex as U21.
        with self.assertRaisesRegex(ValueError, 'connected'):
            sm.coupling_form(4, ((0, 1, 1), (0, 2, 1), (1, 2, 1)))

    def test_duplicate_or_invalid_edges_do_not_silently_select_a_sector(self):
        for edges in (((0, 1, 1), (0, 1, -1)), ((1, 0, 1),),
                      ((0, 2, 1),), ((0, 1, 0),)):
            with self.assertRaises(ValueError):
                sm.coupling_form(2, edges)

    def test_nonunit_lapse_changes_required_scale_derivative(self):
        for sector in (-1, 1):
            good = sm.homogeneous_data(4, 2, 3, Q(2, 3), Q(1, 3), sector)
            bad = sm.homogeneous_data(4, 2, 3, Q(1, 3), Q(1, 3), sector)
            self.assertTrue(all(not any(x) for x in good['torsion'].values()))
            self.assertEqual(bad['torsion']['0,2'][2], Q(-1, 3))

    def test_pointwise_torsion_zero_does_not_establish_patch_solution(self):
        result = controls()['W19_pointwise_zero_is_not_patch_closure']
        self.assertEqual(result['torsion_coefficients'], [-1, 0, 1])

    def test_both_sector_coframes_pass_closure(self):
        result = controls()['W20_coframe_closure_does_not_select_sector']
        self.assertEqual(sm.diagonal_signature(result['K_metric']), (1, 2, 0))
        self.assertEqual(sm.diagonal_signature(result['R_metric']), (3, 0, 0))

    def test_linear_profile_is_flat_in_two_dimensions_but_not_three(self):
        flat = metric_curvature(2, 3, 1, 0, 1)
        self.assertEqual(flat['curvature']['0,1'], gm.scale(gm.identity(2), 0))
        curved = metric_curvature(3, 3, 1, 0, 1)
        self.assertEqual(curved['curvature']['1,2'],
                         gm.matrix(((0, 0, 0), (0, 0, 1), (0, -1, 0))))

    def test_generic_metric_oracle_also_detects_nonzero_second_derivative(self):
        old = warped_chart(Q(1, 2), (1, 0, 1))
        new = metric_curvature(2, Q(5, 4), 1, 2, 1)
        self.assertEqual(new['curvature']['0,1'], old['riemann_uv'])
        self.assertNotEqual(new['curvature']['0,1'], gm.scale(gm.identity(2), 0))

    def test_rotation_sector_reverses_spatial_curvature_sign(self):
        k = metric_curvature(3, 3, 2, 0, 1)['curvature']['1,2']
        r = metric_curvature(3, 3, 2, 0, -1)['curvature']['1,2']
        self.assertEqual(k, gm.scale(r, -1))

    def test_degenerate_adapter_or_zero_selected_coupling_is_rejected(self):
        for a, b, kappa in ((0, 3, 1), (1, 0, 1), (1, 3, 0)):
            with self.assertRaises(ValueError):
                sm.homogeneous_data(3, a, b, kappa, kappa, 1)
        with self.assertRaises(ValueError):
            metric_curvature(3, 0, 1, 0, 1)


if __name__ == '__main__':
    unittest.main()
