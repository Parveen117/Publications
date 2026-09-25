import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import grading_model as gr
import selection_model as sm
from grading_certificate import apply, controls


class NativeGradingTests(unittest.TestCase):
    def test_native_cut_and_metric_sign_involution_give_different_grading(self):
        native = sm.pair_generator(2, 0, 1, 1)
        sign = gm.matrix(((1, 0), (0, -1)))
        self.assertEqual(gm.mul(gm.mul(native, native), native), native)
        self.assertEqual(gm.mul(gm.mul(sign, native), sign), gm.scale(native, -1))

    def test_isolated_pair_has_native_common_grader(self):
        basis = gr.common_graders(2, ((0, 1),))
        self.assertEqual(len(basis), 1)
        self.assertEqual(gm.mul(basis[0], basis[0]), gm.identity(2))

    def test_disjoint_pairs_admit_a_common_involution(self):
        basis = gr.common_graders(4, ((0, 1), (2, 3)))
        self.assertEqual(len(basis), 2)
        j = gm.add(*basis)
        self.assertEqual(gm.mul(j, j), gm.identity(4))

    def test_overlapping_star_and_path_have_no_nonzero_common_grader(self):
        for edges in (((0, 1), (0, 2), (0, 3)), ((0, 1), (1, 2), (2, 3))):
            self.assertEqual(gr.common_graders(4, edges), ())

    def test_hidden_overlap_difference_changes_after_grading(self):
        result = controls()['W22_gluing_loses_grading']
        self.assertEqual(result['initial'], (0, 0, 0))
        self.assertEqual(result['after_J'], (0, 1, -1))

    def test_repaired_observer_is_minimal_and_intertwines_involution(self):
        data = gr.port_datum(4, ((0, 1), (0, 2), (0, 3)))
        repaired = gr.minimum_grading_observer(data['cut'], data['grading'])
        a, j = repaired['observer'], repaired['descended_grading']
        self.assertEqual(repaired['extra_scalar_channels'], 2)
        self.assertEqual(gm.mul(j, a), gm.mul(a, data['grading']))
        self.assertEqual(gm.mul(j, j), gm.identity(6))
        self.assertEqual(gr.rank(gr.grading_memory_lift(data['cut'], data['grading'])), 2)

    def test_cycle_parity_changes_grading_memory_rank(self):
        triangle = gr.port_datum(3, ((0, 1), (0, 2), (1, 2)))
        square = gr.port_datum(4, ((0, 1), (1, 2), (2, 3), (0, 3)))
        self.assertEqual(gr.rank(triangle['target']), 5)
        self.assertEqual(gr.rank(square['target']), 6)
        self.assertEqual(len(triangle['grading'])-gr.rank(triangle['target']), 1)
        self.assertEqual(len(square['grading'])-gr.rank(square['target']), 2)

    def test_grading_blind_triangle_state_is_visible_to_an_addressed_R_step(self):
        data = gr.port_datum(3, ((0, 1), (0, 2), (1, 2)))
        hidden = (1, -1, -1, 1, 1, -1)
        self.assertEqual(apply(data['target'], hidden), (0,)*6)
        for e in range(3):
            out = apply(gm.mul(data['cut'], gr.edge_step(3, e, Q(1, 2))), hidden)
            self.assertNotEqual(out, (0, 0, 0))

    def test_single_step_transcript_observes_all_endpoint_copies(self):
        data = gr.port_datum(4, ((0, 1), (1, 2), (2, 3), (0, 3)))
        steps = [gr.edge_step(4, e, Q(-2, 3)) for e in range(4)]
        self.assertEqual(gr.rank(data['target']), 6)
        self.assertEqual(gr.rank(gr.response_rows(data['cut'], steps)), 8)

    def test_edge_step_formula_matches_full_cayley_and_its_inverse(self):
        for sector, t in ((1, Q(1, 3)), (-1, Q(-2, 3))):
            actual = gr.edge_step(3, 1, t, sector)
            self.assertEqual(actual, gm.cayley(gr.edge_generator(3, 1, sector), t))
            self.assertEqual(gm.mul(actual, gr.edge_step(3, 1, -t, sector)), gm.identity(6))

    def test_declared_K_R_transcripts_are_separated(self):
        result = controls()['W25_declared_sector_transcripts']
        gap = max(abs(x-y) for x, y in zip(result['K_algebraic_candidate'], result['R_native_sector']))
        self.assertEqual(gap, Q(16, 15))

    def test_compression_does_not_licence_a_native_commutator(self):
        result = controls()['W26_compression_invents_commutator']
        self.assertTrue(result['native_commutes'])
        self.assertFalse(result['exact_descent'])
        self.assertNotEqual(result['compressed_commutator'], gm.scale(gm.identity(3), 0))

    def test_non_simple_graph_and_disconnected_port_model_are_rejected(self):
        for edges in (((0, 0),), ((0, 1), (0, 1)), ((1, 0),)):
            with self.assertRaises(ValueError):
                gr.graph_edges(3, edges)
        with self.assertRaisesRegex(ValueError, 'connected'):
            gr.port_datum(3, ((0, 1),))
        self.assertFalse(gr.graph_data(4, ((1, 2), (2, 3), (1, 3)))['bipartite'])

    def test_singular_or_trivial_protocol_cannot_claim_full_recovery(self):
        for sector, t in ((1, 1), (1, -1), (-1, 0)):
            with self.assertRaises(ValueError):
                gr.edge_step(2, 0, t, sector)
        with self.assertRaises(ValueError):
            gr.minimum_grading_observer(((1, 0),), ((2, 0), (0, 1)))


if __name__ == '__main__':
    unittest.main()
