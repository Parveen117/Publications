import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import grading_model as gr
import interaction_model as im
from grading_certificate import apply
from gravity_certificate import in_span
from interaction_certificate import controls


class NativeInteractionTests(unittest.TestCase):
    def test_common_grading_keeps_mixing_even_and_local_R_odd(self):
        j = im.native_grading(3)
        b, r = im.mixing_generator(3, 0, 2), gr.edge_generator(3, 0, -1)
        self.assertEqual(gm.mul(gm.mul(j, b), j), b)
        self.assertEqual(gm.mul(gm.mul(j, r), j), gm.scale(r, -1))
        self.assertEqual(gm.transpose(b), gm.scale(b, -1))
        self.assertEqual(gm.transpose(r), gm.scale(r, -1))

    def test_closed_mixing_step_matches_generic_cayley(self):
        for q in (Q(-2, 3), Q(1), Q(3, 2)):
            u = im.mixing_step(3, 0, 2, q)
            self.assertEqual(u, gm.cayley(im.mixing_generator(3, 0, 2), q))
            self.assertEqual(gm.mul(u, im.mixing_step(3, 0, 2, -q)), gm.identity(6))

    def test_represented_quarter_turn_is_preserved_by_the_selected_family(self):
        z = im.kron(gm.identity(3), ((0, -1), (1, 0)))
        self.assertEqual(gm.mul(z, z), gm.scale(gm.identity(6), -1))
        for u in im.family(3, ((0, 1), (1, 2))):
            self.assertEqual(gm.mul(z, u), gm.mul(u, z))
        # It is not central in the entire EMK algebra: the cut reverses it.
        j = im.native_grading(3)
        self.assertEqual(gm.mul(gm.mul(j, z), j), gm.scale(z, -1))

    def test_genuine_ordered_readouts_differ_before_any_compression(self):
        result = controls()['W27_actual_order_response']
        self.assertEqual(result['R_then_mixing'], 0)
        self.assertEqual(result['mixing_then_R'], Q(16, 25))
        self.assertNotEqual(result['group_commutator'], gm.identity(4))

    def test_proper_lossy_cut_preserves_nonzero_interaction_and_grading(self):
        result = controls()['W28_lossy_interacting_quotient']
        self.assertEqual((result['source_dimension'], result['target_dimension']), (6, 4))
        self.assertTrue(result['grading_and_form_descend'])
        self.assertNotEqual(result['nonzero_target_commutator'], gm.scale(gm.identity(4), 0))

    def test_cut_inside_connected_component_fails_then_repair_is_full_rank(self):
        cut = im.projection(3, (0,))
        links = ((0, 1), (1, 2))
        self.assertIsNone(gm.descend_transport(im.mixing_step(3, 0, 1, Q(1, 2)), cut))
        repaired = im.row_closure(cut, im.family(3, links))
        self.assertEqual(len(repaired['observer']), 6)
        self.assertEqual(repaired['extra_scalar_channels'], 4)

    def test_unobserved_component_stays_out_of_the_closed_target(self):
        cut = (gm.identity(8)[0],)
        links = ((0, 1), (2, 3))
        predicted = im.predicted_observer(4, links, cut)
        actual = im.row_closure(cut, im.family(4, links))
        self.assertEqual(predicted['rank'], 4)
        self.assertEqual(actual['extra_scalar_channels'], 3)
        self.assertEqual(gr.rank(actual['observer']+predicted['observer']), 4)

    def test_cancellation_between_components_does_not_licence_discarding_them(self):
        result = controls()['W30_cancellation_is_not_blind_component']
        self.assertEqual(result['initial'], 0)
        self.assertEqual(result['after_local_R'], (Q(-2, 5),))
        self.assertEqual(result['closed_observer_rank'], 4)

    def test_invariant_form_scales_are_independent_only_between_components(self):
        for links, count in (((), 3), (((0, 1),), 2), (((0, 1), (1, 2)), 1)):
            basis = gm.invariant_forms(6, im.family(3, links))
            expected = im.invariant_form_basis(3, links)
            self.assertEqual(len(basis), count)
            self.assertTrue(all(in_span(h, basis) for h in expected))

    def test_optimal_transcript_decodes_for_nonstandard_parameters_and_roots(self):
        links = ((0, 1), (1, 2), (1, 3))
        state = tuple(Q((-1)**i*(i+1), i+2) for i in range(8))
        for root, phase, mixing in ((0, Q(-2, 3), Q(1)), (3, Q(1), Q(-3, 2))):
            packet = im.optimal_transcript(4, links, root, phase, mixing)
            values = apply(packet['rows'], state)
            self.assertEqual(im.decode_transcript(packet, values), state)
            self.assertEqual(apply(gm.inverse(packet['rows']), values), state)
            self.assertEqual(len(packet['rows']), 8)

    def test_chronological_protocol_order_matches_matrix_action(self):
        links = ((0, 1), (1, 2))
        packet = im.optimal_transcript(3, links, root=2)
        family = im.family(3, links)
        named = {'R0': family[0], 'R1': family[1], 'R2': family[2],
                 'B0,1': family[3], 'B1,2': family[4]}
        for row, word in zip(packet['rows'], packet['chronological_words']):
            whole = gm.identity(6)
            for letter in word:
                whole = gm.mul(named[letter], whole)
            self.assertEqual(whole[4], row)

    def test_single_channel_protocol_is_exactly_two_readings(self):
        packet = im.optimal_transcript(1, ())
        self.assertEqual(len(packet['rows']), 2)
        self.assertEqual(packet['chronological_words'], ((), ('R0',)))
        self.assertEqual(im.decode_transcript(packet, apply(packet['rows'], (2, -3))), (2, -3))

    def test_missing_control_changes_observability_target(self):
        result = controls()['W32_missing_control']
        self.assertEqual(result['mixing_only_rank'], 2)
        self.assertEqual(result['local_phases_only_rank'], 2)
        self.assertEqual(result['full_connected_family_rank'], 4)

    def test_symmetric_mixing_would_violate_native_energy_rule(self):
        wrong = im.kron(((0, 1), (1, 0)), gm.identity(2))
        step = gm.cayley(wrong, Q(1, 2))
        self.assertNotEqual(gm.mul(gm.transpose(step), step), gm.identity(4))

    def test_invalid_graphs_zero_controls_and_disconnected_full_protocol_rejected(self):
        for links in (((0, 0),), ((0, 1), (0, 1)), ((1, 0),)):
            with self.assertRaises(ValueError):
                im.links_for(2, links)
        with self.assertRaises(ValueError):
            im.mixing_step(2, 0, 1, 0)
        with self.assertRaisesRegex(ValueError, 'connected'):
            im.optimal_transcript(2, ())

    def test_readout_and_decoder_dimensions_are_enforced(self):
        with self.assertRaises(ValueError):
            im.predicted_observer(2, ((0, 1),), ((0, 0, 0, 0),))
        with self.assertRaises(ValueError):
            im.row_closure(((1, 0),), (gm.identity(4),))
        with self.assertRaises(ValueError):
            im.decode_transcript(im.optimal_transcript(1, ()), (1,))


if __name__ == '__main__':
    unittest.main()
