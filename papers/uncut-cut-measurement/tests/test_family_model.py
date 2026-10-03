import sys
import unittest
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import finite_model as fm
import family_model as fam

C = (0, 0, 0, 1)
A = (0, 2, 2, 3)
B = (0, 0, 3, 3)


class FamilyContinuationTests(unittest.TestCase):
    def test_empty_family_retains_only_initial_distinctions(self):
        result = fam.family_closure((8, 8, 2), ())
        self.assertEqual(result['cut'], (0, 0, 1))
        self.assertEqual(result['depth'], 0)
        self.assertEqual(result['transitions'], ())
        self.assertIsNone(fam.distinguishing_word((0, 0, 1), (), 0, 1))
        self.assertEqual(fam.distinguishing_word((0, 0, 1), (), 0, 2), ())

    def test_one_arrow_agrees_with_original_history_result(self):
        transition = (1, 2, 3, 3)
        original = fm.stable_history(C, transition)
        family = fam.family_closure(C, (transition,))
        self.assertEqual(family['cut'], original['cut'])
        self.assertEqual(family['depth'], original['depth'])

    def test_separate_closures_do_not_capture_mixed_words(self):
        ha = fam.family_closure(C, (A,))['cut']
        hb = fam.family_closure(C, (B,))['cut']
        joint = fam.joint_cut(ha, hb)
        self.assertEqual(ha, (0, 0, 0, 1))
        self.assertEqual(hb, (0, 0, 1, 2))
        self.assertEqual(joint, hb)
        self.assertTrue(fam.closed_for_family(ha, (A,)))
        self.assertTrue(fam.closed_for_family(hb, (B,)))
        self.assertFalse(fam.closed_for_family(joint, (A, B)))

    def test_family_requires_three_labels_not_product_of_single_banks(self):
        closure = fam.family_closure(C, (A, B))
        self.assertEqual(closure['cut'], (0, 1, 2, 3))
        self.assertEqual([len(set(p)) for p in closure['levels']], [2, 3, 4])
        self.assertEqual(fm.repair_size(C, closure['cut']), 3)
        self.assertEqual(fm.fixed_binary_bits(3), 2)
        for memory in product(range(2), repeat=4):
            self.assertFalse(fam.closed_for_family(fam.joint_cut(C, memory), (A, B)))

    def test_shortest_diagnostic_word_detects_order(self):
        self.assertEqual(fam.distinguishing_word(C, (A, B), 0, 1), (0, 1))
        ab = fam.word_map(4, (A, B), (0, 1))
        ba = fam.word_map(4, (A, B), (1, 0))
        self.assertEqual((C[ab[0]], C[ab[1]]), (0, 1))
        self.assertEqual((C[ba[0]], C[ba[1]]), (0, 0))

    def test_each_unmixed_word_fails_to_separate_hidden_pair(self):
        # Search exhausts the pair graph; this is not a fixed-depth sample.
        self.assertIsNone(fam.distinguishing_word(C, (A,), 0, 1))
        self.assertIsNone(fam.distinguishing_word(C, (B,), 0, 1))

    def test_composite_shortens_depth_without_changing_memory(self):
        composite = fam.word_map(4, (A, B), (0, 1))
        before = fam.family_closure(C, (A, B))
        after = fam.family_closure(C, (A, B, composite))
        self.assertEqual(before['cut'], after['cut'])
        self.assertEqual(before['depth'], 2)
        self.assertEqual(after['depth'], 1)
        self.assertEqual(fm.repair_size(C, before['cut']),
                         fm.repair_size(C, after['cut']))

    def test_duplicate_identity_and_generator_order_preserve_closure(self):
        identity = tuple(range(4))
        base = fam.family_closure(C, (A, B))['cut']
        self.assertEqual(fam.family_closure(C, (B, A, A, identity))['cut'], base)

    def test_joint_initial_cuts_close_consistently_for_same_family(self):
        d = (0, 1, 0, 1)
        left = fam.family_closure(fam.joint_cut(C, d), (A, B))['cut']
        right = fam.joint_cut(fam.family_closure(C, (A, B))['cut'],
                             fam.family_closure(d, (A, B))['cut'])
        self.assertEqual(left, right)

    def test_finer_initial_observer_can_require_more_future_memory(self):
        coarse, fine, transition = (0, 0, 0), (0, 0, 1), (0, 2, 2)
        hc = fam.family_closure(coarse, (transition,))['cut']
        hd = fam.family_closure(fine, (transition,))['cut']
        self.assertTrue(fm.refines(fine, coarse))
        self.assertEqual(fm.repair_size(coarse, hc), 1)
        self.assertEqual(fm.repair_size(fine, hd), 2)

    def test_observer_and_arrow_relabelling_preserve_partition(self):
        rename = (2, 0, 3, 1)
        new_cut = [None] * 4
        new_arrows = [[None] * 4 for _ in (A, B)]
        for x in range(4):
            new_cut[rename[x]] = C[x]
            for j, transition in enumerate((A, B)):
                new_arrows[j][rename[x]] = rename[transition[x]]
        original = fam.family_closure(C, (A, B))['cut']
        changed = fam.family_closure(new_cut, new_arrows)['cut']
        pulled_back = fm.canonical(tuple(changed[rename[x]] for x in range(4)))
        self.assertEqual(original, pulled_back)
        self.assertEqual(fm.repair_size(C, original), fm.repair_size(new_cut, changed))

    def test_shortest_word_bound_and_initially_distinct_candidates(self):
        transition = (1, 2, 3, 3)
        self.assertEqual(fam.distinguishing_word(C, (transition,), 0, 1), (0, 0))
        self.assertEqual(fam.distinguishing_word(C, (transition,), 0, 3), ())
        self.assertIsNone(fam.distinguishing_word(C, (transition,), 1, 1))

    def test_every_generated_transition_has_a_closed_observed_map(self):
        cut = (0, 0, 1, 1)
        transitions = ((1, 0, 3, 2), (2, 3, 0, 1))
        result = fam.family_closure(cut, transitions)
        self.assertEqual(result['cut'], cut)
        self.assertEqual(result['transitions'], ({0: 0, 1: 1}, {0: 1, 1: 0}))
        word = (0, 1, 0, 1)
        actual = fam.word_map(4, transitions, word)
        for x in range(4):
            observed = result['cut'][x]
            for a in word:
                observed = result['transitions'][a][observed]
            self.assertEqual(observed, result['cut'][actual[x]])

    def test_invalid_family_and_words_are_rejected(self):
        with self.assertRaises(ValueError):
            fam.family_closure((), ())
        with self.assertRaises(ValueError):
            fam.family_closure((0, 1), ((1,),))
        with self.assertRaises(ValueError):
            fam.family_closure((0, 1), ((1, 2),))
        with self.assertRaises(ValueError):
            fam.word_map(2, ((1, 0),), (1,))
        with self.assertRaises(ValueError):
            fam.word_map(2, (), (0,))
        with self.assertRaises(ValueError):
            fam.distinguishing_word((0, 1), (), 0, 2)
        with self.assertRaises(ValueError):
            fam.joint_cut((0, 1), (0,))


if __name__ == '__main__':
    unittest.main()
