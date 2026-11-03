import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import finite_model as fm


class CutMeasurementTests(unittest.TestCase):
    def test_readout_alone_cannot_recover_hidden_bit(self):
        cut, target = (0, 0, 1, 1), (0, 1, 0, 1)
        self.assertIsNone(fm.factor_map(cut, target))
        self.assertEqual(fm.target_sets(cut, target),
                         {0: frozenset((0, 1)), 1: frozenset((0, 1))})
        self.assertEqual(fm.repair_size(cut, target), 2)

    def test_existing_target_needs_no_additional_bit(self):
        cut = (0, 0, 1, 1)
        self.assertEqual(fm.repair_size(cut, cut), 1)
        self.assertEqual(fm.fixed_binary_bits(1), 0)

    def test_memory_labels_can_be_reused_across_fibres(self):
        cut, target = (0, 0, 1, 1), (4, 5, 7, 8)
        memory = fm.repair_channel(cut, target)
        self.assertEqual(memory, (0, 1, 0, 1))
        decoder = fm.factor_map(tuple(zip(cut, memory)), target)
        self.assertEqual(decoder, {(0, 0): 4, (0, 1): 5,
                                   (1, 0): 7, (1, 1): 8})

    def test_refinement_product_is_only_an_upper_bound(self):
        coarse, fine, target = (0, 0, 0, 0), (0, 0, 1, 1), (0, 1, 0, 1)
        self.assertTrue(fm.refines(fine, coarse))
        self.assertEqual(fm.repair_size(coarse, target), 2)
        self.assertEqual(fm.repair_size(coarse, fine) *
                         fm.repair_size(fine, target), 4)

    def test_pairwise_consistency_can_fail_globally(self):
        states = ((0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0))
        cuts = tuple(tuple(x[j] for x in states) for j in range(3))
        for i, j in ((0, 1), (0, 2), (1, 2)):
            self.assertTrue(fm.joint_fibre((cuts[i], cuts[j]), (1, 1)))
        self.assertEqual(fm.joint_fibre(cuts, (1, 1, 1)), ())

    def test_separate_common_targets_and_joint_information_differ(self):
        cuts = ((0, 0, 1, 1), (0, 1, 0, 1))
        self.assertEqual(fm.common_components(cuts), (0, 0, 0, 0))
        self.assertEqual(len(set(zip(*cuts))), 4)

    def test_delayed_future_separates_apparently_equal_states(self):
        cut, transition = (0, 0, 0, 1), (1, 2, 3, 3)
        result = fm.stable_history(cut, transition)
        self.assertEqual(tuple(len(set(x)) for x in result["levels"]), (2, 3, 4))
        self.assertEqual(result["depth"], 2)
        self.assertEqual(result["cut"], (0, 1, 2, 3))
        self.assertEqual(fm.repair_size(cut, result["cut"]), 3)
        self.assertEqual(fm.fixed_binary_bits(3), 2)

    def test_constant_cut_already_closes_but_does_not_recover_state(self):
        result = fm.stable_history((0, 0, 0), (1, 2, 0))
        self.assertEqual(result["depth"], 0)
        self.assertEqual(result["cut"], (0, 0, 0))
        self.assertIsNone(fm.factor_map(result["cut"], (0, 1, 2)))

    def test_one_step_descent_can_fail(self):
        self.assertIsNone(fm.descended_transition((0, 0, 1), (0, 2, 2)))

    def test_distinct_states_can_have_all_future_readouts_equal(self):
        result = fm.stable_history((0, 0, 1, 1), (1, 0, 3, 2))
        self.assertEqual(result["cut"], (0, 0, 1, 1))
        self.assertEqual(result["transition"], {0: 0, 1: 1})

    def test_singleton_and_noncanonical_labels(self):
        self.assertEqual(fm.stable_history(("only",), (0,))["depth"], 0)
        self.assertEqual(fm.factor_map((8, 8, 3), (5, 5, 6)), {8: 5, 3: 6})
        self.assertEqual(fm.canonical(("b", "b", "a")), (0, 0, 1))

    def test_invalid_carriers_and_transitions_are_rejected(self):
        with self.assertRaises(ValueError):
            fm.factor_map((), ())
        with self.assertRaises(ValueError):
            fm.factor_map((0, 1), (0,))
        with self.assertRaises(ValueError):
            fm.stable_history((0, 1), (1, 2))
        with self.assertRaises(ValueError):
            fm.stable_history((0, 1), (1, -1))
        with self.assertRaises(ValueError):
            fm.joint_fibre(((0, 1),), ())

    def test_partition_generator_has_known_small_counts(self):
        self.assertEqual([len(tuple(fm.partitions(n))) for n in range(1, 5)],
                         [1, 2, 5, 15])


if __name__ == "__main__":
    unittest.main()
