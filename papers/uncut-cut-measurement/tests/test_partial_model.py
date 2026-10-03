import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import finite_model as fm
import family_model as fam
import partial_model as pm


class PartialContinuationTests(unittest.TestCase):
    def test_domain_mismatch_blocks_quotient_even_with_constant_readout(self):
        self.assertIsNone(pm.descended_partial((0, 0), (0, None)))
        self.assertEqual(pm.descended_partial((0, 1), (0, None)),
                         {0: (0,), 1: ()})

    def test_undefined_is_distinct_from_zero_and_none_readouts(self):
        for value in (0, None, (0,)):
            result = pm.partial_closure((value, value), ((0, None),))
            self.assertEqual(result['cut'], (0, 1))
        self.assertEqual(pm.descended_partial((None,), (0,)), {None: (None,)})

    def test_safe_indistinguishability_is_not_transitive(self):
        cut, arrows = (0, 0, 0, 1), ((0, None, 3, 3),)
        self.assertIsNone(pm.diagnostic(cut, arrows, 0, 1))
        self.assertIsNone(pm.diagnostic(cut, arrows, 1, 2))
        self.assertEqual(pm.diagnostic(cut, arrows, 0, 2)['word'], (0,))

    def test_availability_witness_is_not_an_executable_output_experiment(self):
        cut, arrows = (0, 0), ((0, None),)
        self.assertIsNone(pm.diagnostic(cut, arrows, 0, 1))
        witness = pm.diagnostic(cut, arrows, 0, 1, observe_availability=True)
        self.assertEqual(witness, {'kind': 'availability', 'word': (0,),
                                   'prefix': (), 'arrow': 0})
        self.assertEqual(pm.word_map(2, arrows, witness['word']), (0, None))

    def test_delayed_permission_difference_attains_depth_bound(self):
        cut, arrows = (0, 0, 0), ((1, 2, None),)
        result = pm.partial_closure(cut, arrows)
        self.assertEqual(result['cut'], (0, 1, 2))
        self.assertEqual(result['depth'], 2)
        self.assertEqual(fm.repair_size(cut, result['cut']), 3)
        witness = pm.diagnostic(cut, arrows, 0, 1, observe_availability=True)
        self.assertEqual(witness['kind'], 'availability')
        self.assertEqual(witness['prefix'], (0,))
        self.assertIsNone(pm.diagnostic(cut, arrows, 0, 1))

    def test_search_does_not_return_a_longer_availability_query_too_early(self):
        # After A, a domain difference can be queried; B already separates outputs.
        cut = (0, 0, 0, 0, 1)
        arrows = ((2, 3, None, 3, 4), (0, 4, 2, 3, 4))
        result = pm.diagnostic(cut, arrows, 0, 1, observe_availability=True)
        self.assertEqual(result['kind'], 'output')
        self.assertEqual(result['word'], (1,))

    def test_totalization_matches_domain_closure_without_readout_collision(self):
        cut, arrows = (None, None, (0,)), ((1, None, 2),)
        tc, ta = pm.totalize(cut, arrows)
        total = fam.family_closure(tc, ta)['cut']
        partial = pm.partial_closure(cut, arrows)['cut']
        self.assertEqual(fm.canonical(total[:-1]), partial)
        self.assertNotIn(total[-1], total[:-1])

    def test_total_empty_and_nowhere_defined_families(self):
        cut, arrows = (0, 0, 1), ((1, 2, 2),)
        self.assertEqual(pm.partial_closure(cut, arrows)['cut'],
                         fam.family_closure(cut, arrows)['cut'])
        self.assertEqual(pm.partial_closure(cut, ())['cut'], cut)
        self.assertEqual(pm.partial_closure(cut, ((None,) * 3,))['cut'], cut)

    def test_partial_composite_retains_exact_domain(self):
        arrows = ((1, None, 0), (None, 2, 2))
        self.assertEqual(pm.word_map(3, arrows, (0, 1)), (2, None, None))
        self.assertEqual(pm.word_map(3, arrows, ()), (0, 1, 2))

    def test_invalid_partial_arrows_and_all_word_labels_are_rejected(self):
        for arrows in (((3, None),), ((True, None),), ((-1, 0),), ((0,),)):
            with self.assertRaises(ValueError):
                pm.partial_closure((0, 0), arrows)
        with self.assertRaises(ValueError):
            pm.word_map(1, ((None,),), (0, 1))
        with self.assertRaises(ValueError):
            pm.diagnostic((0,), (), None, 0)


if __name__ == '__main__':
    unittest.main()
