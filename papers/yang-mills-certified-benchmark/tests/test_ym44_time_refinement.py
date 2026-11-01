from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym44_time_refinement as y


class TimeRefinementTests(unittest.TestCase):
    def test_character_moment_endpoints_and_full_grid(self):
        self.assertEqual(y.moment_controls()['moment_enclosures'],136)

    def test_ordered_cubic_residue_and_half_steps(self):
        out=y.ordered_word_controls()
        self.assertEqual(out['word_identities'],30)
        self.assertEqual(out['cubic_residue']['BAB'],'-1/6')
        self.assertEqual(y.split_word_coefficient('BAB'),0)
        self.assertEqual(y.split_word_coefficient('ABA'),F(1,4))

    def test_factorial_tail_and_invalid_gate(self):
        for x,N in ((F(-1),4),(F(4),2),(F(1),-1)):
            with self.assertRaises(ValueError):
                y.factorial_tail(x,N)
        self.assertGreater(y.factorial_tail(F(1),8),0)
        self.assertLess(y.factorial_tail(F(1),12),y.factorial_tail(F(1),8))
        self.assertEqual(y.factorial_tail(F(0),8),0)

    def test_width_sign_and_free_budget(self):
        self.assertEqual(y.refinement_budget(1,F(4),F(1),F(1,4)),0)
        self.assertEqual(y.refinement_budget(5,F(0),F(1),F(1,4)),0)
        self.assertEqual(y.refinement_budget(3,F(-1,7),F(1),F(1,9)),
                         y.refinement_budget(3,F(1,7),F(1),F(1,9)))
        self.assertGreater(y.refinement_budget(4,F(1,16),F(1),F(1,16)),
                           y.refinement_budget(2,F(1,16),F(1),F(1,16)))
        with self.assertRaises(ValueError):
            y.refinement_budget(0,F(1),F(1),F(1,2))

    def test_strict_normalized_gap_gate(self):
        self.assertIsNone(y.gap_gate(F(2),F(1,2),F(1,2)))
        self.assertIsNone(y.gap_gate(F(1),F(1),F(0)))
        self.assertIsNone(y.gap_gate(F(0),F(1,2),F(0)))
        self.assertIsNone(y.gap_gate(F(1),F(1,2),F(-1)))
        self.assertEqual(y.gap_gate(F(1),F(1,2),F(1,8)),F(5,7))
        self.assertEqual(y.gap_gate(F(8),F(1,2),F(1)),F(5,7))

    def test_independent_non_grid_refinement(self):
        t,theta=F(3,4),F(1,7)
        U=y.reference2(t,theta);other=y.taylor_reference2(t,theta)
        S=y.ivpow2(y.split2(t/25,theta),25)
        for i in range(2):
            for j in range(2):
                self.assertFalse(U[i][j].separated_from(other[i][j]))
        residual=[[S[i][j]-U[i][j] for j in range(2)] for i in range(2)]
        self.assertLess(y.row_bound2(residual),y.refinement_budget(2,theta,t,t/25))

    def test_interaction_prevents_exact_subdivision(self):
        a=y.split2(F(1),F(1,16))
        b=y.ivpow2(y.split2(F(1,2),F(1,16)),2)
        self.assertTrue(a[0][0].separated_from(b[0][0]))

    def test_refusals_and_slow_spectator(self):
        controls=y.negative_controls()['groups']
        self.assertEqual(len(controls),10)
        self.assertTrue(all(controls.values()))

    def test_upstream_source_tamper(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()

    def test_record_and_pin_tamper(self):
        data={'verdict':'PASS','runtime':'3.12'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'r.json',Path(folder)/'pin'
            record.write_text(json.dumps(data));pin.write_text(y.canonical_sha(data))
            y.check(data,record,pin)
            record.write_text(json.dumps({**data,'verdict':'FAIL'}))
            with self.assertRaises(ValueError):
                y.check(data,record,pin)
            record.write_text(json.dumps(data));pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(data,record,pin)


if __name__=='__main__':
    unittest.main()
