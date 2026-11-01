"""Independent small refusal/contract tests; the CLI rebuilds the full carrier."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "certificates"))
import ym42_silence_time_ladder as y


class YM42Tests(unittest.TestCase):
    def test_exact_two_channel_model(self):
        B = [[F(1), F(0)], [F(0), F(1)]]
        M = [[F(1), F(1,10)], [F(1,10), F(1,8)]]
        P, Q = y.silence_restriction(M, B, 0)
        self.assertEqual(P, [[F(1,8)]])
        self.assertEqual(Q, [[F(1)]])
        self.assertTrue(y.ceiling_checks(P, Q, F(1,8))[0])
        budget = y.structural_budget(F(1), F(1), F(1,100), F(1,8))
        self.assertEqual(budget["beta"], F(29,196))
        self.assertTrue(budget["accepted"])
        count = y.finite_iterate_controls(M, B, 0, F(1), [F(0), F(1,10)],
                                         F(1,100), F(1), budget, horizon=8)
        self.assertEqual(count, 76)

    def test_diagonal_bounds_do_not_bound_a_complement(self):
        B = [[F(1), F(0)], [F(0), F(1)]]
        M = [[F(1,2), F(1,2)], [F(1,2), F(1,2)]]
        self.assertFalse(y.ceiling_checks(M, B, F(1,2))[0])

    def test_both_signs_are_needed(self):
        self.assertFalse(y.ceiling_checks([[F(-2)]], [[F(1)]], F(1))[0])
        self.assertTrue(y.ceiling_checks([[F(-2)]], [[F(1)]], F(2))[0])

    def test_all_negative_controls_fire(self):
        refused = y.negative_controls()
        self.assertEqual(len(refused), 8)
        self.assertTrue(all(refused.values()))

    def test_space_contraction_does_not_select_time_ratio(self):
        rows = y.time_blindness_control()["finite_positive_models"]
        self.assertEqual({r["beta_space"] for r in rows}, {"0"})
        self.assertEqual(len({r["time_ratio"] for r in rows}), 4)
        self.assertEqual(rows[-1]["time_ratio"], "1023/1024")

    def test_directed_decimals_are_upper_bounds(self):
        for value in (F(1,3), F(-1,3), F(0), F(1), F(17,8)):
            shown = F(y.upper_decimal(value, places=6))
            self.assertTrue(value <= shown < value + F(1,10**6))

    def test_source_tamper_fails(self):
        with patch.object(y, "digest", return_value="tampered"):
            with self.assertRaisesRegex(ValueError, "upstream source pin changed"):
                y.verify_sources()

    def test_result_and_pin_tamper_fail_independently(self):
        cert = {"verdict":"PASS", "value":"1/3"}
        with tempfile.TemporaryDirectory() as folder:
            result, pin = Path(folder)/"result.json", Path(folder)/"pin.txt"
            result.write_text(json.dumps(cert))
            pin.write_text(y.canonical_sha(cert)+"\n")
            y.check_evidence(cert, result, pin)
            result.write_text(json.dumps({**cert, "value":"1/2"}))
            with self.assertRaises(ValueError):
                y.check_evidence(cert, result, pin)
            result.write_text(json.dumps(cert))
            pin.write_text("0"*64)
            with self.assertRaises(ValueError):
                y.check_evidence(cert, result, pin)


if __name__ == "__main__":
    unittest.main()
