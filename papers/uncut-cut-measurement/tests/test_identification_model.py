from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import identification_model as ident
import identification_certificate as cert
import response_model as rm


class IdentificationTests(unittest.TestCase):
    def test_reported_statistics_are_five_raw_readings(self):
        p = ident.protocol(*ident.pair(Q(1, 2), Q(1, 2)))
        self.assertEqual(len(p['raw']), 5)
        self.assertEqual(p['reported'], (Q(-8, 25), Q(16, 25), Q(3, 5)))
        a, b, c, d, e = p['raw']
        self.assertEqual(p['reported'], (a-b, c-d, e))

    def test_order_response_is_blind_to_reciprocal_mixing(self):
        u, v = ident.pair(Q(2, 3), Q(1, 3))
        ur, vr = ident.pair(Q(2, 3), Q(3))
        self.assertNotEqual(u, ur)
        self.assertEqual(rm.order_response(u, v), rm.order_response(ur, vr))

    def test_tensor_does_not_repair_reciprocal_blindness(self):
        a = rm.order_response(*ident.pair(Q(1, 2), Q(1, 2)))
        b = rm.order_response(*ident.pair(Q(1, 2), Q(2)))
        self.assertEqual(rm.tensor(a, rm.quarter_turn(2)), rm.tensor(b, rm.quarter_turn(2)))

    def test_signed_order_recovers_phase_and_mixing_sine(self):
        d, e, _ = ident.forward_formula(Q(-3, 2), Q(-4, 3))
        x = ident.order_invariants(d, e)
        self.assertEqual(x, {'phase': Q(-3, 2), 'mixing_sine': Q(-24, 25)})

    def test_direct_reading_repairs_both_reciprocal_branches(self):
        for s in (Q(1, 3), Q(3), Q(-1, 2), Q(-2), Q(1), Q(-1)):
            r = ident.decode_exact(ident.forward_formula(Q(2, 5), s))
            self.assertEqual((r['phase'], r['mixing']), (Q(2, 5), s))

    def test_arbitrary_triple_requires_circle_consistency(self):
        with self.assertRaises(ValueError):
            ident.decode_exact((1, 1, Q(1, 2)))
        with self.assertRaises(ValueError):
            ident.decode_exact((0, 1, 0))
        with self.assertRaises(ValueError):
            ident.decode_exact((1, 1, -1))

    def test_frame_change_preserves_independently_calibrated_readings(self):
        u, v = ident.pair(Q(-1, 3), Q(2))
        for f in cert.frames():
            fi = gm.inverse(f)
            up, vp = (gm.mul(gm.mul(f, x), fi) for x in (u, v))
            self.assertEqual(ident.protocol(up, vp, f), ident.protocol(u, v))

    def test_reciprocal_probes_fail_augmented_agreement(self):
        a, b = ident.forward_formula(Q(1, 2), Q(1, 2)), ident.forward_formula(Q(1, 2), Q(2))
        self.assertEqual(a[:2], b[:2])
        self.assertEqual(ident.compare_exact(a, b), 'SHARED_PAIR_REJECTED')

    def test_matching_probes_do_not_select_parameter_values(self):
        a = ident.forward_formula(Q(1, 2), Q(1, 2))
        b = ident.forward_formula(Q(1, 3), Q(2, 3))
        self.assertNotEqual(a, b)
        self.assertEqual(ident.compare_exact(a, a), 'EXACT_SHARED_PAIR')
        self.assertEqual(ident.compare_exact(b, b), 'EXACT_SHARED_PAIR')

    def test_posthoc_output_fit_can_erase_a_real_difference(self):
        a = rm.order_response(*ident.pair(Q(1, 2), Q(1, 2)))
        b = rm.order_response(*ident.pair(Q(1, 3), Q(2, 3)))
        self.assertNotEqual(a, b)
        self.assertEqual(gm.mul(gm.mul(a, gm.inverse(b)), b), a)

    def test_zero_error_box_is_exact_model_match(self):
        record = ident.forward_formula(Q(2, 3), Q(-3, 2))
        box = ident.decode_box(record, (0, 0, 0))
        self.assertEqual(box['status'], 'EXACT_MODEL_MATCH')
        self.assertEqual(box['phase'], (Q(2, 3),)*2)
        self.assertEqual(box['mixing'], (Q(-3, 2),)*2)

    def test_noisy_corner_still_contains_true_parameters(self):
        t, s, eps = Q(-1, 3), Q(2), Q(1, 10000)
        true = ident.forward_formula(t, s)
        measured = (true[0]+eps, true[1]-eps, true[2]+eps)
        box = ident.decode_box(measured, (eps,)*3)
        self.assertEqual(box['status'], 'CONDITIONAL_ENCLOSURE')
        self.assertTrue(ident.contains(box['phase'], t))
        self.assertTrue(ident.contains(box['mixing'], s))
        self.assertTrue(ident.candidate_compatible(t, s, measured, (eps,)*3))

    def test_circle_inconsistency_rejects_entire_error_box(self):
        self.assertEqual(ident.decode_box((1, 1, Q(1, 2)), (Q(1, 1000),)*3)['status'], 'MODEL_CLASS_REJECTED')

    def test_error_box_denominator_zero_does_not_create_pass(self):
        weak = ident.forward_formula(Q(1, 100), Q(1, 2))
        self.assertEqual(ident.decode_box(weak, (Q(1, 100),)*3)['status'], 'INSUFFICIENT_RESOLUTION')

    def test_large_mixing_requires_inverse_margin(self):
        record = ident.forward_formula(Q(1, 2), Q(1000))
        self.assertEqual(ident.decode_box(record, (Q(1, 100000),)*3)['status'], 'INSUFFICIENT_RESOLUTION')

    def test_overlap_remains_unresolved_and_narrow_boxes_reject(self):
        a, b = ident.forward_formula(Q(1, 2), Q(1, 2)), ident.forward_formula(Q(1, 2), Q(2))
        self.assertEqual(ident.compare_boxes(a, (0, 0, Q(1, 10)), b, (0, 0, Q(1, 10))), 'SHARED_PAIR_REJECTED')
        self.assertEqual(ident.compare_boxes(a, (0, 0, Q(3, 5)), b, (0, 0, Q(3, 5))), 'UNRESOLVED')

    def test_ratio_bound_controls_perturbed_phase(self):
        d, e, _ = ident.forward_formula(Q(3, 2), Q(2, 3))
        eps = Q(1, 100)
        dh, eh = d+eps, e-eps
        bound = ident.ratio_error_bound(dh, eh, eps, eps)
        self.assertLessEqual(abs(-dh/eh-Q(3, 2)), bound)
        with self.assertRaises(ValueError):
            ident.ratio_error_bound(1, Q(1, 10), 0, Q(1, 10))

    def test_interval_arithmetic_is_exact_and_handles_signs(self):
        self.assertEqual(ident.idiv((1, 2), (-4, -2)), (Q(-1), Q(-1, 4)))
        self.assertTrue(all(isinstance(x, Q) for x in ident.idiv((1, 2), (3, 4))))
        self.assertEqual(ident.isquare((Q(-3), Q(2))), (0, 9))
        with self.assertRaises(ValueError):
            ident.idiv((Q(1), Q(2)), (Q(-1), Q(1)))

    def test_invalid_budget_or_carrier_rejected(self):
        with self.assertRaises(ValueError):
            ident.decode_box((1, 1, 0), (-1, 0, 0))
        with self.assertRaises(ValueError):
            ident.decode_box((1, 1), (0, 0))
        with self.assertRaises(ValueError):
            ident.protocol(gm.identity(2), gm.identity(2))


if __name__ == '__main__':
    unittest.main()
