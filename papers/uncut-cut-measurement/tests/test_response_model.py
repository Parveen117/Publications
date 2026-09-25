from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
import grading_model as gr
import interaction_model as im
import response_model as rm
import response_certificate as rc


class ResponseTests(unittest.TestCase):
    def test_complex_gram_matches_direct_pair_arithmetic(self):
        y = gm.matrix(((1, 2, 0), (Q(1, 2), -1, 3), (0, 2, -2), (4, 0, 1)))
        self.assertEqual(rm.tensor(y, rm.quarter_turn(2)), rc.paired_gram_oracle(y))

    def test_gram_positive_on_gaussian_rational_coefficients(self):
        y = gm.matrix(((1, 2), (3, -1), (-2, Q(1, 3)), (0, 4)))
        g, a = rm.tensor(y, rm.quarter_turn(2))
        # v = p+i q: v*Q v = p^T Gp+q^T Gq-2p^T A q.
        p, q = (Q(2), Q(-1)), (Q(1, 3), Q(2))
        quad = sum(p[i]*g[i][j]*p[j]+q[i]*g[i][j]*q[j]-2*p[i]*a[i][j]*q[j]
                   for i in range(2) for j in range(2))
        expected = sum(sum(y[2*b][j]*p[j]-y[2*b+1][j]*q[j] for j in range(2))**2
                     + sum(y[2*b][j]*q[j]+y[2*b+1][j]*p[j] for j in range(2))**2
                       for b in range(2))
        self.assertEqual(quad, expected)
        self.assertGreaterEqual(quad, 0)

    def test_nonorthogonal_frame_requires_metric_transport(self):
        y, z = gm.identity(2), rm.quarter_turn(1)
        f = gm.matrix(((2, 1), (0, 3)))
        fi = gm.inverse(f)
        h, znew = gm.mul(gm.transpose(fi), fi), gm.mul(gm.mul(f, z), fi)
        self.assertEqual(rm.tensor(f, znew, h), rm.tensor(y, z))
        self.assertNotEqual(gm.mul(gm.transpose(f), f), gm.identity(2))

    def test_grading_conjugates_pairing(self):
        z, y = rm.quarter_turn(2), gm.identity(4)
        g, a = rm.tensor(y, z)
        self.assertEqual(rm.tensor(im.native_grading(2), z), (g, gm.scale(a, -1)))

    def test_half_channel_cut_needs_imaginary_seam(self):
        d = rm.split_tensor(gm.identity(2), ((1, 0),), rm.quarter_turn(1))
        self.assertEqual(d['cross'][0], rm.zero(2))
        self.assertNotEqual(d['cross'][1], rm.zero(2))
        self.assertEqual(d['full'], rm.tensor_add(rm.tensor_add(d['visible'], d['hidden']), d['cross']))

    def test_real_projection_is_not_quantum_monotonicity(self):
        d = rm.split_tensor(gm.identity(2), ((1, 0),), rm.quarter_turn(1))
        g, a = rm.tensor_sub(d['full'], d['visible'])
        self.assertEqual(g[0][0]*g[1][1]-g[0][1]**2-a[0][1]**2, -1)

    def test_one_catalogue_can_miss_an_incompatible_cut(self):
        d = rm.split_tensor(((1,), (0,)), ((1, 0),), rm.quarter_turn(1))
        self.assertEqual(d['cross'][1], ((Q(0),),))
        self.assertFalse(d['global_ledger_closes'])

    def test_noncoordinate_pair_repair_is_minimal_and_closed(self):
        c = gm.matrix(((1, 0, 1, 0),))
        z = rm.quarter_turn(2)
        d = rm.paired_repair(c, z)
        self.assertEqual(d['extra_scalar_channels'], 1)
        self.assertTrue(rm.split_tensor(gm.identity(4), d['observer'], z)['global_ledger_closes'])
        self.assertEqual(gr.rank(d['observer']+gm.mul(d['observer'], z)), 2)

    def test_redundant_and_zero_readouts(self):
        z = rm.quarter_turn(1)
        self.assertEqual(rm.paired_repair(((1, 0), (2, 0)), z)['extra_scalar_channels'], 1)
        self.assertEqual(rm.paired_repair(((0, 0),), z)['rank'], 0)
        self.assertEqual(rm.projector(((0, 0),)), rm.zero(2))

    def test_tensor_closure_is_weaker_than_interaction_closure(self):
        c, z = im.projection(2, (0,)), rm.quarter_turn(2)
        self.assertTrue(rm.split_tensor(gm.identity(4), c, z)['global_ledger_closes'])
        self.assertIsNone(gm.descend_transport(im.mixing_step(2, 0, 1, Q(1, 2)), c))

    def test_signed_response_and_loop_input_conventions(self):
        u, v, ell = rc.loop_datum(2, 0, 1, Q(1, 2), Q(1, 2))
        self.assertEqual(rc.apply(ell, (0, 0, 0, 1)), (Q(-16, 25), Q(-8, 25), 0, 0))
        loop = gm.sub(gm.mul(gm.mul(gm.mul(u, v), gm.transpose(u)), gm.transpose(v)), gm.identity(4))
        self.assertEqual(loop, gm.mul(gm.mul(ell, gm.transpose(u)), gm.transpose(v)))
        self.assertNotEqual(loop, ell)

    def test_response_scale_is_positive_and_target_rank_four(self):
        for phase, mixing in ((Q(1, 2), Q(1, 2)), (Q(-1), Q(2, 3)), (Q(3), Q(-2))):
            _, _, ell = rc.loop_datum(3, 0, 2, phase, mixing)
            c = im.projection(3, (0, 2))
            k = 16*phase**2*mixing**2/((1+phase**2)*(1+mixing**2)**2)
            self.assertEqual(gm.mul(gm.transpose(ell), ell), gm.scale(gm.mul(gm.transpose(c), c), k))
            self.assertEqual(rm.target_repair((gm.identity(6)[0],), ell)['extra_scalar_channels'], 3)

    def test_full_tensor_loses_overall_response_sign(self):
        _, _, ell = rc.loop_datum(2, 0, 1, Q(1, 2), Q(1, 2))
        self.assertNotEqual(ell, gm.scale(ell, -1))
        self.assertEqual(rm.tensor(ell, rm.quarter_turn(2)), rm.tensor(gm.scale(ell, -1), rm.quarter_turn(2)))

    def test_zero_mean_does_not_mean_zero_information(self):
        d = rm.conditional_ledger(((1, -2), (-1, 2)), (Q(1, 2), Q(1, 2)), (0, 0))
        self.assertEqual(d['mean'], (0, 0))
        self.assertEqual(d['discarded'], gm.matrix(((1, -2), (-2, 4))))

    def test_covariance_ledger_matches_independent_pairwise_formula(self):
        points = ((1, 0), (2, 3), (-1, 2), (-1, 2))
        weights = (Q(1, 10), Q(2, 10), Q(3, 10), Q(4, 10))
        labels = (0, 0, 1, 1)
        d = rm.conditional_ledger(points, weights, labels)
        self.assertEqual(d['discarded'], rc.pairwise_discard(points, weights, labels))
        self.assertEqual(d['total'], gm.add(d['recognized'], d['discarded']))

    def test_zero_discard_is_target_faithfulness_on_support(self):
        d = rm.conditional_ledger(((1, 2), (1, 2), (3, 4)), (Q(1, 3),)*3, (0, 0, 1))
        self.assertEqual(d['discarded'], rm.zero(2))
        bad = rm.conditional_ledger(((1, 2), (1, 2), (3, 4)), (Q(1, 3),)*3, (0, 0, 0))
        self.assertNotEqual(bad['discarded'], rm.zero(2))

    def test_nested_variance_tower(self):
        points, p = ((0,), (1,), (3,), (4,)), (Q(1, 4),)*4
        fine = rm.conditional_ledger(points, p, (0, 0, 1, 2))
        coarse = rm.conditional_ledger(points, p, (0, 0, 1, 1))
        extra = rm.conditional_ledger(fine['block_means'], fine['block_weights'], (0, 1, 1))
        self.assertEqual(coarse['discarded'], gm.add(fine['discarded'], extra['discarded']))

    def test_invalid_ensemble_and_incompatible_quarter_turn_rejected(self):
        with self.assertRaises(ValueError):
            rm.mean_covariance(((1,), (2,)), (1, 0))
        with self.assertRaises(ValueError):
            rm.mean_covariance(((1,), (2,)), (1, 1))
        with self.assertRaises(ValueError):
            rm.tensor(gm.identity(2), gm.identity(2))
        with self.assertRaises(ValueError):
            rm.paired_repair(((1, 0),), gm.identity(2))


if __name__ == '__main__':
    unittest.main()
