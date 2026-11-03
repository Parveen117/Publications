import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import gravity_model as gm
from gravity_certificate import in_span, lorentz, warped_chart


class TransportGeometryTests(unittest.TestCase):
    def test_full_star_selects_one_lorentz_type_form_without_metric_input(self):
        arrows = tuple(gm.emk_step(4, i, Q(i, i + 1)) for i in range(1, 4))
        forms = gm.invariant_forms(4, arrows)
        self.assertEqual(len(forms), 1)
        self.assertTrue(in_span(lorentz(4), forms))

    def test_rotation_sector_selects_a_different_signature(self):
        arrows = tuple(gm.emk_step(4, i, Q(1, 2), rotation=True) for i in range(1, 4))
        forms = gm.invariant_forms(4, arrows)
        self.assertEqual(len(forms), 1)
        self.assertTrue(in_span(gm.identity(4), forms))
        self.assertFalse(in_span(lorentz(4), forms))

    def test_K_candidate_is_not_positive_metric_unitary(self):
        u = gm.emk_step(2, 1, Q(1, 2))
        self.assertNotEqual(gm.congruence(u, gm.identity(2)), gm.identity(2))
        self.assertEqual(gm.congruence(u, lorentz(2)), lorentz(2))

    def test_mixing_sectors_can_destroy_all_compatible_forms(self):
        self.assertEqual(gm.invariant_forms(2, (gm.emk_step(2, 1, Q(1, 2)),
                         gm.emk_step(2, 1, Q(1, 2), rotation=True))), ())

    def test_missing_couplings_leave_geometry_underdetermined(self):
        self.assertEqual(len(gm.invariant_forms(4, (gm.emk_step(4, 1, Q(1, 2)),))), 4)
        self.assertEqual(len(gm.invariant_forms(3, ())), 6)

    def test_dilation_is_a_nonmetric_holonomy_obstruction(self):
        self.assertEqual(gm.invariant_forms(2, (gm.scale(gm.identity(2), 2),)), ())

    def test_root_form_transport_uses_the_inverse_path(self):
        path = gm.matrix(((1, 2), (0, 3)))
        source = lorentz(2)
        target = gm.metric_at_root(source, path)
        self.assertEqual(gm.congruence(path, target), source)
        self.assertNotEqual(target, gm.congruence(path, source))

    def test_lossy_cut_cannot_exactly_carry_a_nondegenerate_form(self):
        cut = gm.matrix(((1, 0),))
        self.assertIsNone(gm.descend_form(lorentz(2), cut))
        self.assertEqual(gm.descend_form(((2, 0), (0, 0)), cut), gm.matrix(((2,),)))
        self.assertEqual(gm.descend_form(lorentz(2), gm.identity(2)), lorentz(2))
        # A memory-extended response target can be degenerate before the cut.
        b = gm.emk_step(2, 1, Q(1, 2))
        lift = gm.matrix((b[0] + (0,), b[1] + (0,), (1, 0, 2)))
        memory = gm.matrix(((1, 0, 0), (0, 1, 0), (0, 0, 2)))
        source = gm.matrix(((1, 0, 0), (0, -1, 0), (0, 0, 0)))
        c = gm.matrix(((1, 0, 0), (0, 1, 0)))
        forms = gm.invariant_forms(3, (lift, memory))
        self.assertEqual(len(forms), 1)
        self.assertTrue(in_span(source, forms))
        self.assertEqual(gm.descend_form(source, c), lorentz(2))
        self.assertEqual(gm.descend_transport(lift, c), b)

    def test_compression_artifact_is_rejected_by_exact_descent(self):
        a = gm.matrix(((0, 0, 1), (1, 0, 0), (0, 1, 0)))
        b = gm.mul(a, a)
        cut = gm.matrix(((1, 0, 0), (0, 1, 0)))
        p = gm.transpose(cut)
        self.assertEqual(gm.mul(a, b), gm.mul(b, a))
        ac, bc = gm.mul(gm.mul(cut, a), p), gm.mul(gm.mul(cut, b), p)
        self.assertNotEqual(gm.mul(ac, bc), gm.mul(bc, ac))
        self.assertIsNone(gm.descend_transport(a, cut))
        self.assertEqual(gm.descend_transport(gm.identity(3), cut), gm.identity(2))

    def test_supplied_torsion_free_chart_has_two_route_curvature_agreement(self):
        curved = warped_chart(Q(1, 2), (1, 0, 1))
        self.assertEqual(curved['torsion'], (0, 0))
        self.assertEqual(curved['riemann_uv'], gm.matrix(((0, Q(5, 2)), (Q(8, 5), 0))))
        flat = warped_chart(0, (1, 0, 0))
        self.assertEqual(flat['riemann_uv'], gm.matrix(((0, 0), (0, 0))))

    def test_metric_compatibility_without_torsion_gate_is_insufficient(self):
        k, h = gm.generator(2, 1), lorentz(2)
        zero = gm.matrix(((0, 0), (0, 0)))
        self.assertEqual(gm.add(gm.mul(gm.transpose(k), h), gm.mul(h, k)), zero)
        self.assertEqual(gm.torsion_2d(gm.identity(2), zero, zero, zero, k), (0, -1))

    def test_closed_emk_map_matches_generic_cayley(self):
        for rotation in (False, True):
            g = gm.generator(3, 2, rotation=rotation)
            self.assertEqual(gm.emk_step(3, 2, Q(-2, 3), rotation=rotation),
                             gm.cayley(g, Q(-2, 3)))

    def test_invalid_carriers_and_singular_maps_are_rejected(self):
        with self.assertRaises(ValueError):
            gm.emk_step(2, 1, 1)
        with self.assertRaises(ValueError):
            gm.emk_step(2, 0, Q(1, 2))
        with self.assertRaises(ValueError):
            gm.invariant_forms(2, (((1, 0), (0, 0)),))
        with self.assertRaises(ValueError):
            gm.descend_form(lorentz(2), ((1, 0), (1, 0)))
        with self.assertRaises(ValueError):
            warped_chart(0, (0, 1, 0))


if __name__ == '__main__':
    unittest.main()
