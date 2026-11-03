"""Independent finite obstructions and exact shell-count identities."""
from collections import Counter
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import shell_growth_model as sh


class ShellGrowthTests(unittest.TestCase):
    def test_same_carrier_valence_and_stiffness_different_capacity(self):
        a, b = sh.torus2(7), sh.circulant(49, (1, 2))
        for n, edges in [a, b]:
            self.assertEqual(n, 49)
            self.assertEqual(len(edges), 98)
            self.assertEqual({c for _, _, c in edges}, {1})
            degree = Counter(x for u, v, _ in edges for x in (u, v))
            self.assertEqual(set(degree.values()), {4})
        self.assertEqual(sh.shell_profile(*a, radius=3)['capacities'], (4, 12, 20))
        self.assertEqual(sh.shell_profile(*b, radius=3)['capacities'], (4, 6, 6))

    def test_tree_and_grid_separate_growth_even_with_same_inner_degree(self):
        n, links = sh.tree_ball(4, 3)
        self.assertEqual(sh.shell_profile(n, links, radius=3)['capacities'], (4, 12, 36))
        depth = sh.distances(n, links)
        degree = Counter(x for u, v, _ in links for x in (u, v))
        self.assertTrue(all(degree[i] == 4 for i, k in enumerate(depth) if k < 3))

    def test_cubic_count_has_quadratic_leading_term_and_finite_correction(self):
        n, links = sh.lattice_ball(3, 3)
        self.assertEqual(sh.shell_profile(n, links, radius=3)['capacities'], (6, 30, 78))
        for radius, value in enumerate((6, 30, 78)):
            self.assertEqual(value, 12*radius*radius+12*radius+6)
        self.assertNotEqual(78, 12*2*2)

    def test_exact_lattice_formula_crosses_edge_enumeration(self):
        for dimension in range(1, 5):
            n, links = sh.lattice_ball(dimension, 4)
            values = sh.shell_profile(n, links, radius=4)['capacities']
            self.assertEqual(values, tuple(sh.lattice_shell_formula(dimension, r)
                                           for r in range(4)))

    def test_cubic_count_is_not_pointwise_radial_response(self):
        n, links = sh.lattice_ball(3, 3)
        d = sh.distances(n, links)
        self.assertFalse(sh.radial_equitable(n, links, radius=3))
        field = sh.stationary_shell_means(n, links, radius=3)['field']
        shell_two = {field[i] for i, k in enumerate(d) if k == 2}
        self.assertGreater(len(shell_two), 1)

    def test_even_shell_mean_drop_is_not_inverse_aggregate_count(self):
        n, links = sh.lattice_ball(3, 3)
        profile = sh.shell_profile(n, links, radius=3)
        result = sh.stationary_shell_means(n, links, radius=3)
        self.assertEqual(result['mean_drops'][0], Q(1, 6))
        self.assertEqual(result['mean_drops'][1], Q(13, 378))
        self.assertNotEqual(result['mean_drops'][1], 1/profile['capacities'][1])

    def test_transitive_tree_shells_are_radial(self):
        n, edges = sh.tree_ball(4, 3)
        self.assertTrue(sh.radial_equitable(n, edges, radius=3))
        result = sh.stationary_shell_means(n, edges, radius=3)
        self.assertEqual(result['mean_drops'], (Q(1, 4), Q(1, 12), Q(1, 36)))

    def test_torus_and_circulant_source_profiles_disagree(self):
        a, b = sh.torus2(7), sh.circulant(49, (1, 2))
        pa = sh.stationary_shell_means(*a, radius=3)
        pb = sh.stationary_shell_means(*b, radius=3)
        self.assertNotEqual(pa['means'], pb['means'])
        self.assertNotEqual(pa['mean_drops'][1], Q(1, 12))
        self.assertNotEqual(pb['mean_drops'][1], Q(1, 6))

    def test_upper_area_bound_does_not_force_saturation(self):
        n, links = sh.lattice_ball(1, 4)
        profile = sh.shell_profile(n, links, radius=4)
        self.assertEqual(profile['capacities'], (2, 2, 2, 2))
        self.assertTrue(all(c <= 2*(r+1)**2 for r, c in enumerate(profile['capacities'])))

    def test_invalid_graphs_or_radii_rejected(self):
        for constructor in [lambda: sh.torus2(6), lambda: sh.circulant(8, (1, 4)),
                            lambda: sh.lattice_ball(5, 2), lambda: sh.tree_ball(2, 2)]:
            with self.assertRaises(ValueError):
                constructor()
        n, links = sh.torus2(7)
        with self.assertRaises(ValueError):
            sh.shell_profile(n, links, radius=100)


if __name__ == '__main__':
    unittest.main()
