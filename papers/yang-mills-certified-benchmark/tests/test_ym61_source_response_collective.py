import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym61_source_response_collective as z


class SourceResponseCollectiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()

    def test_fresh_exact_certificate_is_pinned(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_all_cycle_inverses_have_centered_sources_and_zero_residual(self):
        cycles = self.cert['source']['cycles']
        self.assertEqual([F(c['leading']) for c in cycles], [F(1, 2), F(11, 24), F(49, 108), F(4853, 10368)])
        witnesses = [w for c in cycles for w in c['witnesses']]
        self.assertEqual(len(witnesses), 52)
        for w in witnesses:
            self.assertTrue(w['centered'])
            self.assertEqual(F(w['inverse_residual_norm2']), 0)
            self.assertTrue(w['active'])
            self.assertFalse(set(w['active']) & set(w['completed']))

    def test_no_proper_cycle_subset_can_have_a_vacuum_energy_coefficient(self):
        for row in self.cert['parity']['cycle_parity']:
            self.assertEqual(row['even_subsets'], 2)
        self.assertEqual(z.cycle_edges(1), ((0, 1), (0, 1)))
        self.assertEqual(z.cycle_response(1)['leading'], F(1, 2))

    def test_path_parity_forces_minimum_order_and_skips_opposite_parity(self):
        self.assertEqual(self.cert['parity']['path_incidence_cases'], 120)
        for r in range(1, 6):
            edges = z.cycle_edges(r)
            full = (1 << len(edges))-1
            for e in range(len(edges)):
                degrees = z.incidences(edges, full ^ (1 << e), r+1)
                self.assertEqual(sum(d % 2 for d in degrees), 2)
            self.assertEqual(z.previous.secant_error_order(2*r+4, r+2), r+2)

    def test_triangular_response_agrees_with_six_order_closed_formula(self):
        for k in (F(1, 3), F(1), F(2)):
            self.assertEqual(z.cycle_response(2, k)['leading'], z.triangular_response(k))
        self.assertEqual(z.triangular_response(1), F(11, 24))

    def test_higher_correlation_is_neither_factorized_nor_given_by_one_rate(self):
        self.assertNotEqual(z.cycle_response(2)['leading'], z.cycle_response(1)['leading']**2)
        self.assertNotEqual(z.cycle_response(2)['leading'], F(3, 8))
        self.assertEqual(F(self.cert['source']['wrong_factorized_correlation_gives']), F(1, 4))

    def test_safe_site_elimination_is_checked_against_full_source_inverses(self):
        for row in self.cert['source']['inverse_projection_checks']:
            self.assertGreater(row['unprojected_dimension'], 1)
            self.assertEqual(F(row['projection_inverse_residual_norm2']), 0)
        q = z.source.overlap(0, 1, 3)
        r = z.source.overlap(1, 2, 3)
        self.assertEqual(z.source.site_average(z.p.mul(q, r), 1, 3), z.p.scale(z.source.overlap(0, 2, 3), F(1, 4)))

    def test_averaging_an_unfinished_site_destroys_a_real_source(self):
        q = z.source.overlap(0, 1, 2)
        self.assertTrue(q)
        self.assertFalse(z.source.site_average(q, 0, 2))
        self.assertTrue(self.cert['source']['premature_averaging_erases_source'])

    def test_separated_endpoints_have_a_different_native_rate(self):
        self.assertEqual(z.source.omega(1), 1)
        self.assertEqual(z.separated_rate(1), F(3, 4))
        self.assertEqual(len(self.cert['source']['source_susceptibilities']), 12)
        self.assertEqual(1/(2*z.separated_rate(1)), F(2, 3))

    def test_distinct_source_pairs_have_zero_quadratic_cross_response(self):
        x, y = z.source.overlap(0, 1, 3), z.source.overlap(0, 2, 3)
        inverse_y = z.p.scale(y, 1/z.separated_rate(1))
        self.assertEqual(z.pairing(x, inverse_y, 3), 0)
        self.assertGreater(z.pairing(y, inverse_y, 3), 0)

    def test_first_collective_response_cancels_rate_dependence(self):
        rows = self.cert['collective']['first_response_cells']
        self.assertEqual(len({r['free_rate'] for r in rows}), 4)
        self.assertEqual({r['first_matrix_coefficient'] for r in rows}, {'-1/4'})
        k = F(17, 11)
        self.assertEqual((1+3*k)/8/(2*z.source.omega(k)), F(1, 4))

    def test_inverse_square_root_congruence_is_exact_through_cubic_order(self):
        for L in (2, 3, 4, 5, 6):
            data = z.matrix_tower(L)
            I = z.a.identity(L)
            self.assertEqual(data['normalization'], [I]+[z.a.scale(I, 0)]*3)
            self.assertEqual(data['response'][1], z.a.scale(z.band(L, 1), F(-1, 4)))

    def test_replacing_the_anticommutator_by_twice_one_order_breaks_symmetry(self):
        J1, J2 = z.band(4, 1), z.band(4, 2)
        AB, BA = z.a.multiply(J1, J2), z.a.multiply(J2, J1)
        commutator = z.a.add(AB, z.a.scale(BA, -1))
        wrong = z.a.add(z.matrix_tower(4)['response'][3], z.a.scale(commutator, F(77, 768)))
        self.assertNotEqual(wrong, z.a.transpose(wrong))
        self.assertEqual(z.matrix_tower(4)['response'][3], z.a.transpose(z.matrix_tower(4)['response'][3]))

    def test_two_site_series_reproduces_the_previous_result(self):
        self.assertEqual(z.rational_probe_series((1, 1)), tuple(map(F, z.previous.response_controls()['trial_rate_series'])))

    def test_four_site_readout_uses_higher_correlations_and_improves_first_order(self):
        result = z.rational_probe_series((1, 2, 2, 1))
        self.assertEqual(result, (F(3, 8), F(-2, 5), F(73, 400), F(-1293529, 11844000)))
        self.assertEqual(result[1]-z.rational_probe_series((1, 1))[1], F(-3, 20))

    def test_three_site_sine_mode_and_second_order_are_exact_without_rounded_roots(self):
        J = z.band(3, 1)
        v, w = (F(1), F(0), F(1)), (F(0), F(1), F(0))
        mv = lambda M, b: tuple(sum((x*y for x, y in zip(row, b)), F(0)) for row in M)
        self.assertEqual(mv(J, v), tuple(2*x for x in w))
        self.assertEqual(mv(J, w), v)
        B2 = z.matrix_tower(3)['response'][2]
        self.assertEqual(sum(x*y for x, y in zip(v, mv(B2, w))), 0)
        second = (z.quadratic(v, B2)+2*z.quadratic(w, B2))/4
        self.assertEqual(second, F(21, 128))

    def test_zero_padding_embeds_smaller_readout_exactly(self):
        c = (F(1), F(1, 3), F(1, 9), F(1, 27))
        self.assertEqual(z.probe_rate((1, 2, 1), c[:3]), z.probe_rate((1, 2, 1, 0), c))

    def test_coupling_sign_is_implemented_by_staggering(self):
        c = tuple(F(1, 4)**r for r in range(4))
        b = (1, 2, 2, 1)
        self.assertEqual(z.probe_rate(b, c), z.probe_rate(tuple(x*(-1)**i for i, x in enumerate(b)), tuple(x*(-1)**r for r, x in enumerate(c))))

    def test_nonphysical_covariance_and_null_readout_are_refused(self):
        with self.assertRaises(ValueError):
            z.response_matrices((1, F(9, 10), F(-9, 10)))
        with self.assertRaises(ValueError):
            z.probe_rate((1, -1), (1, 1))
        with self.assertRaises(ValueError):
            z.probe_rate((1, 2), (1,))

    def test_invalid_distances_rates_and_inexact_coefficients_are_refused(self):
        for r in (0, -1, True, 2.0):
            with self.assertRaises(ValueError):
                z.cycle_edges(r)
        for k in (0, -1, 1.0, True):
            with self.assertRaises(ValueError):
                z.separated_rate(k)
        for c in ((), (2,), (1, 0.25)):
            with self.assertRaises(ValueError):
                z.response_matrices(c)
        with self.assertRaises(ValueError):
            z.rational_probe_series((0, 0))

    def test_python_312_runtime_is_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_upstream_and_claim_tampering_are_refused(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()
        with tempfile.TemporaryDirectory() as temp:
            result, pin = Path(temp)/'result', Path(temp)/'pin'
            bad = copy.deepcopy(self.cert)
            bad['model']['uniform_interacting_lower_gap'] = True
            result.write_text(json.dumps(bad))
            pin.write_text(z.native.previous.native.chain.old.canonical_sha(bad))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_all_160_upstream_files_and_four_local_inputs_are_declared(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 160)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertIn(str(z.previous.RESULT.relative_to(z.ROOT)), manifest['upstream_sha256'])

    def test_fixed_order_hierarchy_does_not_assert_uniform_convergence_or_mass(self):
        model = self.cert['model']
        for key in ('actual_scalar_operator_retained', 'independent_sources_removed_after_response', 'fixed_distance_uniform_asymptotics', 'collective_first_response_independent_of_k', 'cubic_matrix_response_noncommutative'):
            self.assertTrue(model[key])
        for key in ('uniform_in_distance_series', 'exact_mass_eigenvalue', 'uniform_interacting_lower_gap', 'physical_4d_yang_mills'):
            self.assertFalse(model[key])
        self.assertIn('LOWER_GAP_OPEN', self.cert['claim_status'])


if __name__ == '__main__':
    unittest.main()
