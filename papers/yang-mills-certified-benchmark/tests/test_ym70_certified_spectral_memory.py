import copy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ym70_certified_spectral_memory as z


class CertifiedSpectralMemoryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = z.run()
        cls.atoms = z.leading_atoms(3)
        cls.moments = z.atom_moments(cls.atoms, 11)

    def test_fresh_result_and_pin(self):
        self.assertEqual(z.check(self.cert), z.PIN.read_text().strip())

    def test_python_312_required(self):
        with patch.object(z.native.previous.sys, 'version_info', (3, 11, 0)):
            with self.assertRaises(RuntimeError):
                z.run()

    def test_upstream_sources_preserved(self):
        manifest = json.loads(z.SOURCES.read_text())
        self.assertEqual(len(manifest['upstream_sha256']), 225)
        self.assertEqual(len(manifest['local_inputs']), 4)
        self.assertEqual(manifest['base_tree'], '3e1817f48d1c5018d558ebe3c25a98743db73fb1')
        for stage in ('YM64_VACUUM_ADAPTED_TOWER.md', 'YM69_NEARBY_CONDITIONED_LOCALITY.md'):
            self.assertIn('papers/yang-mills-certified-benchmark/'+stage, manifest['upstream_sha256'])

    def test_changed_upstream_is_rejected(self):
        with patch.object(z.native, 'digest', return_value='0'*64):
            with self.assertRaises(ValueError):
                z.source_checks()

    def test_rehashed_continuum_overclaim_is_rejected(self):
        wrong = copy.deepcopy(self.cert)
        wrong['model']['continuum_R4_QFT'] = True
        with tempfile.TemporaryDirectory() as directory:
            result, pin = Path(directory)/'result', Path(directory)/'pin'
            result.write_text(json.dumps(wrong))
            pin.write_text(z.p.canonical_sha(wrong))
            with self.assertRaises(ValueError):
                z.check(self.cert, result, pin)

    def test_evidence_contains_no_binary_floats(self):
        def walk(value):
            self.assertNotIsInstance(value, float)
            if isinstance(value, dict):
                for child in value.values(): walk(child)
            elif isinstance(value, (list, tuple)):
                for child in value: walk(child)
        walk(self.cert)

    def test_native_single_and_double_shared_paths(self):
        rows = z.native_controls()['native_shared_path_algebra']
        self.assertEqual([r['response_norm_squared'] for r in rows], [F(3, 256), F(3, 64)])
        self.assertTrue(all(r['native_derivative_and_individual_link_Casimir_checks'] for r in rows))

    def test_actual_lattice_incidence_and_static_coefficient(self):
        row = z.native_controls()['actual_cubic_lattice_counts'][0]
        self.assertEqual((row['boundary_links'], row['one_link_plaquettes'], row['two_link_plaquettes']), (12, 40, 4))
        self.assertEqual(row['response_norm_squared'], F(21, 32))
        self.assertEqual(row['static_coefficient'], F(113, 4095))
        self.assertEqual(self.atoms, ((F(9), F(1, 16)), (F(21, 2), F(5, 32)),
                                     (F(25, 2), F(5, 96)), (F(13), F(1, 48))))
        self.assertEqual(sum(weight for rate, weight in self.atoms), F(7, 24))

    def test_declared_block_domain(self):
        for factor in (0, 1, 2, -3, True, 3.0):
            with self.assertRaises(ValueError): z.leading_atoms(factor)

    def test_single_rate_is_exact_at_order_zero(self):
        atoms = ((F(7, 3), F(5, 11)),)
        M = z.atom_moments(atoms, 3)
        for s in (0, F(1, 17), 1000):
            row = z.enclose(M, 0, s, F(7, 3))
            self.assertEqual(row['lower'], z.atom_response(atoms, s))
            self.assertEqual(row['upper'], row['lower'])

    def test_scalar_order_zero_matches_closed_form(self):
        m0, m1, m2 = [value[0][0] for value in self.moments[:3]]
        s = F(17)
        row = z.enclose(self.moments, 0, s, 9)
        self.assertEqual(row['lower'][0][0], m0*m0/(m0+s*m1))
        self.assertEqual(row['upper'][0][0], m0-s*m1*m1/(m1+s*m2))
        self.assertLess(row['lower'][0][0], z.atom_response(self.atoms, s)[0][0])
        self.assertGreater(row['upper'][0][0], z.atom_response(self.atoms, s)[0][0])

    def test_four_atoms_require_four_cyclic_directions(self):
        for order in range(5):
            row = z.enclose(self.moments, order, 1000, 9)
            self.assertEqual(row['rank'], min(order+1, 4))
            self.assertEqual(row['weighted_rank'], row['rank'])
            if order >= 3:
                self.assertEqual(row['lower'], z.atom_response(self.atoms, 1000))
                self.assertEqual(row['upper'], row['lower'])
            else:
                self.assertLess(row['lower'][0][0], row['upper'][0][0])

    def test_zero_parameter_is_exact_at_every_order(self):
        for order in range(5):
            row = z.enclose(self.moments, order, 0, 9)
            self.assertEqual(row['lower'], self.moments[0])
            self.assertEqual(row['upper'], self.moments[0])
            self.assertEqual(row['error_factor'], 0)

    def test_zero_and_redundant_readouts(self):
        zero = [z.a.scale(z.a.identity(3), 0) for _ in range(5)]
        row = z.enclose(zero, 1, 10, 1)
        self.assertEqual(row['rank'], 0)
        self.assertEqual(row['lower'], zero[0])
        self.assertEqual(row['upper'], zero[0])
        gram = z.previous.matrix(((1, 2, 0), (2, 4, 0), (0, 0, 0)))
        M = [z.a.scale(gram, 1/F(3)**n) for n in range(7)]
        row = z.enclose(M, 2, 17, 3)
        self.assertEqual(row['rank'], 1)
        self.assertEqual(row['lower'], z.a.scale(gram, F(3, 20)))
        self.assertEqual(row['upper'], row['lower'])

    def test_noncommuting_ports_and_change_of_readout_basis(self):
        H, C, B, D = z.previous.fixture()
        M = [z.previous.moment(B, D, n) for n in range(5)]
        transform = z.previous.matrix(((1, 0, 0), (0, 1, 2), (0, -1, 3)))
        change = lambda v: z.a.multiply(z.a.multiply(z.a.transpose(transform), v), transform)
        changed = [change(v) for v in M]
        for order in (0, 1):
            first = z.enclose(M, order, 17, 2)
            second = z.enclose(changed, order, 17, 2)
            self.assertEqual(second['lower'], change(first['lower']))
            self.assertEqual(second['upper'], change(first['upper']))

    def test_ordinary_singular_inverse_is_not_used(self):
        G = z.previous.matrix(((1, 1), (1, 1)))
        with self.assertRaises(ValueError): z.a.inverse(G)
        value, rank = z.quotient_form(G, z.a.scale(G, 3), [[F(1)], [F(1)]])
        self.assertEqual((value, rank), ([[F(1, 3)]], 1))

    def test_indefinite_gram_with_zero_diagonal_is_rejected(self):
        for G in (((0, 1), (1, 0)), ((1, 2), (2, 1)), ((1, 2), (0, 1))):
            with self.assertRaises(ValueError): z.psd_basis(G)

    def test_kernel_inconsistent_rhs_and_pencils_are_rejected(self):
        G = z.previous.matrix(((1, 1), (1, 1)))
        for pencil, rhs in ((G, [[1], [0]]), (z.a.identity(2), [[1], [1]]),
                            (z.a.scale(G, 0), [[1], [1]])):
            with self.assertRaises(ValueError): z.quotient_form(G, pencil, rhs)

    def test_moment_dimensions_counts_and_support_are_checked(self):
        with self.assertRaises(ValueError): z.enclose(self.moments[:2], 0, 1, 9)
        with self.assertRaises(ValueError): z.enclose([[[1]], [[1, 0]], [[1]]], 0, 1, 1)
        with self.assertRaises(ValueError): z.enclose(self.moments, 0, 1, 100)
        with self.assertRaises(ValueError): z.enclose([[[-1]], [[1]], [[1]]], 0, 1, 1)

    def test_inputs_remain_exact_and_in_declared_domain(self):
        for s in (-1, True, 0.1):
            with self.assertRaises(ValueError): z.enclose(self.moments, 0, s, 9)
        for gap in (0, -1, True, 9.0):
            with self.assertRaises(ValueError): z.error_factor(0, 1, gap)
        for order in (-1, True, 2.0):
            with self.assertRaises(ValueError): z.enclose(self.moments, order, 1, 9)
        with self.assertRaises(ValueError): z.enclose([[[0.1]]]*3, 0, 1, 1)

    def test_chebyshev_error_matches_explicit_first_orders(self):
        self.assertEqual(z.error_factor(0, 10, 1), F(25, 36))
        self.assertEqual(z.error_factor(1, 10, 1), F(625, 2209))
        self.assertEqual(z.error_factor(3, 0, 1), 0)

    def test_error_improves_with_order_and_gap_and_worsens_with_band(self):
        for order in (0, 1, 3, 15):
            self.assertLess(z.error_factor(order+1, 100, 1), z.error_factor(order, 100, 1))
            self.assertLess(z.error_factor(order, 100, 2), z.error_factor(order, 100, 1))
            self.assertLess(z.error_factor(order, 10, 1), z.error_factor(order, 100, 1))
        self.assertLess(z.error_factor(63, 1000, 1), F(1, 500))

    def test_error_width_uses_one_M0_budget(self):
        for s in (F(1, 7), F(9), F(1000)):
            for order in (0, 1, 2):
                row = z.enclose(self.moments, order, s, 9)
                z.assert_enclosure(row, z.atom_response(self.atoms, s), self.moments[0])

    def test_positive_parameter_bounds_survive_Taylor_failure(self):
        row = z.leading_controls()['outside_Taylor_disk']
        self.assertTrue(row['errors_strictly_increase'])
        exact = z.atom_response(self.atoms, 1000)
        self.assertEqual(z.enclose(self.moments, 3, 1000, 9)['lower'], exact)

    def test_clock_covariance_including_error_factor(self):
        kappa, s = F(2, 7), F(23)
        scaled = [z.a.scale(v, kappa**(1-n)) for n, v in enumerate(self.moments)]
        left = z.enclose(scaled, 1, s, 9*kappa)
        right = z.enclose(self.moments, 1, s/kappa, 9)
        self.assertEqual(left['lower'], z.a.scale(right['lower'], kappa))
        self.assertEqual(left['upper'], z.a.scale(right['upper'], kappa))
        self.assertEqual(left['error_factor'], right['error_factor'])

    def test_nested_order_and_omitted_channel_counterexample(self):
        H, C, B, D = z.previous.fixture()
        M = [z.previous.moment(B, D, n) for n in range(3)]
        s = F(1)
        row = z.enclose(M, 0, s, 2)
        lower = z.previous.effective(z.a.add(C, row['upper'], -1), (0, 1), s)
        upper = z.previous.effective(z.a.add(C, row['lower'], -1), (0, 1), s)
        direct = z.previous.effective(H, (0, 1), s)
        self.assertTrue(z.a.psd(z.a.add(direct, lower, -1)))
        self.assertTrue(z.a.psd(z.a.add(upper, direct, -1)))
        sampled = z.previous.sub(z.previous.effective(H, (0, 1, 2), s), (0, 1))
        self.assertNotEqual(sampled, direct)
        self.assertGreater(sampled[1][1], direct[1][1])

    def test_claims_distinguish_leading_coefficient_from_full_memory(self):
        claim = self.cert['model']
        self.assertTrue(claim['actual_leading_Wilson_memory_computed'])
        self.assertTrue(claim['full_interacting_positive_parameter_enclosure_theorem'])
        for field in ('full_interacting_moments_numerically_computed',
                      'four_rates_for_full_interacting_memory', 'uniform_perturbative_remainder',
                      'enlarged_bare_coupling_window', 'full_spatial_Schur_memory_locality',
                      'continuum_R4_QFT', 'physical_mass_prediction'):
            self.assertFalse(claim[field])


if __name__ == '__main__':
    unittest.main()
