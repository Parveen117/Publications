import copy
from decimal import Decimal, localcontext
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
sys.path.insert(0, str(ROOT/'papers/yang-mills-certified-benchmark/certificates'))
import nsb1_nonlinear_blindness as n
import ym62_gauge_lattice_gap as native


def decimal(value):
    value = F(value)
    return Decimal(value.numerator)/Decimal(value.denominator)


def recorded_polynomial(records):
    return {tuple(row['exponent']): F(row['coefficient']) for row in records}


class GaugeAndSourceCorrectionTests(unittest.TestCase):
    def test_native_quaternion_derivatives_independently_verify_flattened_parent(self):
        p = native.p
        one = {(0, 0, 0, 0): F(1)}
        w = p.COORD[0]
        s = p.add(one, p.scale(p.mul(w, w), -1))
        points = ((F(1), 0, 0, 0), (F(-1), 0, 0, 0),
                  (F(3, 5), F(4, 5), 0, 0),
                  (F(1, 3), F(2, 3), F(2, 3), 0))
        for N in (0, 1, 2, 4):
            power = one
            for _ in range(N+1):
                power = p.mul(power, s)
            f = p.add(one, p.scale(power, -1))
            electric = p.scale(native.electric(f, 1), 4)
            gamma = p.scale(native.energy(f, 1), 4)
            row = n.gauge_parent(N, F(3, 2))
            le = recorded_polynomial(row['electric_f'])
            ga = recorded_polynomial(row['Gamma_f'])
            potential = recorded_polynomial(row['multiplication_potential'])
            for point in points:
                self.assertEqual(sum(F(x)*F(x) for x in point), 1)
                expected_le = native.evaluate(electric, (point,))
                expected_ga = native.evaluate(gamma, (point,))
                self.assertEqual(n.jets.poly_evaluate(le, (F(point[0]),)), expected_le)
                self.assertEqual(n.jets.poly_evaluate(ga, (F(point[0]),)), expected_ga)
                self.assertEqual(n.jets.poly_evaluate(potential, (F(point[0]),)),
                                 -F(3, 4)*expected_le+F(9, 16)*expected_ga)

    def test_arbitrarily_deep_local_jet_matching_does_not_identify_a_gap(self):
        for N in (1, 2, 3, 7):
            row = n.gauge_parent(N, F(5, 3))
            radial = recorded_polynomial(row['potential_in_s'])
            self.assertEqual(min(key[0] for key in radial), N)
            self.assertEqual(radial[(N,)], -F(5, 3)*(N+1)*(2*N+3))
            self.assertEqual(row['matched_free_potential_jet_order'], 2*N-1)
            self.assertFalse(row['original_homogeneous_Yang_Mills_model'])
            self.assertTrue(row['exact_vacuum_requires_one_or_edge_disjoint_active_plaquettes'])
        self.assertEqual(n.gauge_parent(7, 0)['potential_in_s'], [])
        self.assertIsNone(n.gauge_parent(0, 2)['matched_free_potential_jet_order'])

    def test_rayleigh_exponent_and_cap_retain_both_wells(self):
        for N in range(6):
            row = n.gauge_parent(N, 4)
            expected = (3**(N+1)-1)/F(4**(N+1))
            self.assertEqual(row['decay_exponent'], expected)
            self.assertEqual(row['denominator_density_exponent']-row['numerator_density_exponent'], expected)
            self.assertEqual(4/row['Haar_cap_lower'], 162)
            self.assertTrue(row['extra_W_even_restriction_would_exclude_the_trial'])
        self.assertEqual(n.gauge_parent(0, 0)['decay_exponent'], F(1, 2))
        self.assertEqual(n.gauge_parent(2, 0)['decay_exponent'], F(13, 32))

    def test_directed_exponential_encloses_independent_decimal_values(self):
        with localcontext() as context:
            context.prec = 75
            for value in (F(0), F(1, 10), F(13, 2), F(26), F(80)):
                interval = n.exp_negative_interval(value)
                actual = (-decimal(value)).exp()
                self.assertLessEqual(decimal(interval['lower']), actual)
                self.assertGreaterEqual(decimal(interval['upper']), actual)
                self.assertGreater(interval['lower'], 0)
        upper16 = n.gauge_parent(2, 16)['physical_gauge_invariant_gap_upper']
        upper64 = n.gauge_parent(2, 64)['physical_gauge_invariant_gap_upper']
        self.assertLess(upper16, 1)
        self.assertLess(upper64, upper16)
        self.assertLess(upper64, F(1, 10**8))

    def test_negative_compression_supremum_cannot_cancel_the_complement_tail(self):
        diagonal = [F(-100), F(2), F(0)]
        zero_extended = [diagonal[0], F(0), F(0)]
        mu = diagonal[0]
        tail = max(abs(a-b) for a, b in zip(diagonal, zero_extended))
        self.assertEqual(tail, 2)
        self.assertLess(mu+tail, 1)
        self.assertLess(min(1-value for value in diagonal), 0)
        self.assertEqual(n.corrected_compression_bound(mu, tail), max(diagonal))
        self.assertGreater(n.corrected_compression_bound(mu, tail), 1)
        self.assertEqual(n.corrected_compression_bound(F(3, 4), F(1, 8)), F(7, 8))
        self.assertTrue(n.source_erratum_controls()['printed_condition_would_falsely_accept'])

    def test_invalid_exact_parameters_are_rejected(self):
        for value in (True, 0.5, -1):
            with self.assertRaises(ValueError):
                n.exp_negative_interval(value)
        with self.assertRaises(ValueError):
            n.exp_negative_interval(10, terms=2)
        for N, t in ((-1, 0), (True, 1), (2, -1), (2, 0.5)):
            with self.assertRaises(ValueError):
                n.gauge_parent(N, t)
        with self.assertRaises(ValueError):
            n.corrected_compression_bound(0, -1)


class EvidenceTests(unittest.TestCase):
    def test_fresh_certificate_matches_result_and_pin(self):
        cert = n.run()
        self.assertEqual(n.check(cert), n.PIN.read_text().strip())

    def test_source_manifest_binds_the_two_pdfs_and_all_new_proof_inputs(self):
        manifest = json.loads(n.SOURCES.read_text())
        self.assertIn('SPECTRAL 1.pdf', manifest['upstream_sha256'])
        self.assertIn('SPECTRAL_2.pdf', manifest['upstream_sha256'])
        self.assertIn('papers/nonlinear-spectral-blindness/SOURCE_AUDIT.md', manifest['local_inputs'])
        self.assertIn('papers/nonlinear-spectral-blindness/YANG_MILLS_DOCK.md', manifest['local_inputs'])
        self.assertEqual(manifest['base_tree'], n.BASE_TREE)
        checked = n.source_checks()
        self.assertEqual(len(checked), len(manifest['upstream_sha256'])+len(manifest['local_inputs'])+1)

    def test_changed_upstream_or_source_base_is_rejected(self):
        manifest = json.loads(n.SOURCES.read_text())
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'manifest.json'
            bad = copy.deepcopy(manifest)
            bad['upstream_sha256']['SPECTRAL_2.pdf'] = '0'*64
            path.write_text(json.dumps(bad))
            with patch.object(n, 'SOURCES', path), self.assertRaisesRegex(ValueError, 'Changed upstream'):
                n.source_checks()
            bad = copy.deepcopy(manifest); bad['base_tree'] = '0'*40
            path.write_text(json.dumps(bad))
            with patch.object(n, 'SOURCES', path), self.assertRaisesRegex(ValueError, 'base'):
                n.source_checks()

    def test_rehashed_overclaims_are_rejected(self):
        cert = n.run()
        for key in ('Jordan_structure_claimed_blind',
                    'entire_Hessian_eigenvalue_field_determines_flatness',
                    'all_physical_symmetries_imply_all_curvatures_zero',
                    'exponential_pair_has_gauge_inequivalent_selected_connections',
                    'native_counterfamily_shares_the_full_quantum_spectrum',
                    'positive_four_dimensional_continuum_Yang_Mills_gap_proved'):
            bad = copy.deepcopy(cert); bad['model'][key] = True
            with tempfile.TemporaryDirectory() as folder:
                result, pin = Path(folder)/'result.json', Path(folder)/'pin'
                result.write_text(json.dumps(bad)); pin.write_text(n.canonical_sha(bad))
                with self.assertRaisesRegex(ValueError, 'scope'):
                    n.check(bad, result, pin)

    def test_exact_json_and_component_claim_boundaries(self):
        cert = n.run()
        self.assertEqual(json.loads(json.dumps(cert)), cert)
        self.assertFalse(cert['controls']['jet_scope']['infinite_Koopman_spectrum_identified'])
        self.assertFalse(cert['controls']['geometry_scope']['entire_Hessian_spectral_field_determines_flatness'])
        def walk(value):
            self.assertNotIsInstance(value, (float, F))
            if isinstance(value, dict):
                for key, item in value.items():
                    self.assertIsInstance(key, str)
                    walk(item)
            elif isinstance(value, list):
                for item in value:
                    walk(item)
        walk(cert)

    def test_python312_and_active_assertions_required(self):
        with patch.object(n.sys, 'version_info', (3, 11)), self.assertRaises(RuntimeError):
            n.runtime_check()
        with patch.object(n.sys, 'flags', SimpleNamespace(optimize=1)), self.assertRaises(RuntimeError):
            n.runtime_check()


if __name__ == '__main__':
    unittest.main()
