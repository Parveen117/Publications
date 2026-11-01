"""Adversarial finite witnesses for the separately written GE2 theorems."""
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'certificates'))
import ge2_clock_free as g


class ClockFreeRecognitionTests(unittest.TestCase):
    def test_rational_turn_independent_closed_formula_and_inverse(self):
        for r in (F(-12), F(-1, 3), F(0), F(2), F(11)):
            c, s = (16-r*r)/(16+r*r), -8*r/(16+r*r)
            self.assertEqual(g.turn((r, 0, 0)), (c, s, 0, 0))
            self.assertEqual(c*c+s*s, 1)
            self.assertEqual(g.y.qmul(g.turn((r, 0, 0)), g.turn((-r, 0, 0))), g.y.ONE)

    def test_phase_first_jet_and_even_second_jet_by_exact_interpolation(self):
        result = g.step_controls()
        self.assertEqual([r['exact_interpolation_nodes'] for r in result['independent_phase_and_second_jets']], [3, 5, 7])
        self.assertEqual(len(result['exact_operator_bound_controls']), 9)

    def test_weighted_record_lift_uses_the_declared_target_pairing(self):
        G = g.a.diag((1, 4))
        U = [[F(3, 5), F(-8, 5)], [F(2, 5), F(3, 5)]]
        model = g.lift(((F(1, 3), U), (F(2, 3), g.a.identity(2))), G)
        A, W = model['A'], model['W']
        self.assertEqual(g.a.multiply(g.a.multiply(g.a.transpose(A), W), A), G)
        self.assertEqual(g.a.multiply(g.dagger(model['J'], G, W), model['J']), g.a.identity(2))
        self.assertNotEqual(g.a.transpose(model['J']), model['Jdag'])

    def test_memory_quadratic_form_equals_direct_record_variance(self):
        G = g.y.degree_data(1)['G']
        z = g.protocol(1, ((F(1, 3), (1, 0, 0)), (F(2, 3), (0, 2, 0))))
        f = [F(1), F(-2), F(3), F(1, 2)]
        mean = g.a.apply(z['S'], f)
        variance = F(0)
        for w, U in z['actions']:
            difference = [x-y for x, y in zip(g.a.apply(U, f), mean)]
            variance += w*g.a.dot(difference, g.a.apply(G, difference))
        energy = g.a.dot(f, g.a.apply(g.a.multiply(G, g.memory(z['S'], G)), f))
        self.assertEqual(variance, energy)
        self.assertGreater(variance, 0)

    def test_full_record_chain_law_and_no_return_are_typed(self):
        result = g.record_controls()
        self.assertEqual(result['retained_words'], 4)
        self.assertEqual(result['final_event_dimension'], 16)
        self.assertTrue(result['typed_chain_law'])
        self.assertTrue(result['no_return_on_old_zero_mean_records'])

    def test_returning_memory_prevents_general_compressed_composition(self):
        result = g.distinction_controls()
        self.assertEqual(result['returning_memory'], '16/25')
        self.assertEqual(result['compressed_inverse_product'], '9/25')

    def test_ledger_is_from_raw_records_and_scales_quadratically(self):
        records = ((F(1, 3), (F(1, 2), 0, 0)), (F(2, 3), (0, F(1, 3), 0)))
        z = g.protocol(1, records)
        scaled = g.protocol(1, tuple((w, tuple(3*x for x in v)) for w, v in records))
        self.assertEqual(z['h'], F(17, 216))
        self.assertEqual(scaled['h'], 9*z['h'])
        self.assertEqual(scaled['m4'], 81*z['m4'])
        self.assertEqual(scaled['C'], z['C'])

    def test_idle_arrow_has_identity_action_without_a_shape(self):
        z = g.protocol(2, ((1, (0, 0, 0)), (0, (100, 0, 0))))
        self.assertEqual(z['h'], 0)
        self.assertEqual(z['eps2'], 0)
        self.assertEqual(z['S'], g.a.identity(10))
        self.assertIsNone(z['C'])
        with self.assertRaises(ValueError):
            g.normalized_shape(z)

    def test_same_endpoint_does_not_erase_the_history_ledger(self):
        result = g.distinction_controls()
        self.assertEqual(result['closed_endpoint_positive_ledger'], '1/4')

    def test_unequal_step_heat_and_deterministic_phase_enclosures(self):
        result = g.refinement_controls()
        self.assertEqual(len(result['irregular_heat_enclosures']), 8)
        self.assertEqual(len(result['deterministic_phase_enclosures']), 4)
        for row in result['irregular_heat_enclosures']:
            self.assertEqual(F(row['tau']), F(13, 72))
            self.assertLessEqual(F(row['error_upper']), F(row['bound_per_unit_norm']))

    def test_order_survives_equal_budgets_and_integrated_tensors(self):
        result = g.distinction_controls()
        self.assertNotEqual(F(result['ordered_generator_commutator_entry'][2]), 0)
        self.assertTrue(result['rational_word_order_detected'])

    def test_heat_product_defect_is_the_native_energy_form(self):
        C = g.a.diag((0, F(1, 2), F(1, 2)))
        f, h = g.p.COORD[0], g.p.add(g.p.COORD[1], g.p.COORD[3])
        defect = g.p.add(g.p.scale(g.q.lap(g.p.mul(f, h), C), -1),
                        g.p.mul(g.q.lap(f, C), h), g.p.mul(f, g.q.lap(h, C)))
        self.assertEqual(defect, g.p.scale(g.native.gamma(f, h, C), 2))
        self.assertTrue(defect)

    def test_memory_alone_cannot_recover_a_phase_generator(self):
        U = g.action(1, (F(1, 4), 0, 0))
        G = g.y.degree_data(1)['G']
        self.assertEqual(g.memory(U, G), g.old.zero(4))
        self.assertNotEqual(U, g.a.identity(4))
        self.assertEqual(g.memory(g.a.scale(U, -1), G), g.memory(U, G))

    def test_large_symmetric_step_can_be_negative_but_stays_contractive(self):
        z = g.protocol(1, ((1, (8, 0, 0)),))
        self.assertEqual(z['S'], g.a.scale(g.a.identity(4), F(-3, 5)))
        self.assertEqual(g.memory(z['S'], g.y.degree_data(1)['G']),
                         g.a.scale(g.a.identity(4), F(16, 25)))

    def test_response_curvature_and_clock_normalization(self):
        rows = g.curvature_controls()['native_response_fixtures']
        self.assertEqual([r['record_clock_rate'] for r in rows], ['1/4', '1/5', '1/10', '0'])
        self.assertEqual(rows[0]['kappa'], '1')
        self.assertEqual(rows[-1]['kappa'], '0')

    def test_type_and_scope_gates(self):
        self.assertEqual(len(g.refusal_controls()), 12)

    def test_source_pin_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            file = root/'source.py'
            file.write_text('pinned\n')
            pins = {'upstream_sha256': {'source.py': g.old.digest(file)}, 'local_inputs': []}
            g.source_checks(root, pins)
            file.write_text('modified\n')
            with self.assertRaises(ValueError):
                g.source_checks(root, pins)

    def test_certificate_cannot_be_promoted_by_editing_claims(self):
        cert = {'verdict': 'PASS', 'scope': 'declared free record protocol'}
        with tempfile.TemporaryDirectory() as folder:
            result, expected = Path(folder)/'result.json', Path(folder)/'expected.sha256'
            result.write_text(json.dumps(cert))
            expected.write_text(g.y.canonical_sha(cert))
            g.check(cert, result, expected)
            with self.assertRaises(ValueError):
                g.check(dict(cert, scope='physical clock uniquely derived'), result, expected)


if __name__ == '__main__':
    unittest.main()
