from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym51_heat_selection as y


class HeatSelectionTests(unittest.TestCase):
    def test_positive_protocol_energy_on_all_checked_sources(self):
        r=y.protocol_controls()
        self.assertEqual(r['factorization_checks'],25)
        self.assertEqual(r['all_source_energy_checks'],25)
        self.assertTrue(r['YM50_bound_recovered'])

    def test_invalid_weights_vectors_and_inexact_inputs_are_refused(self):
        for records in ([],[(F(1,2),(1,0,0))],
                        [(F(-1),(1,0,0)),(2,(0,1,0))],
                        [(1,(1,0))],[(1.0,(1,0,0))],[(1,(0.1,0,0))]):
            with self.assertRaises(ValueError):y.protocol(records)
        for eta in (F(-1,5),F(6,5),0.25,True):
            with self.assertRaises(ValueError):y.eta_protocol(eta)

    def test_invalid_tensor_positivity_and_symmetry_are_refused(self):
        for C in ([[1,0],[0,1]],[[1,1,0],[0,1,0],[0,0,1]],
                  [[1,2,0],[2,1,0],[0,0,1]]):
            with self.assertRaises(ValueError):y.tensor(C)

    def test_unchanged_reference_and_native_brackets(self):
        r=y.invariance_controls()
        self.assertEqual(r['same_reference_turn_checks'],70)
        self.assertEqual(r['unchanged_native_bracket_checks'],12)

    def test_covariance_is_distinct_from_fixed_law_isotropy(self):
        r=y.covariance_controls()
        self.assertEqual(r['covariant_operator_checks'],60)
        self.assertEqual(r['fixed_law_constraint_inertia'],[5,0,1])
        self.assertTrue(r['covariance_does_not_imply_isotropy'])
        with self.assertRaises(ValueError):
            y.rotation((F(2),F(0),F(0),F(0)))

    def test_full_support_does_not_give_a_uniform_rate(self):
        r=y.selection_controls()
        self.assertEqual(r['linear_blindness_checks'],24)
        self.assertEqual(r['quadratic_witness_checks'],6)
        for row in r['family'][1:]:
            self.assertTrue(row['full_axis_support'])
            self.assertEqual(row['degree_two_energy_inertia'],[9,0,1])
        self.assertEqual(r['family'][0]['degree_two_energy_inertia'],[6,0,4])

    def test_same_linear_generator_and_different_quadratic_action(self):
        first=y.moments(y.eta_protocol(F(1)))[0]
        second=y.moments(y.eta_protocol(F(1,1000)))[0]
        self.assertEqual(y.operator(1,first),y.operator(1,second))
        self.assertNotEqual(y.operator(2,first),y.operator(2,second))
        self.assertEqual(y.lap(y.SLOW,second),y.p.scale(y.SLOW,F(1,500)))
        self.assertEqual(y.y.phi(y.p.mul(y.SLOW,y.SLOW)),F(1,12))

    def test_six_quadratic_probes_recover_off_diagonal_tensor_entries(self):
        r=y.recovery_controls()
        self.assertEqual(r['quadratic_tensor_recovery_checks'],18)
        self.assertTrue(r['diagonal_only_bank_has_positive_counterexample'])
        self.assertEqual(r['minimum_fixed_linear_channels'],6)

    def test_weighted_protocol_matches_raw_repeated_records(self):
        r=y.count_controls()
        self.assertEqual(r['independent_labelled_count_checks'],14)
        self.assertEqual(r['two_step_labels'],5184)

    def test_outward_refinement_bounds_for_independent_scalar_modes(self):
        r=y.refinement_controls()
        self.assertEqual(r['outward_scalar_refinement_checks'],72)
        for row in r['cases']:
            self.assertLessEqual(F(row['error_upper']),F(row['bound']))

    def test_finite_microsteps_differ_despite_equal_second_moment(self):
        r=y.clock_and_micro_controls()
        self.assertTrue(r['same_second_moment_different_microsteps'])
        self.assertEqual(r['microstep_fourth_coefficient_difference'],'1/128')
        first=y.y.cosine_square(F(1,16))
        second=y.y.Iv(F(3,4))+y.y.Iv(F(1,4))*y.y.cosine_square(F(1,4))
        self.assertGreater(second.lo,first.hi)

    def test_clock_rescaling_changes_relative_interaction(self):
        r=y.clock_and_micro_controls()
        self.assertEqual(r['clock_rescaling_checks'],12)
        for row in r['interaction_rescalings']:
            c,theta,step,b,relative=map(F,(row['c'],row['theta'],row['a'],
                                          row['b'],row['relative_interaction']))
            self.assertEqual(c*step,b)
            self.assertEqual(b*relative,step*theta)
            self.assertNotEqual(theta,relative)

    def test_bound_domain_and_YM50_normalization(self):
        self.assertEqual(y.budget(2,F(1),32,F(3),F(9)),F(3,16))
        for args in ((-1,F(1),32,3,9),(2,F(-1),32,3,9),
                     (2,F(1),0,3,9),(2,F(1),32,3,8),(2,0.5,32,3,9)):
            with self.assertRaises(ValueError):y.budget(*args)

    def test_audit_keeps_selection_assumptions_and_curvature_types_explicit(self):
        ledger=json.loads(y.LEDGER.read_text())
        r=y.validate_ledger(ledger)
        self.assertTrue(r['acyclic'])
        self.assertEqual(r['status_counts']['DECLARED_SELECTION'],9)
        self.assertEqual(r['status_counts']['OPEN_EXTENSION'],3)
        by_id={n['id']:n for n in ledger['nodes']}
        self.assertIn('Covariance',by_id['fixed_law_isotropy']['boundary'])
        self.assertIn('physical',by_id['response_metric']['boundary'])

    def test_dependency_graph_and_evidence_mutations_are_refused(self):
        base=json.loads(y.LEDGER.read_text())
        mutations=[]
        m=deepcopy(base);m['nodes'].append(deepcopy(m['nodes'][0]));mutations.append(m)
        m=deepcopy(base);m['nodes'][0]['parents']=['missing'];mutations.append(m)
        m=deepcopy(base);m['nodes'][0]['parents']=['algebra'];mutations.append(m)
        m=deepcopy(base);m['nodes'][0]['sources']=['missing'];mutations.append(m)
        m=deepcopy(base);m['nodes'][0]['boundary']='';mutations.append(m)
        for key in ('fixed_law_isotropy','clock_selector','physical_dictionary'):
            m=deepcopy(base)
            next(n for n in m['nodes'] if n['id']==key)['status']='NATIVE_DERIVED_UNDER_CONTRACT'
            mutations.append(m)
        m=deepcopy(base)
        for n in m['nodes']:
            if 'target_type' in n:n['target_type']='one_curvature'
        mutations.append(m)
        for index,m in enumerate(mutations):
            with self.subTest(mutation=index):
                with self.assertRaises(ValueError):y.validate_ledger(m)

    def test_source_tampering_is_refused(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()

    def test_record_and_pin_tampering_are_refused(self):
        data={'verdict':'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'result.json',Path(folder)/'pin'
            record.write_text(json.dumps(data));pin.write_text(y.y.canonical_sha(data))
            y.check(data,record,pin)
            record.write_text(json.dumps({'verdict':'FAIL'}))
            with self.assertRaises(ValueError):y.check(data,record,pin)
            record.write_text(json.dumps(data));pin.write_text('0'*64)
            with self.assertRaises(ValueError):y.check(data,record,pin)


if __name__=='__main__':unittest.main()
