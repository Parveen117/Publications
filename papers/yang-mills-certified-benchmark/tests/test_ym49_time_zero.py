from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym49_time_zero as y


class TimeZeroTests(unittest.TestCase):
    def test_time_zero_fold_matches_independent_full_paths(self):
        result=y.path_controls()
        self.assertEqual(result['independent_time_zero_path_checks'],48)
        self.assertTrue(result['null_time_zero_action_is_zero'])
        self.assertEqual(result['future_multiplier_null_violation'],'9/16')

    def test_complex_products_daggers_and_operator_bounds(self):
        result=y.observable_controls()
        self.assertEqual(result['complex_product_checks'],9)
        self.assertEqual(result['complex_dagger_checks'],3)
        self.assertEqual(result['complex_norm_bounds'],3)
        # Multiplication by i has square -I and adjoint -i, over exact rationals.
        i=y.complex_multiplier((0,)*4,(1,)*4)
        self.assertEqual(y.multiply(i,i),y.scale(y.identity(8),-1))
        self.assertEqual(y.transpose(i),y.scale(i,-1))

    def test_multiplier_form_bound_is_all_source_psd(self):
        self.assertEqual(y.observable_controls()['all_source_form_multiplier_bounds'],3)
        H=y.y48.generator();values=tuple(map(F,(1,-3,2,4)));M=y.diag(values)
        C=y.jump_gradient_ceiling(values)
        difference=y.add(y.add(y.scale(H,32),y.scale(y.I4,2*C)),
                         y.multiply(y.multiply(M,H),M),-1)
        self.assertTrue(y.psd(difference))

    def test_nonconstant_vacuum_weight_product_energy(self):
        self.assertEqual(y.weighted_form_controls()[
            'nonconstant_weight_product_energy_and_bound_checks'],18)
        x=y.y48.COORD[0]
        actual=y.y48.integral(y.y48.grad_square(y.y48.mul(x,x)))
        no_cross=2*y.y48.integral(y.y48.mul(y.y48.mul(x,x),y.y48.grad_square(x)))
        self.assertEqual(actual,F(1,8))
        self.assertEqual(no_cross,F(1,16))

    def test_exact_collisions_and_cauchy_bounds(self):
        result=y.collision_controls()
        self.assertEqual(result['exact_collision_and_cauchy_checks'],35)
        self.assertEqual([r['n'] for r in result['samples']],[8,32,128,512])
        for row in result['samples']:
            self.assertLessEqual(F(row['source_error_squared_max']),
                                 F(row['collision_bound_upper'])**2)

    def test_collision_errors_shrink_on_independent_samples(self):
        samples=y.collision_controls()['samples']
        for previous,current in zip(samples,samples[1:]):
            self.assertLess(F(current['source_error_squared_max']),
                            F(previous['source_error_squared_max']))
            self.assertLess(F(current['collision_bound_upper']),
                            F(previous['collision_bound_upper']))

    def test_two_time_defect_equals_leakage_square(self):
        result=y.closure_controls()
        self.assertEqual(result['mixed_time_defect_checks'],9)
        self.assertEqual(result['leakage_norm_squared'],'7/288')
        self.assertEqual(result['nonclosing_defect_rank'],1)
        self.assertEqual(result['actual_chain_row_closure'],'OPEN')

    def test_vacuum_probe_and_closing_observer_are_distinct_controls(self):
        U,P,Q,A,B,C,D=y.memory_blocks()
        self.assertEqual(y.apply(C,y.VAC),[F(0)]*4)
        self.assertNotEqual(C,y.scale(y.I4,0))
        _,_,_,_,_,closed_C,_=y.memory_blocks(True)
        self.assertEqual(closed_C,y.scale(y.I4,0))
        self.assertNotEqual(P,y.row_cut(True))

    def test_exact_memory_recursion_and_uniform_tail_bounds(self):
        result=y.memory_controls()
        self.assertEqual(result['positive_memory_coefficient_bounds'],12)
        self.assertEqual(result['independent_recursion_checks'],36)
        self.assertEqual(result['memory_tail_checks'],108)
        self.assertEqual(result['full_inverse_schur_checks'],3)
        self.assertEqual(result['schur_tail_checks'],12)

    def test_hidden_initial_state_changes_same_visible_start(self):
        U,P,Q,A,B,C,D=y.memory_blocks()
        hidden=y.apply(Q,(F(1),F(0),F(0),F(0)))
        self.assertEqual(y.apply(P,hidden),[F(0)]*4)
        response=y.apply(P,y.apply(U,hidden))
        self.assertEqual(response,y.apply(B,hidden))
        self.assertGreater(y.dot(response,response),0)

    def test_full_schur_inverse_and_geometric_domain(self):
        U,P,Q,A,B,C,D=y.memory_blocks();z=F(5,4)
        hidden_inverse=y.inverse(y.add(y.scale(y.I4,z),D,-1))
        sigma=y.multiply(y.multiply(B,hidden_inverse),C)
        augmented=y.add(y.add(y.add(y.scale(P,z),A,-1),sigma,-1),Q)
        source=y.apply(P,(F(2),F(-1),F(3),F(4)))
        direct=y.apply(P,y.apply(y.inverse(y.add(y.scale(y.I4,z),U,-1)),source))
        reduced=y.apply(y.inverse(augmented),source)
        self.assertEqual(direct,reduced)
        self.assertGreater(y.schur_tail_bound(F(7,288),F(1,2),z,4),0)

    def test_row_algebra_invariance_does_not_imply_time_invariance(self):
        U,P,Q,A,B,C,D=y.memory_blocks()
        M=y.diag((1,1,1,-2))
        self.assertEqual(y.multiply(P,M),y.multiply(M,P))
        self.assertNotEqual(y.multiply(P,U),y.multiply(U,P))
        self.assertEqual(y.closure_controls()['observable_time_cyclic_gram_inertia'],[4,0,0])

    def test_future_null_multiplier_and_other_refusals(self):
        controls=y.negative_controls()
        self.assertEqual(len(controls),12)
        self.assertTrue(all(controls.values()))

    def test_noncommuting_time_witness_and_inherited_window(self):
        result=y.noncommuting_controls()
        self.assertEqual(result['exact_fixture_noncommutator_checks'],4)
        self.assertEqual(len(result['inherited_chain_parameter_bounds']),4)
        for row in result['inherited_chain_parameter_bounds']:
            self.assertGreater(F(row['chain_commutator_norm_lower']),0)
            self.assertLess(F(row['chain_commutator_norm_lower'])**2,F(3,8))

    def test_invalid_collision_memory_and_resolvent_bounds(self):
        for delta,epsilon,f,C,M in ((F(0),F(1),F(1),F(1),F(1)),
                (F(1,2),F(1),F(1),F(1),F(1)),
                (F(1,8),F(1),F(1),F(-1),F(1)),
                (F(1,8),F(1),F(1),F(1),F(-1))):
            with self.assertRaises(ValueError):
                y.insertion_budget(delta,epsilon,f,C,M)
        for ell,r,N in ((F(-1),F(1,2),3),(F(1),F(1),3),(F(1),F(1,2),-1),
                        (F(1),F(1,2),F(1,2))):
            with self.assertRaises(ValueError):
                y.memory_tail_bound(ell,r,N)
        with self.assertRaises(ValueError):
            y.schur_tail_bound(F(1),F(1,2),F(1),3)
        with self.assertRaises(ValueError):
            y.inverse([[F(0)]])

    def test_source_tampering_is_refused(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()

    def test_record_and_pin_tampering_are_refused(self):
        data={'verdict':'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'result.json',Path(folder)/'pin'
            record.write_text(json.dumps(data));pin.write_text(y.canonical_sha(data))
            y.check(data,record,pin)
            record.write_text(json.dumps({'verdict':'FAIL'}))
            with self.assertRaises(ValueError):
                y.check(data,record,pin)
            record.write_text(json.dumps(data));pin.write_text('0'*64)
            with self.assertRaises(ValueError):
                y.check(data,record,pin)


if __name__=='__main__':
    unittest.main()
