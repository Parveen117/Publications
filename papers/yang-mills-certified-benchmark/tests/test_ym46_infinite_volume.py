from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym46_infinite_volume as y


class InfiniteVolumeTests(unittest.TestCase):
    def test_spatial_gate_equality_and_sign(self):
        theta,J,R=F(1,4096),8,F(4,3)
        self.assertIsNotNone(y.y45.block_gate(theta,J,R))
        self.assertIsNone(y.joint_gate(theta,J,R,F(2)))
        self.assertIsNotNone(y.joint_gate(theta,J,R,F(199,100)))
        pos=y.joint_gate(theta,J,R,F(3,2));neg=y.joint_gate(-theta,J,R,F(3,2))
        self.assertEqual(pos['margin'],F(7,96))
        self.assertEqual(pos['margin'],neg['margin'])

    def test_entire_open_window_construction_near_endpoint(self):
        theta=F(1,1681);J=5;alpha=F(1,2)+840*theta
        rho=1+(1-alpha)/(112*theta*J);A=(1+alpha)/2;R=(1+1/A)/2
        gate=y.joint_gate(theta,J,R,rho)
        self.assertIsNotNone(gate)
        self.assertEqual(gate['beta'],(1+A)/2)
        self.assertGreater(gate['gamma'].lo,0)

    def test_constructive_cauchy_modulus_off_grid(self):
        a=F(1,7);theta,J,R,rho=y.CELLS[1];tol=F(1,1000)
        box=y.certified_box(a,theta,J,R,rho,tol,support=3)
        bound=y.box_tail(a,theta,J,R,rho,box['space_distance'],box['time_distance_steps'],3)
        self.assertLessEqual(bound,tol)
        self.assertGreater(y.box_tail(a,theta,J,R,rho,1,1,3),tol)

    def test_target_weights_and_all_boundary_faces(self):
        out=y.weighted_incidence_controls()
        self.assertEqual(out['layouts'],27)
        self.assertEqual(out['interior_controls'],216)
        self.assertEqual(out['boundary_controls'],306)
        self.assertTrue(out['multi_target_weights_checked'])

    def test_local_specification_matches_full_conditional(self):
        out=y.specification_controls()
        self.assertEqual(out['conditioned_exterior_configs'],16)
        self.assertTrue(out['local_specification_matches_full_law'])
        self.assertTrue(out['invariance_and_nested_consistency'])

    def test_side_boundaries_and_half_edges_change_readout(self):
        states,p=y.finite_strip(3,2,0,1,False)
        _,q=y.finite_strip(3,2,1,1,False)
        _,half=y.finite_strip(3,2,0,1,True)
        mean=lambda law:sum(prob*x[1] for x,prob in zip(states,law))
        self.assertNotEqual(mean(p),mean(q))
        self.assertNotEqual(mean(p),mean(half))
        self.assertEqual(sum(p),1)

    def test_actual_square_ordering_and_midpoint_paths(self):
        out=y.half_layer_controls()
        self.assertEqual(out['normalized_bridge_and_entry_checks'],84)
        self.assertEqual(out['path_observable_identities'],40)
        self.assertEqual(out['half_edge_conditioned_skeletons'],12)

    def test_positive_time_forms_reflection_and_all_source_ceiling(self):
        out=y.positive_time_controls()
        self.assertEqual(out['reflection_identities'],12)
        self.assertEqual(out['mixed_source_power_controls'],36)
        self.assertEqual(out['all_source_ceiling_inertia'],[2,0,2])

    def test_cutoff_prefactor_cannot_be_silently_uniform(self):
        theta,J,R,rho=y.CELLS[1]
        coarse=y.budgets(F(1),theta,J,R,rho)
        fine=y.budgets(F(1,64),theta,J,R,rho)
        self.assertGreater(fine['space'],coarse['space'])
        self.assertGreater(fine['time'],coarse['time'])
        self.assertEqual(y.budgets(F(1,8),F(0),5,F(3,2),F(2))['space'],0)

    def test_refusals_and_invalid_parameters(self):
        controls=y.negative_controls()
        self.assertEqual(len(controls),12)
        self.assertTrue(all(controls.values()))
        theta,J,R,rho=y.CELLS[1]
        for a in (F(0),F(-1),F(2)):
            with self.assertRaises(ValueError):
                y.budgets(a,theta,J,R,rho)
        with self.assertRaises(ValueError):
            y.certified_box(F(1),theta,J,R,rho,F(0))
        with self.assertRaises(ValueError):
            y.box_tail(F(1),theta,J,R,rho,0,1)

    def test_upstream_source_tampering(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()

    def test_record_and_pin_tampering(self):
        data={'verdict':'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'r.json',Path(folder)/'pin'
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
