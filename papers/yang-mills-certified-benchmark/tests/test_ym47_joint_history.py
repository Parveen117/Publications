from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym47_joint_history as y


class JointHistoryTests(unittest.TestCase):
    def setUp(self):
        self.theta,self.J,self.R,self.rho,self.u=y.CELLS[1]
        self.c=y.constants(*y.CELLS[1])

    def test_entire_window_near_endpoint_and_both_signs(self):
        theta=F(1,1681);J=5;alpha=F(1,2)+840*theta
        rho=1+(1-alpha)/(112*theta*J);weighted=(1+alpha)/2
        R=(1+1/weighted)/2;u=2*rho/(rho+1)
        pos=y.constants(theta,J,R,rho,u);neg=y.constants(-theta,J,R,rho,u)
        self.assertEqual(pos,neg)
        self.assertLess(pos['time_ratio'],1)
        self.assertLess(pos['space_ratio'],1)

    def test_normalization_changes_accumulation_scale(self):
        width,theta,a,t=97,F(1,256),F(1,4096),F(4)
        new=y.normalized_budget(width,theta,a,t)
        raw=y.y44.refinement_budget(width,theta,t,a)
        # Even the old raw bound (before its vacuum division) is worse.
        self.assertLess(new,raw)
        self.assertEqual(y.normalized_budget(width,F(0),a,t),0)
        self.assertEqual(y.normalized_budget(1,theta,a,t),0)

    def test_vacuum_gap_denominator_has_independent_witness(self):
        out=y.vacuum_controls()
        self.assertEqual(len(out['rational_rotated_vacua']),6)
        row=out['rational_rotated_vacua'][2]
        self.assertGreater(F(row['vacuum_distance_squared']),16*F(row['error_squared']))

    def test_independent_three_time_readouts_converge(self):
        out=y.independent_history_controls()
        self.assertEqual(len(out['independent_three_time_histories']),8)
        self.assertEqual(out['exponential_entry_cross_checks'],8)
        for row in out['independent_three_time_histories']:
            self.assertLess(F(row['normalized_operator_error_upper']),
                            F(row['normalized_operator_budget_upper']))

    def test_positive_clock_bound_and_off_grid_times(self):
        out=y.positive_clock_controls()
        self.assertEqual(out['positive_contraction_power_ceilings'],17)
        self.assertEqual(out['off_grid_time_roundings'],27)

    def test_geometric_tail_exact_remainders(self):
        self.assertEqual(y.geometric_controls()['exact_polynomial_and_geometric_tail_identities'],54)
        for q,n in ((F(17,19),57),(F(127,128),513)):
            box=y.power_iv(q,n)
            self.assertLessEqual(box.lo,q**n)
            self.assertGreaterEqual(box.hi,q**n)

    def test_certified_tail_radius_and_monotonicity(self):
        N,bound=y.tail_radius(2,4,F(1),1,F(1),F(1,100000),self.c)
        self.assertLessEqual(bound,F(1,100000))
        self.assertLess(y.tail(2,4,F(1),1,F(1),N+1,self.c),
                        y.tail(2,4,F(1),1,F(1),N,self.c))
        self.assertGreater(y.tail(2,4,F(1),1,F(1),8,self.c),F(1,100000))

    def test_crop_is_least_lawful_radius_and_handles_asymmetry(self):
        a=F(3,1024);D=y.crop_radius(a,self.rho)
        self.assertLessEqual(self.rho**(-D),a*a)
        self.assertGreater(self.rho**(-(D-1)),a*a)
        out=y.cropping_controls()
        self.assertEqual(out['asymmetric_crops'],490)
        self.assertTrue(out['remote_width_1e12_control'])

    def test_joint_modulus_survives_smaller_cutoff(self):
        errors=[]
        for n in (48,96):
            a=F(1,2**n);D=y.crop_radius(a,self.rho)
            error=8*4*self.c['C']*a+y.history_budget(2+2*D,self.theta,a,F(1),1,F(1),self.c)
            error+=3*y.tail(2,4,F(1),1,F(1),D,self.c)
            errors.append(error)
        self.assertLess(errors[1],errors[0])
        self.assertLess(errors[1],F(1,1000000))

    def test_strict_gates_and_time_collisions_are_refused(self):
        for u in (F(1),F(2)):
            with self.assertRaises(ValueError):
                y.constants(self.theta,self.J,self.R,self.rho,u)
        with self.assertRaises(ValueError):
            y.history_budget(2,self.theta,F(1,2),F(1),1,F(1),self.c)
        with self.assertRaises(ValueError):
            y.tail(2,4,F(1),1,F(1),0,self.c)
        with self.assertRaises(ValueError):
            y.normalized_budget(0,self.theta,F(1,16),F(1))
        with self.assertRaises(ValueError):
            y.crop_radius(F(0),self.rho)
        with self.assertRaises(ValueError):
            y.tail_radius(2,4,F(1),1,F(1),F(0),self.c)

    def test_refusal_controls(self):
        controls=y.negative_controls()
        self.assertEqual(len(controls),10)
        self.assertTrue(all(controls.values()))

    def test_source_and_record_tampering(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()
        data={'verdict':'PASS'}
        with tempfile.TemporaryDirectory() as folder:
            record,pin=Path(folder)/'record.json',Path(folder)/'pin'
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
