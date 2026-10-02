from fractions import Fraction as F
import itertools
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym45_temporal_blocks as y


class TemporalBlockTests(unittest.TestCase):
    def test_strict_window_sign_and_endpoint(self):
        theta=F(1,1681);alpha=F(1,2)+840*theta;R=(1+1/alpha)/2
        plus=y.block_gate(theta,5,R);minus=y.block_gate(-theta,5,R)
        self.assertIsNotNone(plus)
        self.assertEqual(plus['alpha'],minus['alpha'])
        self.assertGreater(plus['gamma'].lo,0)
        for J in range(1,40):
            self.assertIsNone(y.block_gate(F(1,1680),J,F(1000001,1000000)))

    def test_weighted_equality_is_refused(self):
        self.assertIsNone(y.block_gate(F(1,4096),8,F(64,41)))
        self.assertIsNotNone(y.block_gate(F(1,4096),8,F(64,41)-F(1,1000)))
        for J,R in ((0,F(2)),(2.5,F(2)),(8,F(1))):
            self.assertIsNone(y.block_gate(F(0),J,R))

    def test_conditioned_skeleton_against_enumerated_paths(self):
        q=F(2,3);k=7;t=3
        for left,right in itertools.product(range(2),repeat=2):
            paths,p=y.bridge_law(k,left,right,q)
            actual=[sum(prob for path,prob in zip(paths,p) if path[t-1]==s) for s in range(2)]
            self.assertEqual(actual,y.bridge_transition(left,right,t,k+1-t,q))
        self.assertNotEqual(y.bridge_transition(0,1,3,1,q),y.power_kernel(3,q)[0])

    def test_free_endpoint_coupling_survives_longer_bridges(self):
        data=y.bridge_controls()
        self.assertEqual(data['skeleton_bridge_cases'],18)
        self.assertLessEqual(F(data['segment_cost_upper']),F(data['uniform_budget']))
        self.assertEqual(data['direct_path_marginal_checks'],12)

    def test_tilt_with_biased_future_and_non_grid_parameters(self):
        k=4;delta=F(1,300);eta=delta/(1-delta)
        env=((0,0),(1,0),(1,1),(0,1));changed=list(env);changed[2]=(0,1)
        paths,p=y.bridge_law(k,0,1,F(4,5),delta,env)
        _,q=y.bridge_law(k,0,1,F(4,5),delta,tuple(changed))
        self.assertLessEqual(sum(abs(a-b) for a,b in zip(p,q))/2,4*eta)
        self.assertLessEqual(y.coupling_cost(paths,p,q),4*eta*k)
        _,free=y.bridge_law(k,0,1,F(4,5))
        floor=((1-delta)/(1+delta))**(2*k)
        self.assertGreaterEqual(min(a/b for a,b in zip(p,free)),floor)

    def test_clipped_duplicates_preserve_short_strip_coverage(self):
        for N,ell in ((1,40),(2,17),(9,3)):
            labels=y.blocks(3,N,ell)
            self.assertEqual(len(labels),3*(N+ell-1))
            for i in range(3):
                for t in range(1,N+1):
                    self.assertEqual(sum(col==i and t in B for col,B in labels),ell)
                self.assertEqual(sum(col==i and B[0]==1 for col,B in labels),ell)
                self.assertEqual(sum(col==i and B[-1]==N for col,B in labels),ell)
        with self.assertRaises(ValueError):
            y.blocks(1,0,4)

    def test_joint_updates_preserve_both_marginals_and_contract(self):
        data=y.joint_update_controls()
        self.assertEqual(data['update_labels'],18)
        for row in data['cases']:
            self.assertTrue(row['both_marginals_exact'])
            self.assertLess(F(row['weighted_mismatch_after']),F(row['weighted_mismatch_before']))

    def test_time_limit_transport_converges_at_fixed_width(self):
        theta,J,R=F(1,4096),8,F(4,3);t=F(3,2);m=7
        gamma=y.block_gate(theta,J,R)['gamma'].lo
        ell=y.ex(-theta*(m-1)*t).lo;q=y.ex(-gamma*t).hi
        bounds=[y.gap_gate(ell,q,y.refinement_budget(m,theta,t,t/n)) for n in (1024,4096,16384)]
        self.assertTrue(all(v is not None and q<=v<1 for v in bounds))
        self.assertTrue(all(b<a for a,b in zip(bounds,bounds[1:])))

    def test_refusal_witnesses(self):
        controls=y.negative_controls()
        self.assertEqual(len(controls),12)
        self.assertTrue(all(controls.values()))

    def test_upstream_and_record_tampering(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):
                y.source_checks()
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
