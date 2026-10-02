from fractions import Fraction as F
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'certificates'))
import ym43_native_time_transfer as y


class NativeTimeBoundTests(unittest.TestCase):
    def test_full_temporal_budget_is_required(self):
        self.assertIsNone(y.gate(F(1,8),F(1,8),F(1,16)))
        self.assertLess(2*F(1,8)+F(1,8)*F(1,16),1)
        self.assertEqual(y.gate(F(1,16),F(1,32),F(1,4)),F(33,128))

    def test_exact_overlap_marginals_including_endpoints(self):
        for p,q in (([F(0),F(1)],[F(1),F(0)]),([F(1,3),F(2,3)],[F(1,3),F(2,3)]),
                    ([F(2,3),F(1,3)],[F(1,5),F(4,5)])):
            J=y.coupling(p,q)
            self.assertEqual([sum(row) for row in J],p)
            self.assertEqual([sum(row[j] for row in J) for j in range(2)],q)
            self.assertEqual(J[0][1]+J[1][0],sum(abs(a-b) for a,b in zip(p,q))/2)

    def test_invalid_kernel_refused_by_iteration(self):
        for T in ([[F(1),F(0)],[F(0),F(1)]],[[F(1),F(2)],[F(3),F(1)]]):
            with self.assertRaises(ValueError):y.vacuum_controls(T)

    def test_square_order_intertwiner_identity(self):
        X=[[F(1),F(1,3)],[F(1,5),F(1,2)]]
        Xd=[list(row) for row in zip(*X)]
        S,T=y.multiply(X,Xd),y.multiply(Xd,X)
        self.assertNotEqual(S,T)
        self.assertEqual(y.multiply(S,S),y.multiply(y.multiply(X,T),Xd))
        self.assertEqual(y.multiply(S,X),y.multiply(X,T))

    def test_directed_rounding(self):
        for x in (F(1,3),F(-1,3),F(17,4)):
            self.assertLessEqual(y.down(x),x)
            self.assertGreaterEqual(y.up(x),x)
            self.assertLessEqual(F(y.decimals(x,False)),x)
            self.assertGreaterEqual(F(y.decimals(x)),x)

    def test_all_refusal_groups(self):
        c=y.negative_controls()
        self.assertEqual(len(c),8)
        self.assertTrue(all(c.values()))

    def test_source_tamper(self):
        with patch.object(y,'digest',return_value='bad'):
            with self.assertRaisesRegex(ValueError,'upstream source changed'):y.source_checks()

    def test_certificate_and_pin_tamper(self):
        data={'verdict':'PASS','q':'1/4'}
        with tempfile.TemporaryDirectory() as folder:
            result,pin=Path(folder)/'r.json',Path(folder)/'p.txt'
            result.write_text(json.dumps(data));pin.write_text(y.canonical_sha(data))
            y.check(data,result,pin)
            result.write_text(json.dumps({**data,'q':'1/8'}))
            with self.assertRaises(ValueError):y.check(data,result,pin)
            result.write_text(json.dumps(data));pin.write_text('0'*64)
            with self.assertRaises(ValueError):y.check(data,result,pin)


if __name__=='__main__':unittest.main()
