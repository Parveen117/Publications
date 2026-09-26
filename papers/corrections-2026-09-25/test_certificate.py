import hashlib, unittest
from fractions import Fraction as F
import certificate as C
class CorrectionTests(unittest.TestCase):
    def test_controls_and_pin(self):
        p=C.build(); self.assertEqual(p['status'],'PASS_FINITE_CONTROLS')
        self.assertTrue(all(p['checks'].values()))
        b=C.canonical(p);self.assertEqual(b,(C.HERE/'CERTIFICATE.json').read_bytes())
        self.assertEqual(hashlib.sha256(b).hexdigest(),(C.HERE/'EXPECTED.sha256').read_text().strip())
    def test_bad_error_budgets_refused(self):
        with self.assertRaises(ValueError):C.spectral_upper(0,-1)
        with self.assertRaises(ValueError):C.pluecker_upper(1,-1)
        with self.assertRaises(ValueError):C.pluecker_upper(-1,1)
    def test_missing_zero_extension_changes_verdict(self):
        self.assertLessEqual(-2+2,1)
        self.assertGreater(C.spectral_upper(F(-2),F(2)),1)
if __name__=='__main__':unittest.main()
