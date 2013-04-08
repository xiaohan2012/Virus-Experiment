import unittest
from ve.fp.fp import *

class BaseResidueFingerprintTest(unittest.TestCase):
    def test_general_case(self):
        fp  = BaseResidueFingerprint(None, 5)
        fp[0] = 1.2
        fp[2] = 2.3
        fp[4] = 4.56
        expected = "1.20,0.00,2.30,0.00,4.56"
        self.assertEqual(fp.fp_body(), expected)
        
    def test_zero_bitlength(self):
        """bit length equals to zero"""
        fp  = BaseResidueFingerprint(None, 0)
        expected = ""
        self.assertEqual(fp.fp_body(), "")

    def test_fp_body_unset(self):
        """all fp elements remain as 0"""
        fp  = BaseResidueFingerprint(None, 5)
        expected = ""
        self.assertEqual(fp.fp_body(), "0.00,0.00,0.00,0.00,0.00")

class GeometricCenterBasedComplexFingerprintTest(unittest.TestCase):
    def test_general_case(self):
        fp  = GCBCFingerprint(5)
        fp[0] = 1.2
        fp[2] = 2.3
        fp[4] = 4.56
        expected = "1.20,0.00,2.30,0.00,4.56"
        self.assertEqual(fp.fp_body(), expected)
        
    def test_zero_bitlength(self):
        """bit length equals to zero"""
        fp  = GCBCFingerprint(0)
        expected = ""
        self.assertEqual(fp.fp_body(), "")

    def test_fp_body_unset(self):
        """all fp elements remain as 0"""
        fp  = GCBCFingerprint(5)
        expected = ""
        self.assertEqual(fp.fp_body(), "0.00,0.00,0.00,0.00,0.00")


if __name__ == "__main__":
    unittest.main()
