import unittest
from ve.fp.fp import *

from common import *

class BaseResidueFingerprintTest(NumericTestCase):
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
        
    def test_fp_array(self):
        fp  = BaseResidueFingerprint(None, 5)
        fp[0] = 1
        fp[1] = 2
        self.assertArrayEqual(fp.fp_array(), [1,2,0,0,0])
        
    def test_append(self):
        """test the append method"""
        fp1 = BaseResidueFingerprint(None, 4)#4 0 0 0
        fp1[0] = 4
        
        fp2 = BaseResidueFingerprint(None, 6)#6 0 0 0 0 0
        fp2[0] = 6
        
        fp3 = fp1.append(fp2)
        
        #ensure fp1 and fp2 are unchanged
        self.assertArrayEqual(fp1.fp_array(), [4,0,0,0])
        self.assertArrayEqual(fp2.fp_array(), [6,0,0,0,0,0])
        
        #fp3 changed accordingly
        self.assertArrayEqual(fp3.fp_array(), [4,0,0,0,6,0,0,0,0,0])
        self.assertEqual(fp3.bitlength, 10)


    def test_add(self):
        """test the add operator"""
        fp1  = BaseResidueFingerprint(None, 3)
        fp1[0] = 1
        fp1[1] = 1
        
        fp2  = BaseResidueFingerprint(None, 3)
        fp2[1] = 1
        fp2[2] = 3
        
        fp3 = fp1 + fp2
        self.assertArrayEqual(fp3.fp_array(), [1,2,3])
        self.assertEqual(fp3.bitlength, 3)
                
        #ensure fp1 and fp2 unchanged
        self.assertArrayEqual(fp1.fp_array(), [1,1,0])
        self.assertArrayEqual(fp2.fp_array(), [0,1,3])
        
        #test  +=
        fp1 += fp2
        self.assertArrayEqual(fp1.fp_array(), [1,2,3])

    def test_to_pickable(self):
        """test for to_pickable method of ResidueFingerPrint"""
        c = load_base_complex()
        res = c.residues[0]
        
        fp  = BaseResidueFingerprint(res, 5)
        fp[0] = 1.2
        fp[2] = 2.3
        fp[4] = 4.56
        #resnum, bit_len, pickled_d
        p= fp.to_pickable()
        self.assertEqual(p.res_id, "D-13 ")
        self.assertEqual(p.bitlength, 5)

        from customcollections import OrderedDefaultDict
        d = OrderedDefaultDict(float)
        d[0] = 1.2
        d[2] = 2.3
        d[4] = 4.56

        self.assertEqual(p.mapping, d)
        
    def test_from_pickable(self):
        from ve.fp.fp import BaseResidueFingerprint as FPClass
        c = load_base_complex()
        res = c.residues[0]
        
        fp  = BaseResidueFingerprint(res, 5)
        fp[0] = 1.2
        fp[2] = 2.3
        fp[4] = 4.56

        p = fp.to_pickable()
        
        identical_fp = FPClass.from_pickable(p, c)

        self.assertEqual(fp, identical_fp)

class HeadlessFingerprintTest(unittest.TestCase):
    def test_general_case(self):
        fp  = HeadlessFingerprint(5)
        fp[0] = 1.2
        fp[2] = 2.3
        fp[4] = 4.56
        expected = "1.20,0.00,2.30,0.00,4.56"
        self.assertEqual(fp.fp_body(), expected)
        
    def test_zero_bitlength(self):
        """bit length equals to zero"""
        fp  = HeadlessFingerprint(0)
        expected = ""
        self.assertEqual(fp.fp_body(), "")

    def test_fp_body_unset(self):
        """all fp elements remain as 0"""
        fp  = HeadlessFingerprint(5)
        expected = ""
        self.assertEqual(fp.fp_body(), "0.00,0.00,0.00,0.00,0.00")

class BaseComplexFingerprintTestCase(NumericTestCase):
    def test_bit_length(self):
        #the residue fingerprint
        res_fp  = HeadlessFingerprint(15)
        
        #the complex fp with the first element set
        fp = BaseComplexFingerprint()
        fp[0] = res_fp

        self.assertEqual(fp.get_bitlength(), 15)

    def prepare_fps(self):
        """
        (TestCase) -> (Complex, BaseComplexFingerPrint)
        
        Prepare the complex fps
        """
        c = load_base_complex()
        fps = BaseComplexFingerprint()
        for r in c.residues[:5]:
            fp  = BaseResidueFingerprint(r, 5)
            fp[0] = 1.2
            fp[2] = 2.3
            fp[4] = 4.56
            fps[r] = fp
        return c, fps
        
    def test_to_pickable(self):
        """test the to_pickable method of complex fp"""
        _, fps = self.prepare_fps()
        p = fps.to_pickable()
        self.assertEqual(p.mapping["D-13 "].mapping.values(), [1.2, 2.3, 4.56])
        self.assertEqual(len(p.mapping), 5)

    def test_from_pickable(self):
        """test the from_pickable method of complex fp"""
        c, fps = self.prepare_fps()
        p = fps.to_pickable()
        new_p = BaseComplexFingerprint.from_pickable(p, c)

        self.assertEqual(fps, new_p)
        
if __name__ == "__main__":
    unittest.main()
