from common import *

from ve.fp.fp_370 import MyComplex as Complex, MyResidue as Residue

class FPTestCase(NumericTestCase):
    def setUp(self):
        c_id = "1SLG_D"

        cs = load_complexes([c_id], complex_cls = Complex, residue_cls = Residue)
        
        self.complex = list(cs)[0]
        
    def test_fp_length(self):
        self.assertEqual(370, len(self.complex.get_fp("atg").fp_array().tolist()))

class FP480TestCase(NumericTestCase):
    """Test for the 480 data set"""
    def setUp(self):
        c_id = "1SLD_B"
        from ve.config import data480_complex_root
        cs = load_complexes([c_id], directory = data480_complex_root, complex_cls = Complex, residue_cls = Residue)
        
        self.complex = list(cs)[0]
        
    def test_fp_length(self):
        self.assertEqual(370, len(self.complex.get_fp("atb").fp_array().tolist()))


if __name__ == "__main__":
    unittest.main()
