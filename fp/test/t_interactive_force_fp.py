from ve.fp.complex_util.interactive_force_fp import InteractiveForceTrait
from ve.fp.residue_util.geom import GeometryTrait

from common import *

#Residue class for interactive force fingerprint
IFResidue = make_residue_class(GeometryTrait)
        
from common import make_complex_class
ComplexClass = make_complex_class(InteractiveForceTrait, residue_class = IFResidue)

class ResidueInteractiveForceTestCase(NumericTestCase):
    """Interactive Force Fingerprint test of the residue fp"""

    def setUp(self):
        self.c = ComplexClass()
        
    def test_atg_based_fp(self):
        """residue-scope antigen as receptor"""
        fp = self.c.gen_if_residue_fp_atg()

        #the length including the residue number should be 16
        self.assertEqual(len(fp.fp_str(int).split("\n")[0].split(",")), 16)

    def test_atb_based_fp(self):
        """residue-scope antibody as receptor"""
        fp = self.c.gen_if_residue_fp_atb()

        #the length including the residue number should be 16
        self.assertEqual(len(fp.fp_str(int).split("\n")[0].split(",")), 16)
        

class ComplexInteractiveForceTestCase(NumericTestCase):
    """Interactive Force Fingerprint test of the complex fp"""
    def setUp(self):
        self.c = ComplexClass()
        
    def test_atg_based_complex_fp(self):
        """complex-scope antigen as receptor"""
        fp = self.c.gen_if_complex_fp_atg()
        
        #fp length should be 150
        self.assertEqual(len(fp.fp_array().tolist()), 150)
        self.assertArrayEqual(fp.fp_array(), [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 6, 6, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0,12, 0, 0,12,12, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0,14, 0, 0,14,10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,15, 0, 0,15,12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_atb_based_complex_fp(self):
        """complex-scope antibody as receptor"""
        fp = self.c.gen_if_complex_fp_atb()
        
        #fp length should be 150
        self.assertEqual(len(fp.fp_array().tolist()), 150)
        self.assertArrayEqual(fp.fp_array(), [1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

if __name__ == "__main__":
    unittest.main()
