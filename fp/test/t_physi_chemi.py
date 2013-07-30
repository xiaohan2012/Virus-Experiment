from common import *

from ve.fp.residue_util.geom import GeometryTrait
ResidueClass = make_residue_class(GeometryTrait)

from ve.fp.complex_util.physi_chemi import PhysiChemiTrait
ComplexClass = make_complex_class(PhysiChemiTrait, residue_class = ResidueClass)

class PhysiChemiTestCase(unittest.TestCase):
    def setUp(self):
        self.c = ComplexClass()

    def test_general_case_atg(self):
        fp =  self.c.get_physi_chemi_atg_fp()
        self.assertEqual(len(fp.fp_str().split(",")), 30)
        self.assertEqual(fp.fp_str(),"0.00,0.00,2.65,7.48,10.51,8.57,10.62,10.00,7.87,3.20,0.00,0.00,0.00,0.18,0.48,0.45,0.68,0.32,0.30,0.04,0.00,0.00,1.00,4.00,8.00,9.00,23.00,7.00,11.00,3.00")

    def test_general_case_atb(self):
        fp =  self.c.get_physi_chemi_atb_fp()
        self.assertEqual(len(fp.fp_str().split(",")), 30)
        self.assertEqual(fp.fp_str(), "0.00,0.00,2.62,0.00,0.05,2.07,0.00,0.00,0.00,0.00,0.05,0.00,0.14,0.00,0.04,0.14,0.00,0.00,0.00,0.00,2.00,0.00,3.00,0.00,1.00,1.00,0.00,0.00,0.00,0.00")

if __name__ == "__main__":
    unittest.main()
