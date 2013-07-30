from ve.fp.complex_util.aaindex import AAIndexTrait
from common import *

ComplexClass = make_complex_class(AAIndexTrait)

class AAIndexTestCase(unittest.TestCase):
    def setUp(self):
        self.c = ComplexClass()

    def test_atg_side_fp(self):
        fp = self.c.get_atg_aaindex_fp()
        fp_str =  fp.fp_str(float)
        self.assertTrue("13,0.61,0.00,-0.01" in fp_str)
        self.assertTrue("37,0.07,0.00,0.00" in fp_str)
        self.assertTrue("62,0.05,1.00,0.11" in fp_str)
        self.assertTrue("133,1.32,0.00,0.01" in fp_str)

    def test_atb_side_fp(self):
        fp = self.c.get_atb_aaindex_fp()
        fp_str =  fp.fp_str(float)
        self.assertTrue("2,0.05,1.00,0.11" in fp_str)
        self.assertTrue("3,0.61,1.00,0.08" in fp_str)
        self.assertTrue("4,1.95,0.00,0.00" in fp_str)
        self.assertTrue("5,0.00,2.00,0.05" in fp_str)

    def test_tri_fp(self):
        fp = self.c.get_tri_aaindex_fp()
        fp_str =  fp.fp_str(float)
                
        self.assertTrue("D-79 +D-87 +D-88 ,0.30,0.39,0.03" in fp_str)
        self.assertTrue("D-78 +D-79 +D-88 ,0.46,0.14,0.00" in fp_str)
        self.assertTrue("D-77 +D-79 +D-89 ,0.43,0.06,0.00" in fp_str)
        self.assertTrue("D-77 +D-79 +D-90 ,0.37,0.13,0.00" in fp_str)
        
if __name__ == '__main__':
    unittest.main()