from ve.fp.complex_util.aaindex import AAIndexTrait
from common import *

ComplexClass = make_complex_class(AAIndexTrait)

class AAIndexTestCase(unittest.TestCase):
    def setUp(self):
        self.c = ComplexClass()

    def test_atg_side_fp(self):
        fp = self.c.get_atg_aaindex_fp()
        fp_str =  fp.fp_str(float)
        self.assertTrue("13,0.61,0.00,-0.01", fp_str)
        self.assertTrue("37,0.07,0.00,0.00", fp_str)
        self.assertTrue("62,0.05,1.00,0.11", fp_str)
        self.assertTrue("133,1.32,0.00,0.01", fp_str)

    def test_atb_side_fp(self):
        fp = self.c.get_atb_aaindex_fp()
        fp_str =  fp.fp_str(float)
        self.assertTrue("2,0.05,1.00,0.11", fp_str)
        self.assertTrue("3,0.61,1.00,0.08", fp_str)
        self.assertTrue("4,1.95,0.00,0.00", fp_str)
        self.assertTrue("5,0.00,2.00,0.05", fp_str)

if __name__ == '__main__':
    unittest.main()