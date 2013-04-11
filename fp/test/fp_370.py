from common import *

from ve.fp.fp_370 import MyComplex as Complex, MyResidue as Residue

class FPTestCase(NumericTestCase):
    def setUp(self):
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)

        c_id = "1SLG_D"

        self.complex = Complex(c_id, atg, atb)

    def test_fp_length(self):
        print self.complex.get_fp().fp_str()
        self.assertEqual(370, len(self.complex.get_fp().fp_array().tolist()))


if __name__ == "__main__":
    unittest.main()
