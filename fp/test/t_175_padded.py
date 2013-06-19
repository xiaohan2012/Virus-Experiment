from common import *

from ve.util.load_pdb import load_pdb_struct

from ve.fp.fp_175_padded import ComplexPlaneBasedComplex, ResiduePlaneBasedComplex
from ve.fp.complex_util.padding import OverallSpatialDistribution

class ComplexPlaneBasedComplexTestCase(unittest.TestCase):
    
    def setUp(self):
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        c_id = "1SLG_D"
        self.c = ComplexPlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)

    def test_fp_length_atg_as_rec(self):
        """
        in case antigen is set as the receptor, test the fp length
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str()
        actual = len(fp_str.split(","))
        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80
        
        self.assertEqual(actual, expected)

    def test_fp_length_atb_as_rec(self):
        """
        in case antibody is set as the receptor, test the fp length
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(atg_as_receptor = False)
        actual = len(fp_str.split(","))
        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)


class ResiduePlaneBasedComplexTestCase(unittest.TestCase):
    
    def setUp(self):
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        c_id = "1SLG_D"
        self.c = ResiduePlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)

    def test_fp_length_atg_as_rec(self):
        """
        in case antigen is set as the receptor, test the fp length
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str()
        actual = len(fp_str.split(","))
        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80

        self.assertEqual(actual, expected)

    def test_fp_length_atb_as_rec(self):
        """
        in case antibody is set as the receptor, test the fp length
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(atg_as_receptor = False)
        actual = len(fp_str.split(","))
        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)
        
if __name__ == '__main__':
    unittest.main()