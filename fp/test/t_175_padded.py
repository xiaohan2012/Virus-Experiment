from common import *

from ve.util.load_pdb import load_pdb_struct

from ve.fp.fp_175_padded import ComplexPlaneBasedComplex, ResiduePlaneBasedComplex
from ve.fp.complex_util.padding import OverallSpatialDistribution

class ComplexPlaneBasedComplexTestCase(unittest.TestCase):
    """using cache and without cache cases altogether"""
    def setUp(self):
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        #c_id = "1SLG_D"
        c_id = "2NLJ_C"
        self.c = ComplexPlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)
        
    def est_fp_length_use_tri_atg_as_rec(self):
        """
        test the fp length
        in case:
        1, antigen is set as the receptor
        2, iterate through antigen residue triangles
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = True, atg_as_receptor = True, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = True, atg_as_receptor = True, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80
        
        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)

    def test_fp_length_use_res_atg_as_rec(self):
        """
        test the fp length
        in case:
        1, antigen is set as the receptor
        2, iterate through antigen residues
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = False, atg_as_receptor = True, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = False, atg_as_receptor = True, use_cache = True)
        actual_using_cache = len(fp_str.split(","))
        
        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80
        
        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)
        
    def est_fp_length_use_res_atb_as_rec(self):
        """
        test the fp length
        in case:
        1, antibody is set as the receptor
        2, iterate through antigen residues
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = False, atg_as_receptor = False, use_cache = False)
        actual = len(fp_str.split(","))
        
        fp_str_using_cache = self.c.gen_fp_str(use_tri = False, atg_as_receptor = False, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)

    def est_fp_length_use_tri_atb_as_rec(self):
        """
        test the fp length
        in case:
        1, antibody is set as the receptor
        2, iterate through antigen residue triangles
        """

        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = True, atg_as_receptor = False, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = True, atg_as_receptor = False, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)


class ResiduePlaneBasedComplexTestCase(unittest.TestCase):
    
    def setUp(self):
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        c_id = "1SLG_D"
        self.c = ResiduePlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)

    def test_fp_length_use_tri_atg_as_rec(self):
        """
        test the fp length
        in case:
        1, antigen is set as the receptor
        2, iterate through antigen residue triangles
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = True, atg_as_receptor = True, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = True, atg_as_receptor = True, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80
        
        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)

    def test_fp_length_use_res_atg_as_rec(self):
        """
        test the fp length
        in case:
        1, antigen is set as the receptor
        2, iterate through antigen residues
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = False, atg_as_receptor = True, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = False, atg_as_receptor = True, use_cache = True)
        actual_using_cache = len(fp_str.split(","))
        
        expected = sum(overall_atg_dist.values()) * (80+15) + sum(overall_atb_dist.values()) * 80
        
        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)
        
    def test_fp_length_use_res_atb_as_rec(self):
        """
        test the fp length
        in case:
        1, antibody is set as the receptor
        2, iterate through antigen residues
        """
        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = False, atg_as_receptor = False, use_cache = False)
        actual = len(fp_str.split(","))
        
        fp_str_using_cache = self.c.gen_fp_str(use_tri = False, atg_as_receptor = False, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)

    def test_fp_length_use_tri_atb_as_rec(self):
        """
        test the fp length
        in case:
        1, antibody is set as the receptor
        2, iterate through antigen residue triangles
        """

        overall_atg_dist,overall_atb_dist =  OverallSpatialDistribution.from_cache()
        fp_str = self.c.gen_fp_str(use_tri = True, atg_as_receptor = False, use_cache = False)
        actual = len(fp_str.split(","))

        fp_str_using_cache = self.c.gen_fp_str(use_tri = True, atg_as_receptor = False, use_cache = True)
        actual_using_cache = len(fp_str.split(","))

        expected = sum(overall_atg_dist.values()) * 80 + sum(overall_atb_dist.values()) * (80 + 15)

        self.assertEqual(actual, expected)
        self.assertEqual(actual_using_cache, expected)

if __name__ == '__main__':
    unittest.main()