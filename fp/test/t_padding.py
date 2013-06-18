from common import *
from ve.fp.complex_util.padding import OverallSpatialDistribution
    
class OverallDistTestCase(unittest.TestCase):
    
    def setUp(self):
        self.atg_dist, self.atb_dist = OverallSpatialDistribution.from_cache("/home/rxzhu/QiuTianyi/code/ve/data/data237/fp/.cache/overall_residue_distribution_test")
        
    def test_ring_count(self):
        actual = self.atg_dist.ring_count
        expected = 50
        self.assertEqual(actual, expected)

        actual = self.atb_dist.ring_count
        expected = 50
        self.assertEqual(actual, expected)
        
    def test_atg_side(self):
        actual = dict(self.atg_dist)
        expected = {1: 1, 4: 1, 6: 2, 7: 2, 8: 1, 9: 1, 11: 1, 12: 1, 13: 2, 14: 2, 16: 2, 17: 3, 18: 3, 19: 1, 20: 1, 21: 2, 22: 1, 23: 2, 24: 3, 26: 1, 27: 2, 28: 2, 29: 1, 30: 1, 31: 2, 32: 3, 33: 3, 34: 2, 35: 1, 36: 2, 37: 2, 38: 2, 39: 3, 40: 1, 41: 2, 42: 2, 43: 2, 44: 3, 45: 1, 46: 1, 47: 4, 48: 2, 49: 4}
        
        self.assertEqual(actual, expected)

    def test_atb_side(self):
        actual = dict(self.atb_dist)
        expected = {1: 1, 2: 1, 3: 2, 4: 1, 5: 1, 6: 1, 7: 1, 8: 3, 9: 2, 10: 1, 11: 2, 12: 2, 13: 3, 14: 2, 16: 2, 17: 2, 18: 2, 19: 2, 20: 1, 21: 2, 22: 1, 23: 2, 24: 3, 25: 1, 26: 2, 27: 2, 28: 5, 29: 3, 31: 1, 32: 3, 33: 3, 34: 4, 35: 1, 36: 1, 37: 3, 38: 2, 39: 3, 40: 1, 41: 1, 42: 2, 43: 2, 44: 2, 46: 1, 47: 2, 48: 2, 49: 4}
        
        self.assertEqual(actual, expected)


from ve.fp.complex_util.res_spat_dist import ResiduePositionDistribution

atg_dist, atb_dist = OverallSpatialDistribution.from_cache()
        
from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait
from ve.fp.complex_util.triangle import ResidueTriangleTrait

from ve.fp.fp_80 import Residue as ResidueClass

TempClass = make_complex_class(SplitCylinderViaComplexPlaneTrait, residue_class = ResidueClass)

from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait
class ComplexClass(TempClass, ResidueTriangleTrait, ResidueSpatialDistributionTrait):
    pass

class PaddingTestCase(unittest.TestCase):
    def setUp(self):
        self.atg_dist = atg_dist

        self.atb_dist = atb_dist
        
        self.c = ComplexClass()

    def test_len_atg_side(self):
        """test the length of finger print of the atg side"""
        from ve.fp.complex_util.padding import PaddedComplexFingerPrint
        from ve.fp.fp import HeadlessFingerprint
        fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                 targets=[(0,self.c.get_triangles())],
                                                 fps = PaddedComplexFingerPrint())
        _, res_dist = self.c.get_atg_res_spat_dist()
        actual = len(fps.fp_str(self.atg_dist, res_dist, number_type=int).split(","))
        expected = sum(self.atg_dist.values()) * fps.get_res_fp_length()

        self.assertEqual(actual, expected)

    def test_len_atb_side(self):
        """test the length of finger print of the atb side"""
        from ve.fp.complex_util.padding import PaddedComplexFingerPrint
        from ve.fp.fp import HeadlessFingerprint
        fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                  targets=[(0,self.c.atb.residues)], 
                                                  fps = PaddedComplexFingerPrint())
        _, res_dist = self.c.get_atb_res_spat_dist()
        actual = len(fps.fp_str(self.atb_dist, res_dist, number_type=int).split(","))
        expected = sum(self.atb_dist.values()) * fps.get_res_fp_length()

        self.assertEqual(actual, expected)

from ve.fp.complex_util.padding import PaddedComplexFingerPrint
class PaddedFingerPrintTest(unittest.TestCase):
    """
    pickle related method of PaddedFingerPrint
    """
    def setUp(self):
        self.c = ComplexClass()

        self.fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                  targets=[(0,self.c.atb.residues)], 
                                                  fps = PaddedComplexFingerPrint())
    def test_to_pickle_mapping_len(self):
        """
        mapping length of to_pickle
        """
        p = self.fps.to_pickable()
        self.assertEqual(len(p.mapping), 7)

    def test_from_pickable(self):
        """
        test whether padded fingerprint loaded from cache works
        """
        p = self.fps.to_pickable()
        fps = PaddedComplexFingerPrint.from_pickable(p, self.c)
        
        self.assertEqual(self.fps, fps)

class FromCachePaddedFingerPrintTest(unittest.TestCase):
    pass
    
if __name__ == '__main__':
    unittest.main()
        
