from common import *
from ve.fp.complex_util.padding import OverallSpatialDistribution
    
class CylinderStatTestCase(unittest.TestCase):
    
    def setUp(self):
        self.atg_dist, self.atb_dist = OverallSpatialDistribution.from_cache()
        
    def test_atg_side(self):
        actual = dict(self.atg_dist)
        expected = {58: 15, 78: 13, 79: 12, 47: 11, 56: 11, 57: 11, 66: 11, 67: 11, 68: 10, 76: 10, 77: 10, 55: 9, 59: 9, 65: 9, 69: 9, 75: 9, 37: 8, 44: 8, 45: 8, 46: 8, 49: 8, 54: 8, 63: 8, 64: 8, 74: 8, 27: 7, 36: 7, 38: 7, 48: 7, 53: 7, 33: 6, 34: 6, 35: 6, 39: 6, 43: 6, 73: 6, 9: 5, 24: 5, 25: 5, 26: 5, 28: 5, 32: 5, 42: 5, 52: 5, 62: 5, 72: 5, 2: 4, 4: 4, 5: 4, 6: 4, 8: 4, 12: 4, 13: 4, 15: 4, 16: 4, 17: 4, 19: 4, 22: 4, 23: 4, 29: 4, 31: 4, 41: 4, 51: 4, 61: 4, 71: 4, 1: 3, 3: 3, 14: 3, 18: 3, 21: 3, 7: 2, 20: 2, 30: 2, 40: 2, 50: 2, 60: 2, 70: 2, 10: 1, 11: 1}
        self.assertEqual(actual, expected)

    def test_atb_side(self):
        actual = dict(self.atb_dist)
        expected = {59: 14, 38: 12, 48: 12, 49: 12, 58: 12, 68: 12, 47: 11, 28: 10, 57: 10, 27: 9, 35: 9, 37: 9, 39: 9, 44: 9, 45: 9, 46: 9, 67: 9, 69: 9, 75: 9, 77: 9, 17: 8, 26: 8, 34: 8, 36: 8, 53: 8, 54: 8, 55: 8, 56: 8, 64: 8, 65: 8, 66: 8, 78: 8, 16: 7, 18: 7, 24: 7, 29: 7, 43: 7, 73: 7, 74: 7, 76: 7, 79: 7, 5: 6, 6: 6, 15: 6, 25: 6, 32: 6, 52: 6, 63: 6, 4: 5, 14: 5, 33: 5, 42: 5, 62: 5, 72: 5, 2: 4, 7: 4, 9: 4, 13: 4, 19: 4, 22: 4, 23: 4, 41: 4, 3: 3, 8: 3, 11: 3, 12: 3, 21: 3, 31: 3, 51: 3, 61: 3, 71: 3, 1: 2, 30: 2, 40: 2, 60: 2, 10: 1, 20: 1, 50: 1, 70: 1}
        self.assertEqual(actual, expected)


from ve.fp.complex_util.res_spat_dist import ResiduePositionDistribution
        
atg_dist = ResiduePositionDistribution({57: 7, 67: 7, 79: 7, 54: 6, 55: 6, 56: 6, 63: 6, 64: 6, 65: 6, 35: 5, 43: 5, 45: 5, 47: 5, 52: 5, 53: 5, 58: 5, 68: 5, 73: 5, 74: 5, 75: 5, 76: 5, 33: 4, 36: 4, 44: 4, 46: 4, 49: 4, 59: 4, 66: 4, 69: 4, 72: 4, 77: 4, 21: 3, 31: 3, 32: 3, 34: 3, 37: 3, 41: 3, 42: 3, 51: 3, 62: 3, 71: 3, 78: 3, 22: 2, 25: 2, 26: 2, 28: 2, 38: 2, 48: 2, 61: 2, 20: 1, 23: 1, 24: 1, 29: 1, 30: 1, 39: 1, 40: 1, 50: 1, 60: 1, 70: 1})

atb_dist = ResiduePositionDistribution({38: 12, 49: 12, 47: 11, 57: 10, 58: 10, 48: 9, 59: 9, 35: 8, 55: 8, 28: 7, 36: 7, 43: 7, 44: 7, 45: 7, 54: 7, 65: 7, 68: 7, 69: 7, 77: 7, 26: 6, 27: 6, 37: 6, 39: 6, 46: 6, 78: 6, 79: 6, 24: 5, 34: 5, 42: 5, 52: 5, 53: 5, 63: 5, 64: 5, 66: 5, 73: 5, 74: 5, 75: 5, 25: 4, 33: 4, 62: 4, 67: 4, 72: 4, 76: 4, 11: 3, 15: 3, 18: 3, 23: 3, 31: 3, 32: 3, 41: 3, 56: 3, 61: 3, 14: 2, 16: 2, 17: 2, 21: 2, 22: 2, 71: 2, 3: 1, 5: 1, 12: 1, 13: 1, 19: 1, 20: 1, 29: 1, 30: 1, 40: 1, 50: 1, 51: 1, 60: 1})

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
                                                 fps = PaddedComplexFingerPrint(self.c.ring_count))
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
                                                  fps = PaddedComplexFingerPrint(self.c.ring_count))
        _, res_dist = self.c.get_atb_res_spat_dist()
        actual = len(fps.fp_str(self.atb_dist, res_dist, number_type=int).split(","))
        expected = sum(self.atb_dist.values()) * fps.get_res_fp_length()

        self.assertEqual(actual, expected)

from ve.fp.complex_util.padding import PaddedComplexFingerPrint
class PaddedFingerPrintTest(unittest.TestCase):
    def setUp(self):
        self.c = ComplexClass()

        self.fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                  targets=[(0,self.c.atb.residues)], 
                                                  fps = PaddedComplexFingerPrint(self.c.ring_count))
    def test_to_pickle_ring_count(self):
        """ring count of to_pickle"""
        p = self.fps.to_pickable()
        self.assertEqual(p.ring_count, self.c.ring_count)
        
    def test_to_pickle_mapping_len(self):
        """mapping length of to_pickle"""
        p = self.fps.to_pickable()
        self.assertEqual(len(p.mapping), 7)

    def test_from_pickable(self):
        p = self.fps.to_pickable()
        fps = PaddedComplexFingerPrint.from_pickable(p, self.c)
        
        self.assertEqual(self.fps, fps)

if __name__ == '__main__':
    unittest.main()
        
