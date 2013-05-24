import os
import logging
import unittest

from ve.util.load_pdb import load_pdb_struct

from ve.machine_setting import base as proj_dir

test_data_dir = os.path.join(proj_dir, "fp/test/data")

from common import *

from fake_class import TestResidue

from ve.fp.residue_util.geom import GeometryTrait

#Residue class for geometric-based-split-cylinder method
GCBResidue = make_residue_class(GeometryTrait)
        
from ve.fp.complex_util.split_cylinder import GBCSplitCylinderTrait

ComplexClass = make_complex_class(GBCSplitCylinderTrait, residue_class = GCBResidue)

class GCBSplitCylinderTestCase(NumericTestCase):
    """Test case for the geometric-center-based cylinder split method"""
    
    def setUp(self):
        self.c = ComplexClass()
        
    def test_general_case_atg_side(self):
        """antigen side test"""
        fp = self.c.get_atg_fp_by_split_cylinder()

        self.assertEqual(len(fp.fp_str().split(",")), 80)
        self.assertEqual(fp.fp_str(int),' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 0, 1, 5, 3, 5, 4, 1, 1, 0, 0, 3, 3, 8, 5, 5, 0, 1, 0, 0, 1, 1, 2, 4, 1, 2, 1, 1, 0, 1, 1, 2, 2, 4, 3, 2, 2, 1, 0, 0, 2, 1, 2, 2, 3, 0, 2, 0, 1')

    def test_general_case_atb_side(self):
        """antibody side test"""
        fp = self.c.get_atb_fp_by_split_cylinder()

        self.assertEqual(len(fp.fp_str().split(",")), 80)
        self.assertEqual(fp.fp_str(int), ' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0')


class ComplexPlaneSplitCylinderTest(NumericTestCase):
    """Test case for the cylinder split method based on complex axial plane"""
    
    def setUp(self):
        from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait
        from ve.fp.complex_util.triangle import ResidueTriangleTrait

        from ve.fp.fp_80 import Residue as ResidueClass
        
        TempClass = make_complex_class(SplitCylinderViaComplexPlaneTrait, residue_class = ResidueClass)
        
        class ComplexClass(TempClass, ResidueTriangleTrait):
            pass
            
        self.c = ComplexClass()
        
    def test_general_case_atg_side(self):
        """antigen side test"""
        fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                 targets=[(0,self.c.get_triangles())])
        
        actual = fp.fp_str()[:404]
        print fp.fp_str()
        expected = "13,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,2.00,1.00,0.00,2.00,0.00,0.00,1.00,5.00,0.00,5.00,3.00,1.00,2.00,2.00,1.00,13.00,11.00,6.00,0.00,1.00,0.00,1.00,1.00,2.00,4.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00"
        self.assertEqual(actual, expected)

    def test_general_case_atb_side(self):
        """antibody side test"""
        fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                 targets=[(0,self.c.atb.residues)])
        
        actual = fp.fp_str()[:401]
        expected = "1,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.00,0.00,0.00,0.00,0.00,1.00,0.00,0.00,0.00,0.00,0.00,2.00,1.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00"
        self.assertEqual(actual, expected)

        
class ResiduePlaneSplitCylinderTest(NumericTestCase):
    """Test case for the cylinder split method based on residue axial plane"""
    
    def setUp(self):
        from ve.fp.complex_util.split_cylinder import SplitCylinderViaResiduePlaneTrait
        from ve.fp.complex_util.triangle import ResidueTriangleTrait

        from ve.fp.fp_80 import Residue as ResidueClass
        
        TempClass = make_complex_class(SplitCylinderViaResiduePlaneTrait, residue_class = ResidueClass)
        
        class ComplexClass(TempClass, ResidueTriangleTrait):
            pass
            
        self.c = ComplexClass()
        
    def test_general_case_atg_side(self):
        """antigen side test"""
        fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                 targets=[(0,self.c.get_triangles())])
        
        actual = fp.fp_str()[:402]
        expected = "13,0.00,0.00,2.00,2.00,4.00,2.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00"
        self.assertEqual(actual, expected)

    def test_general_case_atb_side(self):
        """antibody side test"""
        fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                 targets=[(0,self.c.atb.residues)])
        
        actual = fp.fp_str()[:401]
        expected = "1,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.00,2.00,0.00,0.00,1.00,0.00,0.00,0.00,0.00,0.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00"
        self.assertEqual(actual, expected)

        
if __name__ == "__main__":
    unittest.main()

