import os
import logging
import unittest

from ve.util.load_pdb import load_pdb_struct

from ve.machine_setting import base as proj_dir

from fake_class import TestResidue

test_data_dir = os.path.join(proj_dir, "fp/test/data")


from ve.fp.complex_util.paraepi import init_find_epiparatope_trait
from ve.fp.complex_util.split_cylinder import init_gcb_split_cylinder_trait
from ve.fp.complex_util.axial_plane import init_complex_axial_plane_trait


from common import NumericTestCase

class GCBResidue(TestResidue):
    """Residue class for geometric-based-split-cylinder method"""
    
    def __init__(self, residue):
        TestResidue.__init__(self, residue)
        from ve.fp.residue_util.res_geom import init_geom_trait
        init_geom_trait(self)


class GCBSplitCylinderTestCase(NumericTestCase):
    """Test case for the geometric-center-based cylinder split method"""
    
    def setUp(self):

        #load the antigen and antibody
        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), GCBResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), GCBResidue)

        #for load cache purpose
        self.c_id = "1SLG_D"
        
        #gcb split cylinder trait
        init_gcb_split_cylinder_trait(self)
        
    def test_general_case_atg_side(self):
        """antigen side test"""
        fp = self.get_atg_fp_by_split_cylinder()

        self.assertEqual(len(fp.fp_str().split(",")), 80)
        self.assertEqual(fp.fp_str(int),' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 0, 0, 0, 1, 5, 3, 5, 4, 1, 1, 0, 0, 3, 3, 8, 5, 5, 0, 1, 0, 0, 1, 1, 2, 4, 1, 2, 1, 1, 0, 1, 1, 2, 2, 4, 3, 2, 2, 1, 0, 0, 2, 1, 2, 2, 3, 0, 2, 0, 1')

    def test_general_case_atb_side(self):
        """antibody side test"""
        fp = self.get_atb_fp_by_split_cylinder()

        self.assertEqual(len(fp.fp_str().split(",")), 80)
        self.assertEqual(fp.fp_str(int), ' 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0')

if __name__ == "__main__":
    unittest.main()

