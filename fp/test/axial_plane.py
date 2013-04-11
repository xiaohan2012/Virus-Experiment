import unittest

import os
import logging

from ve.fp.complex_util.axial_plane import init_complex_axial_plane_trait
from ve.fp.complex_util.paraepi import init_find_epiparatope_trait

from ve.util.load_pdb import load_pdb_struct

from ve.machine_setting import base as proj_dir

test_data_dir = os.path.join(proj_dir, "fp/test/data")

from fake_class import TestResidue

from common import NumericTestCase

class ParaEpiAxialPlaneTest(NumericTestCase):
    def setUp(self):
        init_complex_axial_plane_trait(self)

        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), TestResidue)
        self.c_id = "1SLG_D"
        

    def test_general_case(self):
        self.assertArrayAlmostEqual(self.get_para_center().tolist(), [ 19.0855, 1.95487879, 16.1439697 ])
        self.assertArrayAlmostEqual(self.get_epi_center().tolist(), [ 19.49365263,   6.60098947,  12.54066316])
        self.assertArrayAlmostEqual(self.get_paraepi_center().tolist(), [ 19.38842578,   5.40316406,  13.46964063])
    


if __name__ == "__main__":
    unittest.main()
