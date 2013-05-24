import unittest

import os
import logging

from ve.util.load_pdb import load_pdb_struct

from ve.machine_setting import base as proj_dir

test_data_dir = os.path.join(proj_dir, "fp/test/data")

from fake_class import TestResidue

from common import *

from ve.fp.complex_util.axial_plane import HasAxialPlaneTrait

ComplexClass = make_complex_class(HasAxialPlaneTrait)

class ParaEpiAxialPlaneTest(NumericTestCase):
    def setUp(self):
        self.c = ComplexClass()
        

    def test_general_case(self):
        self.assertArrayAlmostEqual(self.c.get_para_center().tolist(), [ 19.0855, 1.95487879, 16.1439697 ])
        self.assertArrayAlmostEqual(self.c.get_epi_center().tolist(), [ 19.49365263,   6.60098947,  12.54066316])
        self.assertArrayAlmostEqual(self.c.get_paraepi_center().tolist(), [ 19.38842578,   5.40316406,  13.46964063])
    


if __name__ == "__main__":
    unittest.main()
