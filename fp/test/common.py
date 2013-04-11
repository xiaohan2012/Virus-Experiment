import unittest

class NumericTestCase(unittest.TestCase):
    def assertArrayEqual(self, first, second):
        self.assertTrue((first == second).all())

    def assertArrayAlmostEqual(self, first, second):
        from itertools import izip
        for a,b in izip(first, second):
            self.assertAlmostEqual(a,b)
        


import os
import logging
import unittest

from ve.machine_setting import base as proj_dir

from fake_class import TestResidue

test_data_dir = os.path.join(proj_dir, "fp/test/data")

from ve.util.load_pdb import load_pdb_struct
