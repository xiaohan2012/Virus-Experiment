import unittest
from ve.util.logger import make_logger
from ve.util.load_pdb import load_pdb_struct, load_complexes, complex_ids

class NumericTestCase(unittest.TestCase):
    def assertArrayEqual(self, first, second):
        self.assertTrue((first == second).all())

    def assertArrayAlmostEqual(self, first, second):
        from itertools import izip
        for a,b in izip(first, second):
            self.assertAlmostEqual(a,b)

            
from fake_class import TestResidue
from ve.util.complex import BaseComplex

def make_complex_class(cls, residue_class = TestResidue):
    """
    create new complex class that extends from `cls`
    """
    class ComplexClass(BaseComplex, cls):
        def __init__(self,**kwargs):
            atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), residue_class)
            atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), residue_class)

            c_id = "1SLG_D"

            super(ComplexClass,self).__init__(complex_id = c_id, antigen = atg, antibody = atb, **kwargs)

        
    return ComplexClass

def load_base_complex():
    """
    () -> BaseComplex
    
    load testcase of basic complex type with residue type of the basic type
    """
    
    atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"))
    atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"))
    
    c_id = "1SLG_D"
    
    return BaseComplex(complex_id = c_id, antigen = atg, antibody = atb)

def make_residue_class(cls):
    class ResidueClass(TestResidue, cls):
        def __init__(self, residue, **kwargs):
            super(ResidueClass,self).__init__(residue, **kwargs)
            
    return ResidueClass
    
from ve.fp.residue_util.geom import GeometryTrait
GeometryResidue = make_residue_class(GeometryTrait)

import os

from ve.machine_setting import base as proj_dir
test_data_dir = os.path.join(proj_dir, "fp/test/data")


