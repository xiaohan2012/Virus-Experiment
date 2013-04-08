import unittest
import os
import logging

from ve.fp.complex_util.paraepi import *

from ve.util.load_pdb import load_pdb_struct

from ve.machine_setting import base as proj_dir

from fake_class import TestResidue

test_data_dir = os.path.join(proj_dir, "fp/test/data")

class FindParatopeTest(unittest.TestCase):
    """ Test case for find paratope method"""
    
    def setUp(self):
        
        init_find_epiparatope_trait(self)

        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), TestResidue)

        self.c_id = "1SLG_D"
        

    def extract_res_ids(self,reses):
        """
        Params:
        reses: list of Residue
        
        Returns:
        list of int, the residue ids

        helper function, extract list of residue their ids
        """
        return map(lambda r: r.resnum, reses)

    def test_paratope_general_case(self):
        self.find_paratope(True)
        self.assertEqual(self.extract_res_ids(self.paratope), [1 ,3 ,4 ,5 ,6 ,7])
        
    def test_epitope_general_case(self):
        self.find_epitope(True)
        self.assertEqual(self.extract_res_ids(self.epitope), [23 ,25 ,27 ,43 ,45 ,46 ,47 ,54 ,79 ,86 ,88 ,89 ,90 ,92 ,108 ,110 ,112 ,124 ,128])
        
if __name__ == "__main__":
    logging.basicConfig( stream=sys.stderr )
    logging.getLogger("paraepi test").setLevel( logging.DEBUG )
    
    unittest.main()
