from ve.fp.complex_util.paraepi import *

from common import *

class FindParatopeTest(unittest.TestCase):
    """ Test case for find paratope method"""
    
    def setUp(self):

        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), TestResidue)

        self.c_id = "1SLG_D"
        
        init_find_epiparatope_trait(self)

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
        paratope = self.find_paratope(True)
        
        self.assertEqual(self.extract_res_ids(paratope), [1 ,3 ,4 ,5 ,6 ,7])
        
    def test_epitope_general_case(self):
        epitope = self.find_epitope(True)
        self.assertGreater(len(epitope), 0)
        self.assertEqual(self.extract_res_ids(epitope), [23 ,25 ,27 ,43 ,45 ,46 ,47 ,54 ,79 ,86 ,88 ,89 ,90 ,92 ,108 ,110 ,112 ,124 ,128])

from ve.fp.complex_util.paraepi import init_io_trait

class IOTestCase(unittest.TestCase):
    def setUp(self):
        #register the trait
        init_io_trait(self)
        
        #load the antibody and antigen
        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), TestResidue)
        
        self.c_id = "1SLG_D"

        #find the paratope and epitope
        self.find_paratope()
        self.find_epitope()

        
    def test_epitope_str_no_empty(self):
        """ensure the epitope string is not empty"""
        self.assertGreater(len(self.epitope_str()), 0)
        
    def test_paratope_str_no_empty(self):
        """ensure the paraitope string is not empty"""
        self.assertGreater(len(self.paratope_str()), 0)
        
    def test_epitope_io(self):
        from ve.config import data237_paraepi_root as paraepi_path
        #write the epitope to file
        self.write_epitope()
        
        #target file path
        output_path = os.path.join(paraepi_path,self.c_id,"epitope.pdb")
        
        #if it exists
        existed = os.path.exists(output_path)
        
        self.assertTrue(existed)
        
        if existed:
            os.remove(output_path)

    def test_paratope_io(self):
        from ve.config import data237_paraepi_root as paraepi_path
        #write the paratope to file
        self.write_paratope()
        
        #target file path
        output_path = os.path.join(paraepi_path,self.c_id,"paratope.pdb")
        
        #if it exists
        existed = os.path.exists(output_path)
        
        self.assertTrue(existed)
        
        if existed:
            os.remove(output_path)

if __name__ == "__main__":
    logging.basicConfig( stream=sys.stderr )
    logging.getLogger("paraepi test").setLevel( logging.DEBUG )
    
    unittest.main()
