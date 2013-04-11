from ve.util.load_pdb import load_pdb_struct


from common import *
class Residue(TestResidue):
    def __init__(self, residue):
        TestResidue.__init__(self, residue)
        
        #center trait
        from ve.fp.residue_util.res_geom import init_geom_trait
        init_geom_trait(self)
    
class PhysChemiTestCase(unittest.TestCase):
    def setUp(self):
        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)

        self.c_id = "1SLG_D"
        
        from ve.fp.complex_util.physi_chemi import init_physi_chemi_fp_trait
        init_physi_chemi_fp_trait(self)
        

    def test_general_case_atg(self):
        fp =  self.get_physi_chemi_atg_fp()
        self.assertEqual(len(fp.fp_str().split(",")), 30)
        self.assertEqual(fp.fp_str(),"0.00,0.00,2.65,7.48,10.51,8.57,10.62,10.00,7.87,3.20,0.00,0.00,0.00,0.18,0.48,0.45,0.68,0.32,0.30,0.04,0.00,0.00,1.00,4.00,8.00,9.00,23.00,7.00,11.00,3.00")

    def test_general_case_atb(self):
        fp =  self.get_physi_chemi_atb_fp()
        self.assertEqual(len(fp.fp_str().split(",")), 30)
        self.assertEqual(fp.fp_str(), "0.00,0.00,2.62,0.00,0.05,2.07,0.00,0.00,0.00,0.00,0.05,0.00,0.14,0.00,0.04,0.14,0.00,0.00,0.00,0.00,2.00,0.00,3.00,0.00,1.00,1.00,0.00,0.00,0.00,0.00")

if __name__ == "__main__":
    unittest.main()
