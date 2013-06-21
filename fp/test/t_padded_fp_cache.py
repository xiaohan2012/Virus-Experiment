import os
from common import *
from ve.fp.complex_util.padding import PaddedComplexFingerPrint
from ve.fp.fp_75_gen import get_15bits

class ComplexPlanePaddedComplexFingerPrintCacheTest(unittest.TestCase):
    """ Test case for padded complex fp cache"""
    
    def setUp(self):
        #cache init
        from ve.fp.complex_util.cache import ComplexFingerPrintCache as C
        self.C = C

        from ve.fp.fp_175_padded import ComplexPlaneBasedComplex
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        c_id = "1SLG_D" 
        self.c = ComplexPlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)
        
    def test_load_atg_tri(self):
        """
        test if cached fingerprint can be loaded as padded fp
        in this case, the atg side iterating through triangles based on complex plane
        and whether the loaded one equals to the computed one
        """
        self.C.set_signature("atg_tri_complex_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                  targets=[(0,self.c.get_triangles())],
                                                  fps = PaddedComplexFingerPrint())

                                                  
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)

    def test_load_atg_res(self):
        """
        the atg side iterating through residues based on complex plane
        """
        self.C.set_signature("atg_res_complex_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                  targets=[(0,self.c.atg.residues)],
                                                  fps = PaddedComplexFingerPrint())
                                                  
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)

    def test_load_atb_res(self):
        """
        the atb side iterating through residues based on complex plane
        """
        self.C.set_signature("atb_res_complex_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                          targets=[(0,self.c.atb.residues)],
                                                  fps = PaddedComplexFingerPrint())
        
                                                  
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)

class InteractiveForceTestCase(unittest.TestCase):

    def setUp(self):
        from ve.fp.complex_util.cache import ComplexFingerPrintCache as C
        self.C = C

        from ve.fp.complex_util.interactive_force_fp import InteractiveForceTrait
        Complex = make_complex_class(InteractiveForceTrait)
        
        self.c = Complex()
    
    def test_load_iforce_atg_as_rec(self):
        """test interactive force whose receptor is set to antigen"""
        expected_fp =  self.c.gen_if_residue_fp_atg(fp = PaddedComplexFingerPrint())
        
        self.C.set_signature("interactive_force_atG_as_rec")
        
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        
        print expected_fp.get_res_fp_by_resid("D-109 ")
        print actual_fp.get_res_fp_by_resid("D-109 ")

        self.assertEqual(expected_fp, actual_fp)
        
    def test_load_iforce_atb_as_rec(self):
        """test interactive force whose receptor is set to antibody"""
        expected_fp =  self.c.gen_if_residue_fp_atb(fp = PaddedComplexFingerPrint())
        
        self.C.set_signature("interactive_force_atB_as_rec")
        
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)
        

class ResiduePlanePaddedComplexFingerPrintCacheTest(unittest.TestCase):
    """ Test case for padded complex fp cache"""
    
    def setUp(self):
        #cache init
        from ve.fp.complex_util.cache import ComplexFingerPrintCache as C
        self.C = C

        from ve.fp.fp_175_padded import ResiduePlaneBasedComplex
        from ve.fp.fp_80 import Residue
        atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), Residue)
        atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), Residue)
        
        c_id = "1SLG_D" 
        self.c = ResiduePlaneBasedComplex(complex_id = c_id, antigen = atg, antibody = atb)
        
    def test_load_atg_tri(self):
        """
        test if cached fingerprint can be loaded as padded fp
        in this case, the atg side iterating through triangles based on complex plane
        and whether the loaded one equals to the computed one
        """
        self.C.set_signature("atg_tri_residue_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                  targets=[(0,self.c.get_triangles())],
                                                  fps = PaddedComplexFingerPrint())

                                                  
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)

    def test_load_atg_res(self):
        """
        the atg side iterating through residues based on complex plane
        """
        self.C.set_signature("atg_res_residue_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atg.residues,
                                                          targets=[(0,self.c.atg.residues)],
                                                          fps = PaddedComplexFingerPrint())
                                                  
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)

    def test_load_atb_res(self):
        """
        the atb side iterating through residues based on complex plane
        """
        self.C.set_signature("atb_res_residue_plane")

        expected_fp = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                                          targets=[(0,self.c.atb.residues)],
                                                          fps = PaddedComplexFingerPrint())
        
        actual_fp = self.C.load(self.c.c_id, self.c, complex_fp_cls =  PaddedComplexFingerPrint)
        self.assertEqual(expected_fp, actual_fp)
            
if __name__ == '__main__':
    unittest.main()