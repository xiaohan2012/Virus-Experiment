"""
16 sets of complex finger print based on cylinder splitting and interacting force
"""
from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait, SplitCylinderViaResiduePlaneTrait
from ve.fp.complex_util.interactive_force_fp import InteractiveForceTrait
from ve.fp.complex import TriangleComplex

from ve.fp.complex_util.cache import ComplexFingerPrintCache as C


class ComplexPlaneBasedComplex(TriangleComplex, SplitCylinderViaComplexPlaneTrait):
    """
    complex-plane based fp generation complex

    1, atg side:
    a, triangle based
    b, residue based

    2, atb side
    a, residue based
    
    """        
    def gen_atg_tri_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antigen side triangle based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atg.residues, targets = [(0, self.get_triangles())])

    def gen_atg_res_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antigen side residue based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atg.residues, targets = [(0, self.atg.residues)])

    def gen_atb_res_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antibody side residue based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atb.residues, targets = [(0, self.atb.residues)])

    def fp_to_cache(self):
        """(Complex) -> NoneType"""
        fp1, fp2, fp3 = self.gen_atg_tri_fp(), self.gen_atg_res_fp(), self.gen_atb_res_fp()

        print fp1.get_bitlength()
        print fp2.get_bitlength()
        print fp3.get_bitlength()
        
        C.set_signature("atg_tri_complex_plane")
        C.dump(self.c_id, fp1)

        C.set_signature("atg_res_complex_plane")
        C.dump(self.c_id, fp2)

        C.set_signature("atb_res_complex_plane")        
        C.dump(self.c_id, fp3)

class ResiduePlaneBasedComplex(TriangleComplex, SplitCylinderViaResiduePlaneTrait):
    """
    residue-plane based fp generation complex

    1, atg side:
    a, triangle based
    b, residue based

    2, atb side
    a, residue based
    
    """
    def gen_atg_tri_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antigen side triangle based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atg.residues, targets = [(0, self.get_triangles())])

    def gen_atg_res_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antigen side residue based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atg.residues, targets = [(0, self.atg.residues)])

    def gen_atb_res_fp(self):
        """
        (Complex) -> BaseComplexFingerprint
        
        antibody side residue based
        """
        return self.gen_fp_by_splitting_cylinder(bases = self.atb.residues, targets = [(0, self.atb.residues)])

    def fp_to_cache(self):
        """(Complex) -> NoneType"""
        fp1, fp2, fp3 = self.gen_atg_tri_fp(), self.gen_atg_res_fp(), self.gen_atb_res_fp()
        
        C.set_signature("atg_tri_residue_plane")
        C.dump(self.c_id, fp1)

        C.set_signature("atg_res_residue_plane")
        C.dump(self.c_id, fp2)

        C.set_signature("atb_res_residue_plane")        
        C.dump(self.c_id, fp3)

from ve.util.complex import BaseComplex

class InteractiveForceComplex(BaseComplex, InteractiveForceTrait):
    def fp_to_cache(self):
        """(Complex) -> NoneType"""
        fp1 = self.gen_if_residue_fp_atg()
        
        fp2 = self.gen_if_residue_fp_atb()

        C.set_signature("interactive_force_atG_as_rec")
        C.dump(self.c_id, fp1)
        C.set_signature("interactive_force_atB_as_rec")
        C.dump(self.c_id, fp2)

from ve.fp.fp_80 import Residue

def main(complex_cls, residue_cls = Residue):
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.config import data480_complex_root
    from ve.fp.complex_util.paraepi  import ParatopeNotFoundError, EpitopeNotFoundError
    
    #ids = complex_ids(data480_complex_root)
    ids = ["3B2U_E","3DVN_XY","2ARJ_R","3RKD_A","3RKD_B","3L5W_J","1BZQ_D","3NFP_I","3DVN_UV","3HI1_J","3B9K_AB","3IU3_K","3HI1_G","3B2U_A","3BN9_A","1BZQ_A","4ETQ_X","3U4E_J","3B2U_B","2ARJ_Q","1BZQ_B","3NFP_K","3IU3_I"]
    
    cs = load_complexes(ids, directory = data480_complex_root, complex_cls = complex_cls, residue_cls = residue_cls)
    
    for c in cs:
        #c.fp_to_cache()
        try:
            c.fp_to_cache()
        except EpitopeNotFoundError:
            print "EpitopeNotFoundError %s" %c.c_id
        except ParatopeNotFoundError:
            print "ParatopeNotFoundError %s" %c.c_id
        except:
            print "other error %s" %c.c_id

        
def usage():
    print """
possible args: res_plane | comp_plane | if 
    """
            
if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        usage()
    elif sys.argv[1] == "res_plane":
        main(ResiduePlaneBasedComplex)
    elif sys.argv[1] == "comp_plane":
        main(ComplexPlaneBasedComplex)
    elif sys.argv[1] == "if":
        main(InteractiveForceComplex)
    else:
        usage()
                                 
