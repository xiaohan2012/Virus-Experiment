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
    #complex_ids()
    #ids = ["1FJ1_E", "1FJ1_F", "0HEZ_E", "1GC1_G", "3B2U_I", "3L5W_I", "1G9M_G"]
    ids = ["1SLG_D"]
    cs = load_complexes(ids, complex_cls = complex_cls, residue_cls = residue_cls)
    for c in cs:
        c.fp_to_cache()
        """
        try:
            c.fp_to_cache()
        except:
            print "error! %s" %c.c_id
        """
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
                                 
