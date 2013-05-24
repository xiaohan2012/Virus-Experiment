from ve.util.residue import BaseResidue
from ve.util.load_pdb import load_pdb_struct
#from ve.util.dist import FP105DistCache

from ve.fp.residue_util.neighbour import FindNeighbourTrait
from ve.fp.residue_util.axial_plane import AtgAtbBasedAxialPlaneTrait
from ve.fp.residue_util.physi_chemi import PhysiChemiTrait
from ve.fp.residue_util.geom import GeometryTrait


from complex import TriangleComplex
from res_triangle import ResTriangle

from ve.config import *

class Residue(BaseResidue, FindNeighbourTrait,#triangle requires it, a little bit akward
              PhysiChemiTrait, #the last 30 bits fp
              AtgAtbBasedAxialPlaneTrait,#the first 50 requires it
              GeometryTrait
):
    def gen_last30_fp(self, others):
        first = others[0]
        if isinstance(first, ResTriangle):
            return self.gen_electric_fp_triangle(others)
        elif isinstance(first, Residue):
            return self.gen_electric_fp(others)
        else:
            raise TypeError("Only ResTriangle and Residue is ok.")

from ve.fp.complex_util.split_cylinder import SplitCylinderTrait

class GenericComplex(TriangleComplex):
    def __init__(self,complex_id, antigen, antibody):
        super(GenericComplex,self).__init__(complex_id, antigen, antibody)
        
        #self.distcache = FP105DistCache()

    def gen_antigen_fp(self):
        """
        antigen side residue finger print
        """
        fp = self.get_fp_generic(bases=self.atg.residues,
                                 targets=[(0,self.triangles)]) #from split_cylinder trait
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp(self.get_triangles()))
        return fp
        
    def gen_antibody_fp(self):
        """
        antibody side residue finger print
        """
        fp = self.get_fp_generic(name="antibody",
                                 bases=self.atb.residues,
                                 targets=[(0,self.atb.residues)])
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp(self.atb.residues))
        return fp

from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait

class ComplexUsingComplexPlane(GenericComplex, SplitCylinderViaComplexPlaneTrait):
    pass

from ve.fp.complex_util.split_cylinder import SplitCylinderViaResiduePlaneTrait

class ComplexUsingAtgAtbPlane(GenericComplex, SplitCylinderViaResiduePlaneTrait):
    pass


def single_test(complex_id):            
    data_dir = os.path.join(data237_complex_root,complex_id)

    print data237_complex_root

    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue)
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue)

    c = ComplexDual(complex_id, antigen, antibody)
    atg_fp = c.gen_antigen_fp()
    atb_fp = c.gen_antibody_fp()
    
    type = "double"
    atb_fp.tofile("%s_atb_%s.fp" %(complex_id, type))
    atg_fp.tofile("%s_atg_%s.fp" %(complex_id, type))

    return atg_fp, atb_fp

if __name__ == "__main__":
    single_test("1SLG_D")
    single_test("1N4X_L")
    single_test("1JV5_A")
    single_test("1STS_B")



