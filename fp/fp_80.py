import numpy as np

from ve.util.residue import BaseResidue
from ve.util.load_pdb import load_pdb_struct
from ve.util.dist import ResDistCache,FP105DistCache

from residue_util import init_resdist_util, init_neighbour_util, init_axial_plane_util, init_electric_fp_util
from complex_util import init_split_cylinder_util,init_triangle_util
from complex import TriangleComplex
from res_triangle import ResTriangle

from ve.config import *

class Residue(BaseResidue):
    def __init__(self,*args,**kw):
        BaseResidue.__init__(self,*args,**kw)
        self.center =np.average(np.array([a.xyz for a in self.atom]),0)
        
        #adding the required methods
        init_resdist_util(self)
        init_neighbour_util(self)
        
        init_axial_plane_util(self)        

        init_electric_fp_util(self)
        
    def gen_last30_fp(self, others):
        first = others[0]
        if isinstance(first, ResTriangle):
            return self.gen_last30_fp_triangle(others)
        elif isinstance(first, Residue):
            return self.gen_last30_fp_residue(others)
        else:
            raise TypeError("Only ResTriangle and Residue is ok.")
       
            
    def gen_last30_fp_triangle(self,others):
        self.gen_electric_fp_triangle(others)
        return self.fp

    def gen_last30_fp_residue(self,others):
        self.gen_electric_fp(others)
        return self.fp

class GenericComplex(TriangleComplex):
    def __init__(self,complex_id, antigen, antibody, atg_params_wrapper, atb_params_wrapper):
        TriangleComplex.__init__(self,complex_id, antigen,antibody)
        
        init_split_cylinder_util(self)

        self.atg_param_50, self.atg_param_30 = atg_params_wrapper()#use wrapper to delay expression evaluation
        self.atb_param_50, self.atb_param_30 = atb_params_wrapper()

        self.set_atg_atb_center()

        self.distcache = FP105DistCache()

    def gen_antigen_fp(self):
        fp = self.get_fp_generic(** self.atg_param_50)
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp(self.atg_param_30))
        return fp
        
    def gen_antibody_fp(self):
        fp = self.get_fp_generic(** self.atb_param_50)
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp(self.atb_param_30))
        return fp

    def set_atg_atb_center(self):
        """antigen and antibody center"""
        pts = [a.xyz for r in self.atg.residues for a in r.atom]
        self.atg_center = np.average(np.array(pts),0)

        pts = [a.xyz for r in self.atb.residues for a in r.atom]
        self.atb_center = np.average(np.array(pts),0)

        for r in self.atg.residues:
            r.set_axial_plane(self.atg_center)
        for r in self.atb.residues:
            r.set_axial_plane(self.atb_center)

class ComplexSingle(GenericComplex):
    def __init__(self,complex_id, antigen, antibody):
        GenericComplex.__init__(self,complex_id, antigen, antibody,
                                lambda :
                                    (dict(name="antigen",
                                          bases=self.atg.residues,
                                          targets=[(0,self.triangles)]),
                                     
                                     self.triangles),
                                lambda :
                                    (dict(name="antibody",
                                          bases=self.atb.residues,
                                          targets=[(0,self.atb.residues)]),
                                     self.atb.residues))

class ComplexDual(GenericComplex):
    def __init__(self,complex_id, antigen, antibody):
        GenericComplex.__init__(self,complex_id, antigen, antibody,
                                lambda :
                                    (dict(name="antigen",
                                          bases=self.atg.residues,
                                          targets=[(0,self.atb.residues)]),
                                     self.atb.residues),
                                lambda :
                                    (dict(name="antibody",
                                          bases=self.atb.residues,
                                          targets=[(0,self.triangles)]),
                                     self.triangles))

    

def single_test(complex_id):            
    data_dir = os.path.join(data237_complex_root,complex_id)

    print data237_complex_root

    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue)
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue)

    c = ComplexDual(complex_id, antigen, antibody)
    atg_fp = c.gen_antigen_fp()
    atb_fp = c.gen_antibody_fp()
    
    
    atb_fp.tofile("%s_atb.fp" %complex_id)
    atg_fp.tofile("%s_atg.fp" %complex_id)

    return atg_fp, atb_fp

if __name__ == "__main__":
    single_test("1SLG_D")
    single_test("1N4X_L")
    single_test("1JV5_A")
    single_test("1STS_B")



