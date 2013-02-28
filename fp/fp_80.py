import numpy as np

from ve.util.residue import BaseResidue
from ve.util.load_pdb import load_pdb_struct
from ve.util.dist import ResDistCache,FP105DistCache

from residue_util import init_resdist_util, init_neighbour_util, init_axial_plane_util, init_electric_fp_util
from complex_util import init_split_cylinder_util,init_triangle_util
from complex import TriangleComplex

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
        

    def gen_last30_fp_triangle(self,others):
        self.gen_electric_fp_triangle(others)
        return self.fp

    def gen_last30_fp_residue(self,others):
        self.gen_electric_fp(others)
        return self.fp

class Complex(TriangleComplex):
    def __init__(self,complex_id, antigen, antibody):
        TriangleComplex.__init__(self,complex_id, antigen,antibody)
        
        init_split_cylinder_util(self)

        self.atg_param = dict(name="antigen",
                         bases=self.atg.residues,
                         targets=[(0,self.triangles)])
        self.atb_param = dict(name="antibody",
                         bases=self.atb.residues,
                         targets=[(0,self.atb.residues)])
        self.set_atg_atb_center()

        self.distcache = FP105DistCache()

    def gen_antigen_fp(self):
        fp = self.get_fp_generic(** self.atg_param)
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp_triangle(self.triangles))
        return fp
        
    def gen_antibody_fp(self):
        fp = self.get_fp_generic(** self.atb_param)
        for r in fp.residues():
            fp[r].append(r.gen_last30_fp_residue(self.atb.residues))
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


def single_test(complex_id):            
    data_dir = os.path.join(data237_complex_root,complex_id)

    print data237_complex_root

    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue)
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue)

    c = Complex(antigen, antibody)
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



