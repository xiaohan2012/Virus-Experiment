import os
import sys
import math
import numpy as np

from itertools import chain
from collections import defaultdict
from glob import glob
from customcollections import OrderedDefaultDict

from complex import TriangleComplex

from geom import *
from fp import FP0112
from residue_util import init_resdist_util,init_neighbour_util

from ve.util.load_pdb import load_pdb_struct
from ve.util.file import fpstr2file
from ve.util.residue import BaseResidue
from ve.util.complex import BaseComplex


from ve.config import *


class Residue(BaseResidue):
    def __init__(self,*args,**kw):
        BaseResidue.__init__(self,*args,**kw)
        self.center =np.average(np.array([a.xyz for a in self.atom]),0)
        
        #adding the required methods
        init_resdist_util(self)
        init_neighbour_util(self)

    

class Complex(TriangleComplex):
    def __init__(self,antigen, antibody):
        TriangleComplex.__init__(self,antigen,antibody)
        
        self.cylinder_radius = 10
        self.cylinder_radius_step = 2
        self.cylinder_height = 20
        self.cylinder_height_step = 2

        self.cylinder_height_range = (self.cylinder_height/2 - self.cylinder_height , self.cylinder_height/2)
        print "height range", self.cylinder_height_range

        self.layer_size = self.cylinder_radius / self.cylinder_radius_step
        self.layer_count = self.cylinder_height / self.cylinder_height_step


        print "setting axial plane"
        self._set_pe_center()
        self._set_axial_plane()
        


    def get_idx(self,dist2plane,dist2center):
        min_height,max_height = self.cylinder_height_range
        
        #in range
        if (dist2plane >= min_height and dist2plane <= max_height) and dist2center < self.cylinder_radius:
            layer_idx = self.layer_count/2 + int(math.floor(dist2plane / self.cylinder_height_step))
            ring_idx = int(math.floor(dist2center / self.cylinder_radius_step))

            return layer_idx * self.layer_size + ring_idx
        #return None            

    def get_fp_generic(self,name="",bases=[],targets=[]):
        """
        sample input data
            "name":"epitope triangles",
            "bases":triangles,
            "targets":[(0,triangles),(50,atb.residues)]
       """
        fps = FP0112(targets)

        print "processing %s"  %(name)
        for base in bases:
            base_center = self.axial_plane.get_perp_point(base.center)
            for index_base,iterables in targets:
                for other in iterables:
                    perp_point = self.axial_plane.get_perp_point(other.center)

                    dist2plane = self.axial_plane.dist2point(other.center)
                    dist2center = perp_point.dist2point(base_center)
                    
                    #same as the base, ignore
                    if base == other:
                        continue
                    
                    idx = self.get_idx(dist2plane, dist2center)
                    if idx:
                        fps[base][index_base + idx] += 1
        return fps
                       
    def gen_fps(self):
        fp1 = self.get_fp_generic(name="antigen,triangle",
                              bases=self.triangles,
                              targets=[(0,self.triangles),(50,self.atb.residues)])
        fp2 = self.get_fp_generic( name="antigen,residue",
                              bases=self.triangles,
                              targets = [(0,self.triangles),(50,self.atb.residues)])        
        fp3 = self.get_fp_generic( name="antigen,combined",
                              bases=self.triangles,
                              targets = [(0,self.triangles),(50,self.atg.residues),(100,self.atb.residues)])        
        fp4 = self.get_fp_generic(name="antibody,triangle",
                              bases=self.atb.residues,
                              targets=[(0,self.triangles),(50,self.atb.residues)])
        fp5 = self.get_fp_generic( name="antibody,residue",
                              bases=self.atb.residues,
                              targets = [(0,self.triangles),(50,self.atb.residues)])        
        fp6 = self.get_fp_generic( name="antibody,combined",
                              bases=self.atb.residues,
                              targets = [(0,self.triangles),(50,self.atg.residues),(100,self.atb.residues)])        
        return [fp1,fp2,fp3,fp4,fp5,fp6]

    def _set_pe_center(self):
        pts = []
        for r in chain(self.epitope):
            for a in r.atom:
                pts.append(a.xyz)

        self.epi_center = np.average(np.array(pts),0)

        for r in chain(self.paratope):
            for a in r.atom:
                pts.append(a.xyz)

        self.pe_center = np.average(np.array(pts),0)

        self.epi_center, self.pe_center
    
    def _set_axial_plane(self):
        norm_vec = self.epi_center - self.pe_center
        point = self.pe_center
        plane = get_perp_plane(norm_vec, point)
        self.axial_plane = plane

### this time, we mean it ###
def data237_fp_gen(refresh=False):
    for path in glob(os.path.join(data237_complex_root ,"*")):
        complex_id = os.path.basename(path)

        output_path = os.path.join(data237_fp6_root,complex_id)
        if not refresh and os.path.exists(output_path):
            print("%s already processed" %complex_id)
            continue

        print "start processing %s" %complex_id

        data_dir = os.path.join(data237_complex_root,complex_id)
        antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue)
        antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue)
        try:
            c = Complex(antigen,antibody)
            fps = c.gen_fps()
        except:
            sys.stderr.write("complex %s encountered error.\n" %complex_id)
        
        if not os.path.exists(output_path):
            os.makedirs(output_path)
        
        for i,fp in enumerate(fps):
            path = os.path.join(output_path,"fp%d.txt" %(i+1))
            print "to path", path
            fp.tofile(path)

def single_test(complex_id):            
    data_dir = os.path.join(data237_complex_root,complex_id)
    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue)
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue)
    c = Complex(antigen,antibody)
    fps = c.gen_fps()
    print fps[2].fp_str()

if __name__ == "__main__":
    #data237_fp_gen()
    single_test("1SLG_D")
