# -*- coding=utf8 -*-
"""
100,100以及150位指纹

对于抗原侧三角形以及抗体侧氨基酸两种情况的100,100,150指纹

共6种指纹

它们如何计算的见该链接:  https://docs.google.com/spreadsheet/ccc?key=0Avn4Kj4hHHDNdGdHdVYyR3N6Z1U5UXBiX2tFczNNY1E#gid=0
"""
import os
import sys
import math
import numpy as np

from itertools import chain
from collections import defaultdict
from glob import glob
from customcollections import OrderedDefaultDict

from types import MethodType

from complex import TriangleComplex

from geom import *

from residue_util import init_resdist_util,init_neighbour_util
from complex_util import init_split_cylinder_util

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
        
        init_split_cylinder_util(self)

        print "setting axial plane"
        self._set_pe_center()
        self._set_axial_plane()
        
        
        params =[
            dict(name="antigen,triangle",
                 bases=self.triangles,
                 targets=[(0,self.triangles),(50,self.atb.residues)]),
            dict( name="antigen,residue",
                  bases=self.triangles,
                  targets = [(0,self.atg.residues),(50,self.atb.residues)]),
            dict( name="antigen,combined",
                  bases=self.triangles,
                  targets = [(0,self.triangles),(50,self.atg.residues),(100,self.atb.residues)]),
            dict(name="antibody,triangle",
                 bases=self.atb.residues,
                 targets=[(0,self.triangles),(50,self.atb.residues)]),
            dict( name="antibody,residue",
                  bases=self.atb.residues,
                  targets = [(0,self.atg.residues),(50,self.atb.residues)]),
            dict( name="antibody,combined",
                  bases=self.atb.residues,
                  targets = [(0,self.triangles),(50,self.atg.residues),(100,self.atb.residues)])
            ]

        #delay and force, lazy computing
        self.delay = map(lambda (x): (False,x), params)
        
        def force(n):
            try:
                mark,other = self.delay[n-1]
            except IndexError:
                print "unknown fp id"
                return None
            if mark:#in cache
                print "hit cache",n
                return other
            else:
                fp = self.get_fp_generic(**other)
                self.delay[n-1] = (True,fp)
                return fp
                
        self.nth_fp = force
    
    def gen_fps(self):
        return [self.nth_fp(i) for i in xrange(1,7)]

    def _set_pe_center(self):
        """epitope and paratope center"""
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
    c.nth_fp(1)
    fps = c.gen_fps()
    print fps[2].fp_str()
    return fps

if __name__ == "__main__":
    #data237_fp_gen()
    fps = single_test("1SLG_D")
1
