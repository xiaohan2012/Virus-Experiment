import math
import os

import numpy as np

from types import MethodType

from fp import FP0112, BaseComplexFingerprint
from geom import *
from res_triangle import ResTriangle

def split_cylinder_method(self, name="",bases=[],targets=[]):
    """
    sample input data
    "name":"epitope triangles",
    "bases":triangles,
    "targets":[(0,triangles),(50,atb.residues)]
    """
    res_fp_length=50
    def get_idx(self,dist2plane,dist2center):
        min_height,max_height = self.cylinder_height_range
        
        #in range
        if (dist2plane >= min_height and dist2plane <= max_height) and dist2center < self.cylinder_radius:
            layer_idx = self.layer_count/2 + int(math.floor(dist2plane / self.cylinder_height_step))
            ring_idx = int(math.floor(dist2center / self.cylinder_radius_step))
            
            return layer_idx * self.layer_size + ring_idx

    self.get_idx = MethodType(get_idx, self)

#    fps = FP0112(targets)
    fps = BaseComplexFingerprint()

    print "processing %s"  %(name)
    for base in bases:
        if hasattr(self,"axial_plane"): #has axial plane
            axial_plane = self.axial_plane            
        else:
            axial_plane = base.axial_plane#use residue axial_plane
        base_center = axial_plane.get_perp_point(base.center)
        for index_base,iterables in targets:
            for other in iterables:
                perp_point = axial_plane.get_perp_point(other.center)
                
                dist2plane = axial_plane.dist2point(other.center)
                dist2center = perp_point.dist2point(base_center)
                
                    #same as the base, ignore
                if base == other:
                    continue
                
                idx = self.get_idx(dist2plane, dist2center)
                if idx:
                    if not fps.has_res(base):#not registered yet
                        fps.add_res(base,res_fp_length)
                    fps[base][index_base + idx] += 1
    return fps

def init_split_cylinder_util(self,cylinder_radius = 10, 
                             cylinder_radius_step = 2,
                             cylinder_height = 20,
                             cylinder_height_step = 2):

        self.cylinder_radius  = cylinder_radius
        self.cylinder_radius_step = cylinder_radius_step
        self.cylinder_height = cylinder_height
        self.cylinder_height_step = cylinder_height_step

        self.cylinder_height_range = (self.cylinder_height/2 - self.cylinder_height , self.cylinder_height/2)

        self.layer_size = self.cylinder_radius / self.cylinder_radius_step
        self.layer_count = self.cylinder_height / self.cylinder_height_step
                
        self.get_fp_generic = MethodType(split_cylinder_method, self)



from pickle import load, dump
def _find_paratope(self):
    """antibody side"""
    attr_name = "paratope"

    if self.cache_exists(attr_name):
        setattr(self, attr_name, self.load_cache(attr_name))
        return 

    for b_r in self.atb.residues:
        for g_r in self.atg.residues:
            dist =b_r.dist2residue(g_r)
            if dist <=self.paratope_threshold:
                self.paratope.append(b_r)
                break
    self.save_cache(attr_name)

def _find_epitope(self):
    """antigen side"""
    attr_name = "epitope"

    if self.cache_exists(attr_name):
        setattr(self, attr_name, self.load_cache(attr_name))
        return 

    for g_r in self.atg.residues:
        for b_r in self.atb.residues:
            if g_r.dist2residue(b_r) <= self.paratope_threshold:
                self.epitope.append(g_r)
                break
    self.save_cache(attr_name)


def cache_exists(self,cache_type):
    par_dir = ".cache/%s" %cache_type
    if not os.path.exists(par_dir):
        os.makedirs(par_dir)

    cache_path = os.path.join(".cache/%s" %cache_type,
                              self.id)
    mark = os.path.exists(cache_path)
    if mark:
        print "cache exists", self.id,cache_type
    return mark

def paratope_cl_process(self,content):
    return filter(lambda r: r.resnum in content, 
                  self.atb.residues)

def epitope_cl_process(self,content):
    return filter(lambda r: r.resnum in content, 
                  self.atg.residues)

def epi_para_cd_process(self,obj):
    return map(lambda r: r.resnum, obj)


def load_cache(self,cache_type):
    cache_path = os.path.join(".cache/%s" %cache_type,
                              self.id)
    cache_content = load(open(cache_path,"r"))
    return self.cache_load_methods[cache_type](cache_content)

def save_cache(self,cache_type):
    cache_path = os.path.join(".cache/%s" %cache_type,
                              self.id)
    val = self.cache_dump_methods[cache_type](getattr(self,cache_type))
    print val
    dump(val, open(cache_path,"w"))
    
def init_find_epiparatope_util(self):
    self.paratope_threshold = 5
    self.epitope_threshold = 5

    self._find_paratope = MethodType(_find_paratope, self)
    self._find_epitope = MethodType(_find_epitope, self)

    self.cache_exists = MethodType(cache_exists, self)
    self.load_cache = MethodType(load_cache, self)
    self.save_cache = MethodType(save_cache, self)

    self.paratope_cl_process = MethodType(paratope_cl_process, self)
    self.epitope_cl_process = MethodType(epitope_cl_process, self)

    self.epi_para_cd_process = MethodType(epi_para_cd_process, self)

    
    self.cache_load_methods = {
        "paratope": self.paratope_cl_process,
        "epitope": self.epitope_cl_process,
        }
    self.cache_dump_methods = {
        "paratope": self.epi_para_cd_process,
        "epitope": self.epi_para_cd_process,
        }

def _find_triangles(self):
    for r in self.epitope:
        r.find_neighbours(self.epitope,self.neighbour_threshold)
        
    triangles = []
    byproduct = []#stores residue tuples that are possible to form triangle with external residue

    for parent in self.epitope:
        for child in parent.get_neighbours():
            temp_list = []
            for grandchild in child.get_neighbours():
                # and child != grandchild and parent != child
                if parent != grandchild:
                    tri = ResTriangle([parent,child,grandchild])
                    temp_list.append(tri)
            if not temp_list:#it is empty
                byproduct.append((parent,child))
            triangles += temp_list
    remaining_residues = set(self.atg.residues)-set(self.epitope)
    for r1,r2 in byproduct:
        for r3 in remaining_residues:
            if r1.dist2residue(r3) <= self.neighbour_threshold and\
                    r2.dist2residue(r3) <= self.neighbour_threshold:
                tri = ResTriangle([r1,r2,r3])
                triangles.append(tri)
    print "%d more triangles" %len(byproduct)
    triangles = list(set(triangles))
    print "found %d triangles in total" %len(triangles)
    
    for t in triangles:
        t.cal_center()
        
    self.triangles = triangles

def init_triangle_util(self):
    self.neighbour_threshold = 4
    self._find_triangles = MethodType(_find_triangles, self)


#calculate axial plane utility
def _set_pe_center(self):
    """epitope and paratope center"""
    pts = []
    for r in self.epitope:
        for a in r.atom:
            pts.append(a.xyz)

    self.epi_center = np.average(np.array(pts),0)

    for r in self.paratope:
        for a in r.atom:
            pts.append(a.xyz)

    self.pe_center = np.average(np.array(pts),0)

    self.epi_center, self.pe_center
    
def _set_axial_plane(self):
    norm_vec = self.epi_center - self.pe_center
    point = self.pe_center
    plane = get_perp_plane(norm_vec, point)
    self.axial_plane = plane

def set_axial_plane(self):
    self._set_pe_center()
    self._set_axial_plane()

def init_complex_axial_plane_util(self):
    self._set_axial_plane = MethodType(_set_axial_plane, self)
    self._set_pe_center = MethodType(_set_pe_center, self)

    self.set_axial_plane = MethodType(set_axial_plane, self)


#set antigen and antibody centers utility
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

def init_atg_atb_axial_plane_util(self):
    self.set_atg_atb_center = MethodType(set_atg_atb_center, self)
