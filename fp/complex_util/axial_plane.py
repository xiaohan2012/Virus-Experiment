# -!- coding=utf8 -!-
import numpy as np

from types import MethodType

from geom import *

import sys
import logging

logging.basicConfig( stream=sys.stderr )
logging.getLogger("Axial Plane").setLevel( logging.DEBUG )
logger = logging.getLogger("Axial Plane")

#calculate axial plane utility

def get_epi_center(self):
    """get the epitope center"""
    if hasattr(self, "epi_center"):
        return self.epi_center

    #iterate through all in atoms in the epitope part and calculate the average sum of their x,y,z coordinates    
    self.epi_center = Point(np.average(np.array([a.xyz for r in self.find_epitope() for a in r.atom]),0))    
    
    #return it
    return self.epi_center

def get_para_center(self):
    """get the paratope center"""
    if hasattr(self, "para_center"):
        return self.para_center
    
    self.para_center = Point(np.average(np.array([a.xyz for r in self.find_paratope() for a in r.atom]),0))    
    return self.para_center

def get_paraepi_center(self):
    """
    set the geometric center of the epitope and paratope part of the complex
    as well as the center of epitope *plus* paratope
    """
    if hasattr(self, "paraepi_center"): return self.paraepi_center

    from ve.fp.geom import Point
        
    #the center of epitope and paratope altogether
    from itertools import chain
    self.paraepi_center = Point(np.average(np.array([a.xyz for r in chain(self.paratope, self.epitope) for a in r.atom]),0))
    
    return self.paraepi_center

def get_axial_plane(self):
    """get the axial plane(中轴面)of the compelx"""
    
    if hasattr(self,"axial_plane"): return self.axial_plane

    #first set the epitope and paratope center
    
    #vector between epitope center and paratope center
    norm_vec = self.get_epi_center() - self.get_para_center()
    
    #the point that the plane passes through
    point = self.get_paraepi_center()
    
    #calculate the axial plane based on the point and the vector
    from ve.fp.geom import get_perp_plane
    plane = get_perp_plane(norm_vec, point)

    #set it
    self.axial_plane = plane
    
    return self.axial_plane

def init_complex_axial_plane_trait(self):
    """
    Dependency:
    1, find paratope and epitope
    """
    
    #find_epiparatope trait
    from ve.fp.complex_util.paraepi import init_find_epiparatope_trait
    init_find_epiparatope_trait(self)

    self.get_axial_plane = MethodType(get_axial_plane, self)
    
    #getter method for three types of center
    self.get_epi_center = MethodType(get_epi_center, self)
    self.get_para_center = MethodType(get_para_center, self)
    self.get_paraepi_center = MethodType(get_paraepi_center, self)

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

def init_atg_atb_axial_plane_trait(self):
    self.set_atg_atb_center = MethodType(set_atg_atb_center, self)

