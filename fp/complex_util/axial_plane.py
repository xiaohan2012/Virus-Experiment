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
def set_paraepi_center(self):
    """
    set the geometric center of the epitope and paratope part of the complex
    as well as the center of epitope *plus* paratope
    """

    #iterate through all in atoms in the epitope part and calculate the average sum of their x,y,z coordinates
    self.epi_center = np.average(np.array([a.xyz for r in self.epitope for a in r.atom]),0)

    #similar, but iterate through the paratope part
    self.pe_center = np.average(np.array([a.xyz for r in self.paratope for a in r.atom]),0)
    
    #the center of epitope and paratope altogether
    from itertools import chain
    self.paraepi_center = np.average(np.array([a.xyz for r in chain(self.paratope, self.epitope) for a in r.atom]),0)
    
    logger.debug("\nepitope center %s.\nparatope center: %s. \nparaepi center: %s" %(repr(self.epi_center), repr(self.pe_center), repr(self.paraepi_center)))

def set_axial_plane(self):
    """get the axial plane(中轴面)of the compelx"""
    
    #first set the epitope and paratope center
    self.set_paraepi_center()

    #vector between epitope center and paratope center
    norm_vec = self.epi_center - self.pe_center
    
    #the point that the plane passes through
    point = self.pe_center
    
    #calculate the axial plane based on the point and the vector
    plane = get_perp_plane(norm_vec, point)

    #set it
    self.axial_plane = plane


def init_complex_axial_plane_trait(self):
    self.set_paraepi_center = MethodType(set_paraepi_center, self)

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

def init_atg_atb_axial_plane_trait(self):
    self.set_atg_atb_center = MethodType(set_atg_atb_center, self)

