# -!- coding=utf8 -!-
import numpy as np

from types import MethodType

from geom import *

from common import *
logger = make_logger("Axial Plane")

from ve.fp.complex_util.paraepi import FindParaEpiTrait

class HasAxialPlaneTrait(FindParaEpiTrait):
    """
    calculate axial plane of the complex
    """
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


class AtgAtbAxialPlaneTrait(object):
    def get_atg_center(self):
        if not hasattr(self, "atg_center"):
            pts = [a.xyz for r in self.atg.residues for a in r.atom]
            self.atg_center = np.average(np.array(pts),0)
        return self.atg_center

    def get_atb_center(self):
        if not hasattr(self, "atb_center"):
            pts = [a.xyz for r in self.atb.residues for a in r.atom]
            self.atb_center = np.average(np.array(pts),0)
        return self.atb_center
        
