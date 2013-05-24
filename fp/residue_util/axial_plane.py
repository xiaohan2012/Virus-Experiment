from ve.fp.geom import get_perp_plane

from common import *

class AxialPlaneTrait(object):
    """Residue Axial Plane Trait"""
    def calc_axial_plane(self, other_center):
        """
        (Residue, Point) -> Plane
        
        given another point, calculate the axial plane based on self.center
        """

        vec = self.get_center() - other_center
        return get_perp_plane(vec, self.center)

class AtgAtbBasedAxialPlaneTrait(AxialPlaneTrait):
    """ axial plane definition based on atg or atb center"""
    
    def get_axial_plane(self):
        if not hasattr(self, "atg_atb_axial_plane"):
            self.atg_atb_axial_plane = (self.calc_axial_plane(self.c.get_atg_center())
                                        if self.belongs_to_atg()
                                        else self.calc_axial_plane(self.c.get_atb_center()))
                
        return self.atg_atb_axial_plane