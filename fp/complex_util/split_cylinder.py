
import math
import logging
import sys

from ve.fp.fp import BaseComplexFingerprint

from common import *

logger = make_logger("Split Cylinder]")

class Cylinder(object):#the cylinder concept can be factored out
    """
    Splitted Cylinder used in the split cylinder functionality
    """
    
    def __init__(self, center, plane, radius = 10, radius_step = 2, height = 20, height_step = 2):
        self.center = center
        self.plane = plane
        
        self.cylinder_radius  = radius
        self.cylinder_radius_step = radius_step
        self.cylinder_height = height
        self.cylinder_height_step = height_step

        self.cylinder_height_range = (self.cylinder_height/2 - self.cylinder_height , self.cylinder_height/2)
        self.layer_size = self.cylinder_radius / self.cylinder_radius_step
        self.layer_count = self.cylinder_height / self.cylinder_height_step

    def in_which(self, pt):
        """
        (Cylinder, Point) => int

        calculate the point's distance to the axial plane and its distance to the axis, 
        return this point's position index within the cylinder
        """
        dist2axis = self.plane.get_perp_point(pt).dist2point(self.center)

        dist2plane = self.plane.dist2point(pt)
        
        min_height,max_height = self.cylinder_height_range

        #test if it is within the cylinder
        if (dist2plane > min_height and dist2plane < max_height) and dist2axis < self.cylinder_radius:
            layer_idx = self.layer_count/2 + int(math.floor(dist2plane / self.cylinder_height_step))
            ring_idx = int(math.floor(dist2axis / self.cylinder_radius_step))

            #ring index as the first dimension 
            #layer index as the second
            return layer_idx * self.layer_size + ring_idx
            
    def __str__(self):
        return "center: %s, plane: %s, radius: %d by %d, height: %d by %d" %(self.center, self.plane, self.cylinder_radius, self.cylinder_radius_step, self.cylinder_height, self.cylinder_height_step)

class SplitCylinderUtility(object):
    """
    Generic trait class for cylinder splitting method
    """
    
    def __init__(self,  radius = 10, radius_step = 2, height = 20, height_step = 2):
        self.radius = radius
        self.radius_step = radius_step
        self.height = height
        self.height_step = height_step
        self.ring_count = self.radius / self.radius_step * self.height / self.height_step
        
        super(SplitCylinderUtility,self).__init__()

    def make_cylinder(self, center, plane):
        return Cylinder(center, plane, self.radius, self.radius_step, self.height, self.height_step)
        
def make_split_cylinder_method(get_cylinder_func):
    """
    the axial plane used may vary by the models applied, use a function wrapper
    """
    def gen_fp_by_splitting_cylinder(self, bases=[],targets=[], fps = None):
        """
        (list of Residues, list of (int, list of Residue / Triangle)) => ComplexFingerprint
        """
        if fps is None:#use userdefine fingerprint
            fps = BaseComplexFingerprint()
        
        for base in bases:
            cylinder = get_cylinder_func(self, base)
            fp_length = cylinder.layer_size * cylinder.layer_count

            for index_base,iterables in targets:
                for other in iterables:
                    #same as the base, ignore
                    if base == other:
                        continue
                    
                    idx = cylinder.in_which(other.get_center())

                    #if within range
                    if idx:
                        if not fps.has_res(base):#not registered yet
                            fps.add_res(base, fp_length)
                        fps[base][index_base + idx] += 1
        return fps
        
    return gen_fp_by_splitting_cylinder

class SplitCylinderTrait(SplitCylinderUtility):
    def __init__(self,  radius = 20, radius_step = 2, height = 40, height_step = 5, **kwargs):
        super(SplitCylinderTrait,self).__init__(radius = radius, radius_step = radius_step, height = height, height_step = height_step, **kwargs)

    def gen_fp_by_splitting_cylinder(self, bases=[],targets=[]):
        raise NotImplementedError

from types import MethodType
from ve.fp.complex_util.axial_plane import HasAxialPlaneTrait

class SplitCylinderViaComplexPlaneTrait(SplitCylinderTrait, HasAxialPlaneTrait):
    def __init__(self, *args, **kwargs):
        super(SplitCylinderViaComplexPlaneTrait,self).__init__(*args, **kwargs)

        self.gen_fp_by_splitting_cylinder = MethodType(make_split_cylinder_method(lambda self, base: self.make_cylinder(base.get_center(), self.get_axial_plane())), self)

class SplitCylinderViaResiduePlaneTrait(SplitCylinderTrait, HasAxialPlaneTrait):
    def __init__(self, *args, **kwargs):
        super(SplitCylinderViaResiduePlaneTrait,self).__init__(*args, **kwargs)

        self.gen_fp_by_splitting_cylinder = MethodType(make_split_cylinder_method(lambda self, base:
                                                                                  self.make_cylinder(base.get_center(), base.get_axial_plane())), self)

#halt here, a new chapter unfolds
from ve.fp.complex_util.paraepi import FindParaEpiTrait
from ve.fp.complex_util.axial_plane import HasAxialPlaneTrait

class GBCSplitCylinderTrait(SplitCylinderUtility, HasAxialPlaneTrait, FindParaEpiTrait):
    """
    Geometric-center-based complex split cylinder method

    this trait depends on
    1, find paratope, epitope
    2, axial plane(for epitope and paratope center calculation)
    """
            
    def __init__(self, **kwargs):
        logger.info("Register gcb split cylinder trait")
        super(GBCSplitCylinderTrait,self).__init__(**kwargs)

    def gcb_gen_cylinder(self):
        """
        Generate the cylinder geometry data(center and axial plane)

        Point: middle point of the line segment between `paratope geometric center` and `epitope geometric center`
        Plane: the plane 1, which is perpendicular to the above line segment and 2, in which the `Point` lies in

        """
        from ve.fp.geom import Vector, get_perp_plane

        #the line between epitope center and paratope center
        line = Vector(self.get_epi_center() - self.get_para_center())

        #get the paratope and epitope center as the cylinder center

        self.cylinder_center = self.get_paraepi_center()

        #get the axial plane
        self.axial_plane = get_perp_plane(line, self.cylinder_center)

    def get_cylinder_center(self):
        """as the name suggests"""
        if not hasattr(self, "cylinder_center"):
            self.gcb_gen_cylinder()
        return self.cylinder_center

    def get_axial_plane(self):
        """as the name suggests"""
        if not hasattr(self, "axial_plane"):
            self.gcb_gen_cylinder()
        return self.axial_plane

    def gen_fp_by_splitting_cylinder(self, targets):
        """
        geometric-center-based split cylinder method

        (list of Redisue) => Fingerprint

        """
        #init the finger print object
        from ve.fp.fp import HeadlessFingerprint

        fp = HeadlessFingerprint(self.radius / self.radius_step * self.height / self.height_step)

        cylinder = self.make_cylinder(self.get_cylinder_center(), self.get_axial_plane())

        for other in targets:

            idx = cylinder.in_which(other.get_center())
            
            #it's within the cylinder
            if idx:
                fp[idx] += 1
        return fp

    def get_atg_fp_by_split_cylinder(self):
        return self.gen_fp_by_splitting_cylinder(self.atg.residues)

    def get_atb_fp_by_split_cylinder(self):
        return self.gen_fp_by_splitting_cylinder(self.atb.residues)



