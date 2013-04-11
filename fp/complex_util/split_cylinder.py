import math
import logging
import sys

from types import MethodType

from ve.fp.fp import BaseComplexFingerprint

logging.basicConfig( stream=sys.stderr )
logging.getLogger("Split Cylinder").setLevel( logging.DEBUG )
logger = logging.getLogger("Split Cylinder")

def get_idx(self,dist2plane,dist2center):
    """
    (float, float) => int

    given the point's distance to the axial plane and its distance to the axis, 
    return this point's position index within the cylinder
    """
    min_height,max_height = self.cylinder_height_range
    
    #test if it is within the cylinder
    if (dist2plane > min_height and dist2plane < max_height) and dist2center < self.cylinder_radius:
        layer_idx = self.layer_count/2 + int(math.floor(dist2plane / self.cylinder_height_step))
        ring_idx = int(math.floor(dist2center / self.cylinder_radius_step))
        
        #ring index as the first dimension 
        #layer index as the second
        return layer_idx * self.layer_size + ring_idx

def split_cylinder_method(self, name="",bases=[],targets=[]):
    """
    sample input data
    "name":"epitope triangles",
    "bases":triangles,
    "targets":[(0,triangles),(50,atb.residues)]
    """
    res_fp_length=50


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

def config(self, radius = 10, radius_step = 2, height = 20, height_step = 2):#verbose naming
    self.cylinder_radius  = radius
    self.cylinder_radius_step = radius_step
    self.cylinder_height = height
    self.cylinder_height_step = height_step

    self.cylinder_height_range = (self.cylinder_height/2 - self.cylinder_height , self.cylinder_height/2)
    self.layer_size = self.cylinder_radius / self.cylinder_radius_step
    self.layer_count = self.cylinder_height / self.cylinder_height_step

def init_split_cylinder_util(self,radius = 10, 
                             radius_step = 2,
                             height = 20,
                             height_step = 2):
    config(self, radius, radius_step, height, height_step)
    self.get_fp_generic = MethodType(split_cylinder_method, self)
    self.get_idx = MethodType(get_idx, self)

def gcb_gen_cylinder(self):
    """
    Generate the cylinder data(center and axial plane) related to the geometric-center-based method

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
    
def gcb_split_cylinder(self, targets):
    """
    geometric-center-based split cylinder method

    (list of Redisue) => Fingerprint
    
    """
    self.get_idx = MethodType(get_idx, self)

    #init the finger print object
    from ve.fp.fp import HeadlessFingerprint
    
    fp = HeadlessFingerprint(self.cylinder_radius / self.cylinder_radius_step * self.cylinder_height / self.cylinder_height_step)

    for other in targets:
        perp_point = self.get_axial_plane().get_perp_point(other.center)
            
        #distance between target's center and cylinder axial plane
        dist2plane = self.get_axial_plane().dist2point(other.center)
        
        #distance between target's center and the cylinder axis
        dist2center = perp_point.dist2point(self.get_cylinder_center())
        
        #get the space index of the target's center
        idx = self.get_idx(dist2plane, dist2center)
        
        #it's within the cylinder
        if idx:
            fp[idx] += 1
    return fp

def get_atg_fp_by_split_cylinder(self):
    return self.gen_fp_by_splitting_cylinder(self.atg.residues)

def get_atb_fp_by_split_cylinder(self):
    return self.gen_fp_by_splitting_cylinder(self.atb.residues)

def init_gcb_split_cylinder_trait(self, radius = 20, radius_step = 2, height = 40, height_step = 5):
    """
    Geometric-center-based complex split cylinder method

    this trait depends on
    1, find paratope, epitope
    2, axial plane(for epitope and paratope center calculation)
    """
    #find paretope and epitope trait
    from ve.fp.complex_util.paraepi import init_find_epiparatope_trait
    init_find_epiparatope_trait(self)
    
    #axial plane trait
    from ve.fp.complex_util.axial_plane import init_complex_axial_plane_trait
    init_complex_axial_plane_trait(self)
    
    
    #cylinder parameter configuration
    config(self, radius, radius_step, height, height_step)

    
    logger.info("Register gcb split cylinder trait")
    
    self.gen_fp_by_splitting_cylinder = MethodType(gcb_split_cylinder, self)
    self.gcb_gen_cylinder = MethodType(gcb_gen_cylinder, self)
    
    self.get_atg_fp_by_split_cylinder = MethodType(get_atg_fp_by_split_cylinder, self)
    self.get_atb_fp_by_split_cylinder = MethodType(get_atb_fp_by_split_cylinder, self)
    
    #getter method for axial plane and cylinder center
    self.get_cylinder_center = MethodType(get_cylinder_center, self)
    self.get_axial_plane = MethodType(get_axial_plane, self)
