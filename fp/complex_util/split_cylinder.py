import math

from types import MethodType

from ve.fp.fp import BaseComplexFingerprint


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
        self.get_idx = MethodType(get_idx, self)

def gcb_gen_cylinder(self):
    """
    Generate the cylinder data(center and axial plane) related to the geometric-center-based method

    Point: middle point of the line segment between `paratope geometric center` and `epitope geometric center`
    Plane: the plane 1, which is perpendicular to the above line segment and 2, in which the `Point` lies in

    """
    #the line between epitope center and paratope center
    line = Line(self.epi_center, self.pe_center)
    
    #get the paratope and epitope center as the cylinder center

    self.cylinder_center =self.paraepi_center
    
    #get the axial plane
    from ve.fp.geom import get_perp_plane
    self.axial_plane = get_perp_plane(line, self.cylinder_center)
    
    
def gcb_split_cylinder(self, targets):
    """
    geometric-center-based split cylinder method

    (Plane, Point, list of Redisue) => Fingerprint
    
    """
    from ve.fp.fp import GCBCFingerprint

    for index_base,iterables in targets:
        for other in iterables:
            perp_point = self.axial_plane.get_perp_point(other.center)
            
            #distance between target's center and cylinder axial plane
            dist2plane = axial_plane.dist2point(other.center)

            #distance between target's center and the cylinder axis
            dist2center = perp_point.dist2point(self.cylinder_center)
            
            #get the space index of the target's center
            idx = self.get_idx(dist2plane, dist2center)
            
            #if it's within the cylinder
            if idx
                fp[index_base + idx] += 1
        return fp


def init_gcb_split_cylinder_trait(self):
    """
    Geometric-center-based complex split cylinder method

    this trait depends on
    1, find paratope, epitope
    2, axial plane(for epitope and paratope center calculation)
    """

    self.split_cylinder = MethodType(gcb_split_cylinder, self)
    self.gen_cylinder = MethodType(gcb_gen_cylinder, self)
    
    #init the finger print object
    self.fp = GCBCFingerprint(50)

    
