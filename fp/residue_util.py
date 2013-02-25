
from collections import defaultdict

from geom import get_perp_plane
from fp import BaseResidueFingerprint

from types import MethodType

def find_neighbours(self,others,dist=4):
    """find neighbours"""
    for other in others:
        if  self.dist2residue(other) <= dist and self != other:
            self.nbs.append(other)
    return self.nbs

def get_neighbours(self):
    return self.nbs

def dist2residue(self,other):
    return self.distcache.get(self,other)

def init_resdist_util(self):
    self.dist_cache = None
    self.dist2residue = MethodType(dist2residue,self)

def init_neighbour_util(self):
    self.nbs = []
    self.find_neighbours = MethodType(find_neighbours,self)
    self.get_neighbours = MethodType(get_neighbours,self)

def set_axial_plane(self, other_center):
    vec = self.center - other_center
    self.axial_plane = get_perp_plane(vec, self.center)

def init_axial_plane_util(self):
    self.set_axial_plane = MethodType(set_axial_plane,self)


#electric finger print calculation

def gen_electric_fp(self,others):
    bitlength=30
    layer_res_list = defaultdict(list)
    for other in others:
        dist = self.distcache.get(other,self)
        dist_ind =  int(dist / self.dist_step)
        
        #within range
        if dist >= self.dist_min and dist <= self.dist_max:
            layer_res_list[dist_ind].append(other)
    
    fp = BaseResidueFingerprint(self,bitlength)

    for i in xrange(self.layer_count):# at layer i
        h_bond,charged,hydro = 0,0,0

        for res in layer_res_list[i]:# for res in layer i
            code = res.getCode()

            hydro += self.hydro_dict[code]
            charged += self.charged_dict[code]
            h_bond += self.h_bond_dict[code]
        
        #fp for layer i,in the 3 aspects
        fp[i], fp[10+i], fp[20+i] = hydro, charged, h_bond

    self.fp = fp

def gen_electric_fp_triangle(self,others):
    bitlength=30

    layer_tri_list = defaultdict(list)
    for other in others:
        dist = self.distcache.get(other,self)
        dist_ind =  int(dist / self.dist_step)
        
        #within range
        if dist >= self.dist_min and dist <= self.dist_max:
            layer_tri_list[dist_ind].append(other)
    
    fp = BaseResidueFingerprint(self,bitlength)
    
    #for layer,lst in layer_tri_list.items():
        #print layer,lst

    for i in xrange(self.layer_count):# at layer i
        h_bond,charged,hydro = 0., 0., 0.

        for tri in layer_tri_list[i]:# for tri in layer i
            temp1,temp2,temp3=0. , 0. , 0.
            for res in tri:#each residue in the triangle
                code = res.getCode()

                temp1 += self.hydro_dict[code]
                temp2 += self.charged_dict[code]
                temp3 += self.h_bond_dict[code]
            h_bond += temp1 / 3
            charged += temp2 / 3
            hydro += temp3 / 3
        #fp for layer i,in the 3 aspects
        fp[i], fp[10+i], fp[20+i] = hydro, charged, h_bond

    self.fp = fp

def init_electric_fp_util(self,dist_max=20, dist_step=2):
    self.hydro_dict = {'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88}
    self.charged_dict={'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03}
    self.h_bond_dict={'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1}

    self.gen_electric_fp = MethodType(gen_electric_fp,self)
    self.gen_electric_fp_triangle = MethodType(gen_electric_fp_triangle,self)
    
    self.dist_min = 0
    self.dist_max = dist_max
    self.dist_step = dist_step
    self.layer_count = self.dist_max / self.dist_step
