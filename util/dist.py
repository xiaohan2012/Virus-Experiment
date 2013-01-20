import numpy as np
from ve.util.cache import BaseCache

def coord_dist(p1,p2):
    xyz1,xyz2 = np.array(p1),np.array(p2)
    return np.sqrt(np.sum((xyz1 -xyz2)**2))


def atom_dist(a1,a2):
    return coord_dist(a1.xyz, a2.xyz)

def res_ca_dist(r1,r2):
    return coord_dist(r1.ca.xyz, r2.ca.xyz)

def min_dist_func(r1,r2):
    min_dist = float("inf")

    for a1 in r1.atom:
        for a2 in r2.atom:
            a2a_dist = atom_dist(a1,a2)
            if a2a_dist < min_dist:
                min_dist = a2a_dist
    return min_dist

class ResDistCache(BaseCache):
    def __init__(self,dist_func = min_dist_func):
        BaseCache.__init__(self)
        self.dist_func = dist_func

    def set(self,res1,res2,value):
        key = (res1,res2)
        self[key] = value

    def get(self,res1,res2,default = None):
        key1 = (res1,res2)
        key2 = (res2,res1)

        if key1 in self:
            return self[key1]
        if key2 in self:
            return self[key2]
        #print "first time"
        dist = self.dist_func(res1,res2)
        self.set(res1,res2, dist)
        return dist

class FP105DistCache(ResDistCache):
    def __init__(self):
        ResDistCache.__init__(self)

    def get(self,obj1,obj2):
        key1 = (obj1,obj2)
        key2 = (obj2,obj1)

        if key1 in self:
            return self[key1]
        if key2 in self:
            return self[key2]
        
        #no cache found, we need to calculate it
        if hasattr(obj1,"center"):
            p1 = obj1.center
        elif hasattr(obj1,"ca"):
            p1 = obj1.ca.xyz
        else:
            raise AttributeError("%s has neither center nor ca." %obj1)

        if hasattr(obj2,"center"):
            p2 = obj2.center
        elif hasattr(obj2,"ca"):
            p2 = obj2.ca.xyz
        else:
            raise AttributeError("No center nor ca attribute.")

        dist = coord_dist(p1,p2)
        self.set(obj1,obj2, dist)
        return dist
