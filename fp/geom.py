from __future__ import division
import numpy as np

class Array(object):
    def __init__(self,data):
        self.d = np.array(data)
        
    def __add__(self,other):
        if isinstance(other,Array):
            return self.d + other.d
        elif isinstance(other,ndarray):
            return self.d + other
        else:
            raise TypeError("`other` is not a Array instance")

    def __sub__(self,other):
        if isinstance(other,Array):
            return self.d - other.d
        elif isinstance(other,np.ndarray):
            return self.d - other
        else:
            raise TypeError("`other` is not a Array instance")

    def __mul__(self,other):            
        if isinstance(other,Array):
            return self.d * other.d
        elif isinstance(other,np.ndarray):
            return self.d * other
        else:
            raise TypeError("`other` is not a Array instance")

    def __pow__(self,other): 
        return Array(self.d ** other)
    
    def __len__(self):
        return len(self.d)

    def __getitem__(self,k):
        return self.d[k]





class Point(Array):
    def __init__(self,data):
        Array.__init__(self,data)
        
    def dist2point(self,point):
        return np.sqrt(np.sum((self - point) ** 2))

    def __repr__(self):
        string = ",".join( "%.2f" %d for d in self)
        return "%dD Point (%s)" %(len(self), string)

class Line(Array):
    """currently 3D only"""
    def __init__(self,start_p,end_p):
        self.sp = start_p
        self.ep = end_p

    def dist2point(self,point):
        return np.linalg.norm(np.cross(point - self.sp, point - self.ep)) /\
               np.linalg.norm(self.sp - self.ep)

class Vector(Array):
    pass

class Plane(Array):
    """currently 3D plane only"""
    def __init__(self,data):
        if len(data) != 4:
            raise ValueError("length `4` required")
        Array.__init__(self,data)
        self.A, self.B, self.C, self.D = self

    def dist2point(self,point):
        return (sum(point * self[:3]) - self[3]) /\
                np.sqrt(np.sum(self[:3] ** 2))

    def get_perp_point(self,ext_point):
        """get perpendicular point of a external point"""
        A,B,C,D = self.A, self.B, self.C, self.D
        abc = np.array([A,B,C])
        x,y,z = ext_point
        k = (np.sum(abc**2)) / (np.sum(abc * ext_point) - D)
        return Point([x - A/k, y - B/k, z - C/k])

    def __repr__(self):
        return "%.2fx + %.2fy + %.2fz = %.2f" %(self.A, self.B, self.C, self.D)

def get_perp_plane(norm_vec,point): 
    A,B,C = norm_vec
    D = np.sum(point * norm_vec)
    return Plane([A,B,C,D])

def get_perp_point(line_point1,line_point2,ext_point):
    x1,y1,z1 = line_point1
    x2,y2,z2 = line_point2
    m,n,p    = ext_point
    t = ((x1-m)*(x1-x2) + (y1-n)*(y1-y2) + (z1-p)*(z1-z2)) /\
        ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    print 
    cx,cy,cz = x1+t*(x2-x1),\
               y1+t*(y2-y1),\
               z1+t*(z2-z1)
    return Point([cx,cy,cz])

if __name__ == "__main__":
    p1 = Point([1,2,3])
    p2 = Point([2,3,4])
    print p1,p2
    print p1 - p2
    print p1 + p2
    print p1 * p2
    print p1.dist2point(p2)

    norm_vec = Vector([1,3,-2])
    point = Point([1,2,3])
    plane = get_perp_plane(norm_vec, point)
    
    print plane

    plane = Plane([1,-2,3,5])
    point = Point([2,3,1])
    print plane
    print "point to plane distance"
    print plane.dist2point(point)

    print "plane perpendicular point"
    perp_point = plane.get_perp_point(point)
    print perp_point
    print perp_point - point
    
    p1 = Point([1,0,1])
    p2 = Point([0,1,1])
    line = Line(p1,p2)
    point = Point([1,1,1])
    print "distance from point to line is",line.dist2point(point)
