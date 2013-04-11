from __future__ import division
import numpy as np

class Array(object):
    def __init__(self,data):
        self.d = np.array(data)
        
    def __add__(self,other):
        return self.__class__(self.d + other.d)
    
    def __sub__(self,other):
        return self.__class__(self.d - other.d)
    
    def __mul__(self,other):            
        if isinstance(other, float) or isinstance(other, int):
            return self.__class__(self.d * other)
        elif isinstance(other, np.ndarray):
            return self.__class__(self.d * other)
        else:
            return self.__class__(self.d * other.d)
        
    def __getattr__(self, name):
        if hasattr(self.d, name):
            return getattr(self.d, name)

    def __pow__(self,other): 
        return self.__class__(self.d ** other)

    def __len__(self):
        return self.d

    def __getitem__(self,k):
        return self.d[k]

    def __eq__(self, other):
        return self.d == other.d



class Point(Array):
    def __init__(self,data):
        Array.__init__(self,data)
        
    def dist2point(self,point):
        return np.sqrt(np.sum((self - point) ** 2))

    def __repr__(self):
        string = ",".join( "%.7f" %d for d in self)
        return "Point (%s)" %(string)

class Line(Array):
    """
    Line class in the 3D space
    """
    def __init__(self,start_p,end_p):
        """
        (Point, Point) => Line
        
        given the start point and end point, construct a line
        """
        self.sp = start_p
        self.ep = end_p

    def dist2point(self,point):
        """
        (Point) => float
        the distance between a point and itself
        """
        return np.linalg.norm(np.cross(point - self.sp, point - self.ep)) /\
               np.linalg.norm(self.sp - self.ep)

class Vector(Array):
    pass

class Plane(Array):
    """
    Plane representation in the 3D space
    """

    def __init__(self,data):
        """
        ([float, float, float, float]) => Plane
        given the a b c d parameter of a plane and construct a plane
        """
        if len(data) != 4:
            raise ValueError("length `4` required")
        Array.__init__(self,data)
        self.A, self.B, self.C, self.D = self

    def dist2point(self,point):
        """(Point) => float
        distance between a Point and the plane"""
        return (sum(point * self[:3]) - self[3]) /\
                np.sqrt(np.sum(self[:3] ** 2))

    def get_perp_point(self,ext_point):
        """
        (Point) => Point
        
        given a point, get perpendicular point of it
        """
        A,B,C,D = self.A, self.B, self.C, self.D
        abc = np.array([A,B,C])
        x,y,z = ext_point
        den = (np.sum(abc * ext_point) - D)
        if den == 0.: #the point is on the plane already!
            return Point(ext_point)
        else:
            k = (np.sum(abc**2)) / den
            return Point([x - A/k, y - B/k, z - C/k])

    def __repr__(self):
        return "%.2fx + %.2fy + %.2fz = %.2f" %(self.A, self.B, self.C, self.D)

def get_perp_plane(norm_vec,point): 
    """
    (Vector, Point) => Plane
    calculate the plane is perpendicular to norm_vec as its normal vector and also passes through `point`
    """

    A,B,C = norm_vec
    D = np.sum(point * norm_vec)
    return Plane([A,B,C,D])

def get_perp_point(line_point1,line_point2,ext_point):
    x1,y1,z1 = line_point1
    x2,y2,z2 = line_point2
    m,n,p    = ext_point
    t = ((x1-m)*(x1-x2) + (y1-n)*(y1-y2) + (z1-p)*(z1-z2)) /\
        ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

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

    plane = Plane([1,-2,3,100])
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
