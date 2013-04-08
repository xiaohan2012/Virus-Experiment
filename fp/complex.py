"""
A bunch of Complex class variants
"""
from ve.util.complex import BaseComplex

from ve.fp.complex_util.paraepi import init_find_epiparatope_util
from ve.fp.complex_util.triangle import init_triangle_util
from ve.fp.complex_util.propagate_distcache import init_propagate_distcache_util

class TriangleComplex(BaseComplex):
    """
    Complex with triangle structure detection functionality
    """
    
    def __init__(self,complex_id, antigen, antibody):
        BaseComplex.__init__(self,complex_id, antigen,antibody)
        
        #the distance matrix cache needs to be propagated
        init_propagate_distcache_util(self)

        #register find epitope and paratope util
        init_find_epiparatope_util(self)
        
        #register find triangle util
        init_triangle_util(self)

        self.neighbour_threshold = 4

        #find paratope and epitope
        print "finding paratope"
        self.find_paratope()
        print "finding epitope"
        self.find_epitope()
        print len(self.paratope),len(self.epitope)
        
        #find triangle
        self.find_triangles()
