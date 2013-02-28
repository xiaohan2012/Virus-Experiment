from ve.util.complex import BaseComplex
from ve.util.dist import ResDistCache

from complex_util import init_find_epiparatope_util, init_triangle_util




class TriangleComplex(BaseComplex):
    def __init__(self,complex_id, antigen, antibody):
        BaseComplex.__init__(self,complex_id, antigen,antibody)

        self.distcache = ResDistCache()

        init_find_epiparatope_util(self)
        init_triangle_util(self)

        self.neighbour_threshold = 4

        self.paratope,self.epitope = [],[]
        print "finding paratope"
        self._find_paratope()
        print "finding epitope"
        self._find_epitope()
        print len(self.paratope),len(self.epitope)

        self._find_triangles()

    def __setattr__(self,name,value):
        self.__dict__[name] = value
        if name == "distcache":
            self._set_dist_cache()

    def _set_dist_cache(self):
        for r in self.atg.residues:
            r.distcache = self.distcache
        for r in self.atb.residues:
            r.distcache = self.distcache

