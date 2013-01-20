from ve.util.complex import BaseComplex
from ve.util.dist import ResDistCache

from res_triangle import ResTriangle

class TriangleComplex(BaseComplex):
    def __init__(self,antigen, antibody):
        BaseComplex.__init__(self,antigen,antibody)

        self.distcache = ResDistCache()

        self.paratope_threshold = 5
        self.epitope_threshold = 5
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

    def _find_paratope(self):
        """antibody side"""
        for b_r in self.atb.residues:
            for g_r in self.atg.residues:
                #print b_r.dist2residue(g_r)
                dist =b_r.dist2residue(g_r)
                if dist <=self.paratope_threshold:
                    self.paratope.append(b_r)
                    break
        #print len(self.paratope)
                            
    def _find_epitope(self):
        """antigen side"""
        for g_r in self.atg.residues:
            for b_r in self.atb.residues:
                #print g_r.dist2residue(b_r)
                if g_r.dist2residue(b_r) <= self.paratope_threshold:
                    self.epitope.append(g_r)
                    break
        #print len(self.epitope)                    

    def _set_dist_cache(self):
        for r in self.atg.residues:
            r.distcache = self.distcache
        for r in self.atb.residues:
            r.distcache = self.distcache

    def _find_triangles(self):
        self.neighbour_threshold = 4
        for r in self.epitope:
            r.find_neighbours(self.epitope,self.neighbour_threshold)

        triangles = []
        byproduct = []#stores residue tuples that are possible to form triangle with external residue

        for parent in self.epitope:
            for child in parent.get_neighbours():
                temp_list = []
                for grandchild in child.get_neighbours():
                    # and child != grandchild and parent != child
                    if parent != grandchild:
                        tri = ResTriangle([parent,child,grandchild])
                        temp_list.append(tri)
                if not temp_list:#it is empty
                    byproduct.append((parent,child))
                triangles += temp_list
        remaining_residues = set(self.atg.residues)-set(self.epitope)
        for r1,r2 in byproduct:
            for r3 in remaining_residues:
                if r1.dist2residue(r3) <= self.neighbour_threshold and\
                   r2.dist2residue(r3) <= self.neighbour_threshold:
                   tri = ResTriangle([r1,r2,r3])
                   triangles.append(tri)
        print "%d more triangles" %len(byproduct)
        triangles = list(set(triangles))
        print "found %d triangles in total" %len(triangles)

        for t in triangles:
            t.cal_center()

        self.triangles = triangles
