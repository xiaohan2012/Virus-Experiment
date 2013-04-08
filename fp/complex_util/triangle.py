from res_triangle import ResTriangle

def find_triangles(self):
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

def init_triangle_util(self):
    self.neighbour_threshold = 4
    self.find_triangles = MethodType(find_triangles, self)

