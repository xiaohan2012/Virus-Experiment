from itertools import combinations
import numpy as np

class ResTriangle(object):
    def __init__(self,reslist):
        if len(reslist) != 3:
            raise ValueError("The residue number should be 3")
        self.pts = reslist
        self.resnums = sorted([res.resnum for res in reslist])
        self.fingerprint = "-".join(map(str,self.resnums))

    def cal_center(self):
        a = np.array([a.xyz for res in self.pts for a in res.atom])
        self.center = np.average(a,0)

    def __eq__(self,other):
        if isinstance(other,ResTriangle):
            return self.fingerprint  == other.fingerprint
        return False

    def __hash__(self):
        return hash(self.fingerprint)
    
    def __repr__(self):
        #sides_str = ",".join(map(lambda n: "%.2f" %n, (r1.dist2residue(r2) for (r1,r2) in combinations(self.pts,2))))
        return "triangle-%s" %(self.fingerprint)
    
    def __iter__(self):
        return iter(self.pts)


def make_triangles(residues):
    pass

if __name__ == "__main__":
    from ve.util.load_pdb import load_pdb_struct
    from ve.config import *
    data_dir = os.path.join(data237_complex_root,"1SLG_D")
    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"))
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"))

    r1,r2,r3,r4 = list(antigen.residue)[:4]

    rt1 = ResTriangle([r1,r2,r3])
    rt2 = ResTriangle([r3,r1,r2])
    rt3 = ResTriangle([r1,r2,r4])
    rt4 = ResTriangle([r3,r2,r4])
    
    s = set([rt1,rt2,rt3,rt4])
    print s
