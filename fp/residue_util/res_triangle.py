from itertools import combinations
import numpy as np

class ResTriangle(object):
    def __init__(self,reslist):
        if len(reslist) != 3:
            raise ValueError("The residue number should be 3, however %s" %(repr(reslist)))
        self.pts = reslist
        self.res_ids = sorted([res.res_id for res in reslist])
        self.fingerprint = "+".join(self.res_ids)

    def get_center(self):
        if not hasattr(self, "center"):
            self.center = np.average(np.array([a.xyz for res in self.pts for a in res.atom]), 0)
        return self.center
        
    def __eq__(self,other):
        return  self.fingerprint  == other.fingerprint if isinstance(other, self.__class__) else False

    def __hash__(self):
        return hash(self.fingerprint)
    
    def __repr__(self):
        return "triangle-%s" %(self.fingerprint)
    
    def __iter__(self):
        return iter(self.pts)


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
