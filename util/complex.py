import numpy as np
from ve.fp.geom import Point
class BaseComplex(object):
    def __init__(self,complex_id = "", antigen = None, antibody = None, **kwargs):
        self.atg = antigen
        self.atb = antibody
        self.c_id = complex_id
        
        #set residue's belonging complex
        from itertools import chain
        for r in chain(self.atg.residues, self.atb.residues):
            r.c = self

        self.residues = list(self)
        
        super(BaseComplex,self).__init__(**kwargs)
        
    def __iter__(self):
        from itertools import chain
        return chain(self.atg.residues, self.atb.residues)

    def get_atg_center(self):
        if not hasattr(self, "atg_center"):
            self.atg_center = Point(np.average(np.array([a.xyz for r in self.atg.residues for a in r.atom]),0))
        return self.atg_center
        
    def get_atb_center(self):
        if not hasattr(self, "atb_center"):
            self.atb_center = Point(np.average(np.array([a.xyz for r in self.atb.residues for a in r.atom]),0))
        return self.atb_center

    def get_res_from_resids(self, res_ids):
        """(CacheTrait, list of str)
        return the matching list of residues
        """
        from itertools import chain
        return [r for r in self if r.res_id in res_ids]

    def get_res_from_resid(self, res_id):
        """(CacheTrait, str) -> Residue
        return the first matching Residue
        """
        
        from itertools import chain
        return filter(lambda r: r.res_id == res_id, self)[0]

    def __str__(self):
        return "%s: %s" %(self.__class__, self.c_id)

    def __repr__(self):
        return str(self)