import numpy as np

from schrodinger.structure import _Residue

from ve.util.dist import *

class BaseResidue(object):
    """wrapper class for _Residue class"""
    def __init__(self,residue=None, **kwargs):
        self.residue = residue
        self.res_id = "%s-%d%s" %(self.chain, self.resnum, self.inscode)
        super(BaseResidue, self).__init__()
        
    def __getattr__(self,name):
        if hasattr(self.residue, name):
            return getattr(self.residue,name)
        else:
            raise AttributeError("No such attribute `%s`" %name)
            
    def __str__(self):
        return "Residue-%s" %(self.res_id)

    def __hash__(self):
        return hash(self.res_id)
        
    def __repr__(self):
        return str(self)

    def belongs_to_atg(self):
        if not hasattr(self, "belongs_to_atg"):
            self.belongs_to_atg = self in self.c.atg
        return self.belongs_to_atg

    def belongs_to_atb(self):
        return not self.belongs_to_atg()