import numpy as np

from schrodinger.structure import _Residue

from ve.util.dist import *

class BaseResidue(object):
    """wrapper class for _Residue class"""
    def __init__(self,residue):
        self.residue = residue
    
    def __getattr__(self,name):
        if hasattr(self.residue, name):
            return getattr(self.residue,name)
        else:
            raise AttributeError("\"%s\" method not found." %name)
    def __repr__(self):
        return "Residue-%d" %self.resnum


if __name__ == "__main__":
    pass
