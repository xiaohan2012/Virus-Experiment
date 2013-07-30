from common import *

class ResDistTrait(object):
    """residue pair-wise distance trait"""
    def __init__(self, **kwargs):
        super(ResDistTrait, self).__init__()
        
    def close_enough_to(self, distance, other):
        """(Residue, Residue) -> bool

        return if the two residues are close enought regarding to the distance
        """
        return self.dist_to(other) <= distance

    def dist_to(self,other):
        """
        (Residue, Residue) -> float

        distance to other
        """
        return self.distcache.get(self,other)

