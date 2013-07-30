from dist import ResDistTrait
from common import *

class FindNeighbourTrait(ResDistTrait):
    """
    find neighboring residues trait
    
    Dependency:
    1, distance between residue trait
    """

    def __init__(self, triangle_threshold = 4, **kwargs):
        self.triangle_threshold = triangle_threshold
        super(FindNeighbourTrait,self).__init__(**kwargs)

    def get_neighbours(self, others, threshold):
        """get the residue neighbours from others's set'"""
        return [other 
                for other in others
                if self.dist_to(other) <= threshold and self != other]

    def get_neighbours_in_epitope(self):
        """find neighbouring residues in epitope region"""
        if not hasattr(self, "epi_nbs"):
            self.epi_nbs = self.get_neighbours(self.c.epitope, self.triangle_threshold)
        return self.epi_nbs

    def may_form_triangle_with(self, other):
        """(Residue, Residue) -> bool

        return if other is within the triangle formation range with the current residue
        """
        return self.close_enough_to(self.triangle_threshold, other)


