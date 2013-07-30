from ve.fp.geom import Point
from ve.fp.complex_util.paraepi import FindParaEpiTrait
import numpy as np

class GeometryTrait(FindParaEpiTrait):
    """
    Dependency:
    1, find paratopi and epitope
    """

    def get_geom_center(self):
        """
        get the geometric center of the compelx
        the mean of x y z of all the atoms of the paratope plus the epitope
        """
        from itertools import chain

        #if not computed
        if not hasattr(self, "geom_center"):
            #save the result
            self.geom_center = Point(np.average(np.array([a.xyz for r in chain(self.find_paratope(), self.find_epitope()) for a in r.atom]), axis = 0))

        #return it
        return self.geom_center

    def get_epitope_center(self):
        """
        the epitope center
        """
        if not hasattr(self, "epitope_center"):
            self.epitope_center = Point(np.average(np.array([a.xyz for r in self.find_epitope() for a in r.atom]), axis = 0))
        return self.epitope_center
        
    def get_paratope_center(self):
        """
        the paratope center
        """
        if not hasattr(self, "paratope_center"):
            self.paratope_center = Point(np.average(np.array([a.xyz for r in self.find_paratope() for a in r.atom]), axis = 0))
        return self.paratope_center
