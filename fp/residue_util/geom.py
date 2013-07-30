import numpy as np

from types import MethodType

from ve.fp.geom import Point

class GeometryTrait(object):

    def get_center(self):
        if not hasattr(self,"center"):
            self.center =Point(np.average(np.array([a.xyz for a in self.atom]),0))
        return self.center
    

