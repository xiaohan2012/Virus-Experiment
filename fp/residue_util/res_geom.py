import numpy as np

from types import MethodType

from ve.fp.geom import Point

def set_center(self):
    self.center =Point(np.average(np.array([a.xyz for a in self.atom]),0))

def init_geom_trait(self):
    self.set_center = MethodType(set_center, self)
    self.set_center()

