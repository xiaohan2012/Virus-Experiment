from types import MethodType

from ve.fp.geom import Point

import numpy as np

def get_geom_center(self):
    """
    get the geometric center of the compelx
    the mean of x y z of all the atoms of the complex
    """
    from itertools import chain
    
    #if not computed
    if self.geom_center is None:
        #save the result
        self.geom_center = Point(np.average(np.array([a.xyz for r in chain(self.find_paratope(), self.find_epitope()) for a in r.atom]), axis = 0))

    #return it
    return self.geom_center

def init_geom_trait(self):
    """
    Dependency:
    1, find paratopi and epitope
    """
    self.get_geom_center = MethodType(get_geom_center, self)
    self.geom_center = None
    
    from ve.fp.complex_util.paraepi import init_find_epiparatope_trait
    
    init_find_epiparatope_trait(self)

    

