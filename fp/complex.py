"""
A bunch of Complex class variants
"""
from ve.util.complex import BaseComplex

from ve.fp.complex_util.triangle import ResidueTriangleTrait

class TriangleComplex(BaseComplex, ResidueTriangleTrait):#seems like redundant
    """
    Complex with triangle structure detection functionality
    """
    
    def __init__(self,complex_id, antigen, antibody, **kwargs):
        super(TriangleComplex,self).__init__(complex_id, antigen, antibody, **kwargs)