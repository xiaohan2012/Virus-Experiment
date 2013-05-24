from ve.util.residue import BaseResidue
from ve.fp.residue_util.dist import ResDistTrait

class TestResidue(BaseResidue, ResDistTrait):
    def __init__(self, residue, **kwargs):
        super(TestResidue, self).__init__(residue, **kwargs)
