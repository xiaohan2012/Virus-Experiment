from ve.util.residue import BaseResidue
from ve.fp.residue_util import init_resdist_util

class TestResidue(BaseResidue):
    def __init__(self, residue):
        BaseResidue.__init__(self, residue)
        init_resdist_util(self)
