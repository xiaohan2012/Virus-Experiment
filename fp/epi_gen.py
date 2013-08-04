from ve.fp.complex_util.paraepi import ParaEpiIOTrait
from ve.util.complex import BaseComplex
from ve.util.load_pdb import *

from ve.fp.test.fake_class import TestResidue

class ParaEpiGen(BaseComplex, ParaEpiIOTrait):
    def gen_paraepi(self, complex_dir, paraepi_dir):
        self.find_paratope()
        self.find_epitope()

        self.write_epitope(complex_dir, paraepi_dir + "/epitope")
        self.write_paratope(complex_dir, paraepi_dir + "/paratope")
    
if __name__ == "__main__":
    import sys
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.fp.complex_util.paraepi  import ParatopeNotFoundError, EpitopeNotFoundError

    var = "HIV"

    complex_dir = "../data/three-groups/split-complexes/%s" %var
    
    #get the complex ids
    ids = complex_ids(complex_dir)
    print ids

    cs = load_complexes(ids, directory =  complex_dir, complex_cls = ParaEpiGen)
    
    for c in cs:
        c.gen_paraepi(complex_dir, "../data/three-groups/paraepi/%s" %var)
