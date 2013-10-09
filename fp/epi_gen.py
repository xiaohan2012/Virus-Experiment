from ve.fp.complex_util.paraepi import ParaEpiIOTrait
from ve.util.complex import BaseComplex
from ve.util.load_pdb import *

from ve.fp.test.fake_class import TestResidue

class ParaEpiGen(BaseComplex, ParaEpiIOTrait):
    def gen_paraepi(self, complex_dir, paraepi_dir, only_paratope = False):
        self.find_paratope()
        self.write_paratope(complex_dir, paraepi_dir)
        
        if only_paratope:
            self.find_epitope()
            self.write_epitope(complex_dir, paraepi_dir)
        
            
    
if __name__ == "__main__":
    import sys
    from ve.util.load_pdb import complex_ids, load_complexes
    from ve.fp.complex_util.paraepi  import ParatopeNotFoundError, EpitopeNotFoundError
    from ve.config import data_root
    
    var = sys.argv[1]

    complex_dir = os.path.join(data_root,"three-groups/split-complexes", var)
    
    #get the complex ids
    ids = complex_ids(complex_dir)
    
    print ids

    cs = load_complexes(ids, directory =  complex_dir, complex_cls = ParaEpiGen, residue_cls =  TestResidue)
    
    for c in cs:
        try:
            c.gen_paraepi(complex_dir, os.path.join(data_root,"three-groups/paraepi", var), only_paratope = True)
        except ParatopeNotFoundError:
            print "paratope not found"
        except EpitopeNotFoundError:
            print "epitope not found"
        