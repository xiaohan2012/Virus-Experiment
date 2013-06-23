from ve.fp.complex_util.paraepi import init_io_trait
from ve.util.load_pdb import *

from ve.fp.test.fake_class import TestResidue

class ParaEpiGen(object):
    def __init__(self, complex_dir, c_id):
        
        init_io_trait(self)
        
        #load the antibody and antigen
        self.atg = load_pdb_struct(os.path.join(complex_dir, c_id, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(complex_dir, c_id, "antibody.pdb"), TestResidue)

        self.c_id = c_id
        
        #find the paratope and epitope
        self.find_paratope()
        self.find_epitope()
        
    def gen_paraepi(self):
        self.write_epitope()
        self.write_paratope()
    
if __name__ == "__main__":
    import sys

    from ve.config import data237_raw_complex as raw_complex_path, data237_complex_root as splitted_complex_path

    from ve.fp.complex_util.paraepi import ParatopeNotFoundError, EpitopeNotFoundError
    #get the complex ids
    cids = complex_ids(raw_complex_path)
    print cids
    
    path = splitted_complex_path
    for c_id in cids:
        try:
            obj = ParaEpiGen(path, c_id)
            logger.info("gen paraepi for %s" %c_id)
            obj.gen_paraepi()
        except ParatopeNotFoundError, EpitopeNotFoundError:
            logger.error("Abnormal data %s" %c_id)
    
