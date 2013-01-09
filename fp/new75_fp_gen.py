import os

from itertools import chain
from collections  import OrderedDict
from glob import glob

from fp_75_gen import get_15bits
from fp_gen2 import Residue, Complex

from ve.util.load_pdb import load_pdb_struct
from ve.util.file import fpstr2file

from ve.config import *
class FP75(object):
    def __init__(self,antigen,antibody):
        self.antigen = antigen
        self.antibody = antibody
        self.fp_15 = {}
        self.fp_first_30 = OrderedDict()
        self.fp_last_30 = OrderedDict()

    def gen_15_bits(self):
        antigen_fp_15 = get_15bits(self.antigen, self.antibody)
        antibody_fp_15 = get_15bits(self.antibody, self.antigen)
        self.fp_15.update(antigen_fp_15)
        self.fp_15.update(antibody_fp_15)

    def gen_first_30(self):
        """result from iterating through the `antigen` side"""
        complex = Complex(self.antigen)
        for res in chain(self.antigen.residue,self.antibody.residue):
            res = Residue(res, complex)
            self.fp_first_30[res.resnum] = res.last_30_bits()

    def gen_last_30(self):
        """result from iterating through the `antibody` side"""
        complex = Complex(self.antibody)
        for res in chain(self.antigen.residue,self.antibody.residue):
            res = Residue(res, complex)
            self.fp_last_30[res.resnum] = res.last_30_bits()

    def gen_fp_str(self):
        numtostr = lambda s:  s if isinstance(s,basestring) else "%.3f" %s

        self.gen_15_bits()
        self.gen_first_30()
        self.gen_last_30()
        self.concat_fps()
        string = ""
        return "\n".join("%d: %s" %(res,','.join(map(numtostr,fp))) \
                for res,fp in self.fp.items())
        

    def gen_fp_dict(self):
        self.gen_15_bits()
        self.gen_first_30()
        self.gen_last_30()
        self.concat_fps()
        return self.fp

    def concat_fps(self):
        self.fp = OrderedDict()

        print self.fp_15.keys()
        print self.fp_first_30.keys()
        print self.fp_last_30.keys()

        for res,fp in self.fp_15.items():
            self.fp[res] = fp + self.fp_first_30[res] + self.fp_last_30[res]


####Real life date processing part#####
def data237_fp_gen(refresh = True):
    path = os.path.join(data237_root ,"splitted_complex","*")
    print path
    for fp in glob(path):
        complex_id = os.path.basename(fp)
        
        #prepare neccessary directory
        if not os.path.exists(data237_fp75_root): os.makedirs(data237_fp75_root)

        fp_path = os.path.join(data237_fp75_root,"%s.fp" %complex_id)

        if not refresh and os.path.exists(fp_path):
            print "%s processed" %complex_id
        else:
            print "processing %s" %complex_id
            atb_path = os.path.join(fp,"antibody.pdb")
            atg_path = os.path.join(fp,"antigen.pdb")
            antibody = load_pdb_struct(atb_path)
            antigen = load_pdb_struct(atg_path)
        
            
            fp = FP75(antigen,antibody)
            fp_str = fp.gen_fp_str()

            fpstr2file(fp_str, fp_path)
            print "fingerprint saved"

if __name__ == "__main__":
    data237_fp_gen(refresh = False)

