import os

from fp_75_gen import get_15bits
from fp_80 import ComplexSingle as ComplexSingle80, ComplexDual as ComplexDual80, Residue as Residue80

from ve.util.load_pdb import load_pdb_struct, complex_ids

from ve.config import data237_fps808015_root, data237_complex_root

def gen_15_bits(self):
    return get_15bits(self.atg, self.atb)

def gen_fp_to_file(self, complex_id, file_dir=data237_fps808015_root, fp_type="single"):
    atg_fp = self.gen_antigen_fp()
    atb_fp = self.gen_antibody_fp()
    bit_15 = self.gen_15_bits()
    
    file_dir = os.path.join(file_dir, fp_type, complex_id)
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)

    atb_fp.tofile(os.path.join(file_dir,"atb.fp" ))
    atg_fp.tofile(os.path.join(file_dir,"atg.fp"))
    bit_15.tofile(os.path.join(file_dir,"15bits.fp"))

    print file_dir, "saved"

#building new classes
ComplexSingle = type("ComplexSingle", 
                     (ComplexSingle80,),
                     {"gen_15_bits":gen_15_bits,
                      "gen_fp_to_file":gen_fp_to_file,})

ComplexDual = type("ComplexDual", 
                     (ComplexDual80,),
                     {"gen_15_bits":gen_15_bits,
                      "gen_fp_to_file":gen_fp_to_file,})
        

def gen_fp(complex_id):            
    data_dir = os.path.join(data237_complex_root,complex_id)

    antigen = load_pdb_struct(os.path.join(data_dir,"antigen.pdb"),residue_cls = Residue80)
    antibody = load_pdb_struct(os.path.join(data_dir,"antibody.pdb"),residue_cls = Residue80)

    c = ComplexSingle(complex_id, antigen, antibody)
    
    c.gen_fp_to_file(complex_id, fp_type="single")

    c = ComplexDual(complex_id, antigen, antibody)

    c.gen_fp_to_file(complex_id, fp_type="double")

def already_processed(complex_id):
    file_types = ["atg.fp","15bits.fp","atb.fp"]
    calc_types = ["single", "double"]
    for ct in calc_types:
        for ft in file_types:
            file_dir = os.path.join(data237_fps808015_root, ct, complex_id, ft)
            print file_dir
            if not os.path.exists(file_dir): return False#not processed completely
    return True#processed completely
    
if __name__ == "__main__":
    ids = ["1SLG_D", "1N4X_L", "1JV5_A", "1STS_B"]
    
    for c_id in complex_ids(path = data237_complex_root):
        if already_processed(c_id):
            print "already processed", c_id
            continue
        else:
            print "processing", c_id
            try:
                gen_fp(c_id)
            except:
                import sys
                sys.stderr.write("error encountered processing %s" %c_id)
                continue
    
