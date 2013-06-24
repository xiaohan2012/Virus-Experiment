import os
from source import make_dataloader, make_single_line_converter
from ve.config import data237_fp175_padded_root


from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

from sim_metric import corr_coef

def callback(c1,c2,value):
    print c1,c2,value

def main(plane_type, atg_or_atb, res_or_tri):
    data_dir = os.path.join(data237_fp175_padded_root, "%s-%s-%s" %(plane_type, atg_or_atb, res_or_tri))

    c_id_list = complex_ids(data_dir)

    dataloader = make_dataloader(data_dir, make_single_line_converter(None))

    complex_pairwise_calc(c_id_list, dataloader, corr_coef, callback = callback)

def usage():
    print """
Usage:
    
    python fp_175_padded.py  complex|residue atg|atb tri|res
    """
    
if __name__ == '__main__':
    import sys
    try:
        plane_type, atg_or_atb, res_or_tri = sys.argv[1:]
    except ValueError:
        usage()
        sys.exit(-1)

    if plane_type not in ["complex", "residue"] or \
       atg_or_atb not in ["atg", "atb"] or \
       res_or_tri not in ["res", "tri"]:
        usage()
        sys.exit(-1)
        
    main(plane_type, atg_or_atb, res_or_tri)