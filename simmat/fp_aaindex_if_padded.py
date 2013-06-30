import os
from source import make_dataloader, make_single_line_converter
from ve.config import data480_root


from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

from sim_metric import corr_coef

def callback(c1,c2,value):
    print c1,c2,value
    

def gen_pairwise_dist():
    data_dir = os.path.join(data480_root, "fp_aaindex_if_padded")
    dataloader = make_dataloader(data_dir, make_single_line_converter(None))

    c_id_list = complex_ids(data_dir)
    
    complex_pairwise_calc(c_id_list, dataloader, corr_coef, callback = callback)


if __name__ == '__main__':
    gen_dist_matrix()
