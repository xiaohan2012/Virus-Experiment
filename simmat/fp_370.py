"""
Pairwise finger print calculator for 370-bits finger print data set
"""

from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

import source as src
from sim_metric import corr_coef


def callback(c1,c2,value):
    print c1,c2,value

def main():
    filter_list = set(["1UWX_AP"])
    c_id_list = set(complex_ids()) - filter_list

#    complex_pairwise_calc(c_id_list, src.fp370_atg_dataloader, corr_coef, callback = callback)
    complex_pairwise_calc(c_id_list, src.fp370_atb_dataloader, corr_coef, callback = callback)
    
#    complex_pairwise_calc(c_id_list, src.first_110_atg_dataloader, corr_coef, callback = callback)
#    complex_pairwise_calc(c_id_list, src.second_110_atg_dataloader, corr_coef, callback = callback)
#    complex_pairwise_calc(c_id_list, src.last_150_atg_dataloader, corr_coef, callback = callback)

#    complex_pairwise_calc(c_id_list, src.first_110_atb_dataloader, corr_coef, callback = callback)
#    complex_pairwise_calc(c_id_list, src.second_110_atb_dataloader, corr_coef, callback = callback)
#    complex_pairwise_calc(c_id_list, src.last_150_atb_dataloader, corr_coef, callback = callback)

if __name__ == "__main__":
    main()
