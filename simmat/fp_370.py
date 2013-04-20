"""
Pairwise finger print calculator for 370-bits finger print data set
"""

from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

from io import fp370_atg_dataloader, fp370_atb_dataloader
from sim_metric import corr_coef


def callback(c1,c2,value):
    print c1,c2,value

def main():
    c_id_list = complex_ids()
    complex_pairwise_calc(c_id_list, fp370_atb_dataloader, corr_coef, callback = callback)


if __name__ == "__main__":
    main()
