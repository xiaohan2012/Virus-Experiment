import os
from source import make_dataloader, make_single_line_converter
from ve.config import data480_root


from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

from sim_metric import corr_coef

fp_dir = os.path.join(data480_root, "fp_aaindex_if_padded")

first_3_dataloader = make_dataloader(fp_dir, make_single_line_converter(slice(0, 429)))
second_3_dataloader = make_dataloader(fp_dir, make_single_line_converter(slice(429, 429+447)))
last_15_dataloader = make_dataloader(fp_dir, make_single_line_converter(slice(429+447, 3021)))


def callback(c1,c2,value):
    print c1,c2,value
    

def gen_pairwise_dist():
    data_dir = os.path.join(data480_root, "fp_aaindex_if_padded")
    dataloader = make_dataloader(data_dir, make_single_line_converter(None))

    c_id_list = complex_ids(data_dir)
    
    complex_pairwise_calc(c_id_list, dataloader, corr_coef, callback = callback)


def print_splitted_fp_in_csv(which = 0):
    """
    split 360 bits into parts
    """
    data_dir = os.path.join(data480_root, "fp_aaindex_if_padded")

    delimiter = ","
    cids = complex_ids(data_dir)
    
    name = ("atg 3 bits", "atb 4 bits", "15 bits")
    dls = [first_3_dataloader, second_3_dataloader, last_15_dataloader]

    print("cid, %s" %name[which])
    for cid in sorted(cids):
        print("%s,%s" %(cid,
                        delimiter.join(map(lambda d: "%.2f" %d, dls[which](cid)))))

def usage():
    print """
python fp_aaindex_if_padded.py printfp 0|1|2
Or
python fp_aaindex_if_padded.py distmat
    """
if __name__ == '__main__':
    import sys
    oper = sys.argv[1]
    if oper in ("printfp", "distmat"):
        if oper is "distmat":
            gen_dist_matrix()
        else:
            which = int(sys.argv[2])
            print_splitted_fp_in_csv(which)
    else:
        usage()
