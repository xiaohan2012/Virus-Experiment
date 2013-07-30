import os
from source import make_dataloader, make_single_line_converter
from ve.config import data480_root


from ve.util.calculator import complex_pairwise_calc
from ve.util.load_pdb import complex_ids

from sim_metric import corr_coef



def callback(c1,c2,value):
    print c1,c2,value
    

def gen_pairwise_dist(res_or_tri):
    data_dir = os.path.join(data480_root, "fp_aaindex_if_padded" if res_or_tri == "res" else "fp_aaindex_if_padded_tri")
    dataloader = make_dataloader(data_dir, make_single_line_converter(None))

    c_id_list = complex_ids(data_dir)
    
    complex_pairwise_calc(c_id_list, dataloader, corr_coef, callback = callback)

def print_splitted_fp_in_csv(res_or_tri, which = 0):
    """
    split 360 bits into parts
    """
    data_dir = os.path.join(data480_root, "fp_aaindex_if_padded" if res_or_tri == "res" else "fp_aaindex_if_padded_tri")

    delimiter = ","
    cids = complex_ids(data_dir)
    
    name = ("%s 3 bits" %(res_or_tri), "atb 4 bits", "15 bits")

    from ve.fp.complex_util.padding import PaddedComplexFingerPrint, OverallSpatialDistribution

    atg_dist, atb_dist, tri_dist =  OverallSpatialDistribution.from_cache()

    first_3_count = sum(atg_dist.values() if res_or_tri == "res" else tri_dist.values()) * 3
    second_3_count = sum(atb_dist.values()) * 3
    last_15_count = sum(atg_dist.values()) * 15
    
    first_3_dataloader = make_dataloader(data_dir, make_single_line_converter(slice(0, first_3_count)))
    second_3_dataloader = make_dataloader(data_dir, make_single_line_converter(slice(first_3_count, first_3_count + second_3_count)))
    last_15_dataloader = make_dataloader(data_dir, make_single_line_converter(slice(first_3_count + second_3_count, first_3_count + second_3_count + last_15_count)))

    dls = [first_3_dataloader, second_3_dataloader, last_15_dataloader]

    print("cid, %s" %name[which])
    for cid in sorted(cids):
        print("%s,%s" %(cid,
                        delimiter.join(map(lambda d: "%.2f" %d, dls[which](cid)))))

def usage():
    print """
python fp_aaindex_if_padded.py printfp tri|res 0|1|2
Or
python fp_aaindex_if_padded.py distmat
    """
if __name__ == '__main__':
    import sys
    oper = sys.argv[1]
    tri_or_res = sys.argv[2]

    if oper in ("printfp", "distmat"):
        if oper == "distmat":
            gen_pairwise_dist(tri_or_res)
        else:
            which = int(sys.argv[3])
            print_splitted_fp_in_csv(tri_or_res, which)
    else:
        usage()
