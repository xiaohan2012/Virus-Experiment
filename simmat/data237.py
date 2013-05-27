import os
from source import make_dataloader, make_single_line_converter, make_fp_str_saver, load_cids
from ve.config import data237_root

fp370_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(0, 370)))
fp370_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(0, 370)))

first_110_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(0, 110)))
second_110_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(110, 220)))
last_150_atg_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atg"), make_single_line_converter(slice(220, 370)))

first_110_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(0, 110)))
second_110_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(110, 220)))
last_150_atb_dataloader = make_dataloader(os.path.join(data237_root, "fp_370_atb"), make_single_line_converter(slice(220, 370)))


first_110_atg_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atg_first_110"))
second_110_atg_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atg_second_110"))
last_150_atg_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atg_last_150"))

first_110_atb_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atb_first_110"))
second_110_atb_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atb_second_110"))
last_150_atb_datasaver = make_fp_str_saver(os.path.join(data237_root, "fp_370_atb_last_150"))

def print_splitted_fp_in_csv(use_atg = True):
    """
    split 360 bits into parts
    """
    delimiter = " "
    cids = load_cids("data/fp_370_atb.txt")
    print("cid, first 110, second 110, last 150")
    for cid in sorted(cids):
        if use_atg:
            print("%s,%s,%s,%s" %(cid,
                                  delimiter.join(map(lambda d: "%.2f" %d, first_110_atg_dataloader(cid))),
                                  delimiter.join(map(lambda d: "%.2f" %d, second_110_atg_dataloader(cid))),
                                  delimiter.join(map(lambda d: "%.2f" %d, last_150_atg_dataloader(cid))),
                              ))
        else:
            print("%s,%s,%s,%s" %(cid,
                                  delimiter.join(map(lambda d: "%.2f" %d, first_110_atb_dataloader(cid))),
                                  delimiter.join(map(lambda d: "%.2f" %d, second_110_atb_dataloader(cid))),
                                  delimiter.join(map(lambda d: "%.2f" %d, last_150_atb_dataloader(cid))),
                              ))
if __name__ == '__main__':
    print_splitted_fp_in_csv(True)