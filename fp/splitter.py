import os
import re

from glob import glob

from config import *

def read_237_config(fp = os.path.join(data237_root,"list_237.txt")):
    cfg = {}
    with open(fp,"r") as f:
        for l in f.readlines():
            name,g_chain,b_chain,_ = re.split(r"\s+",l,3)
            key = "%s_%s" %(name,g_chain)
            cfg[key] = (list(g_chain),list(b_chain))
    return cfg

def complex_chain_cfg(complex_id,cfg):
    return  cfg[complex_id]

def complex_file_paths(f_dir = os.path.join(data237_root,"cp_pdb","*")):
    for fp in glob(f_dir):
        yield fp

def split():
    cfg = read_237_config()
    for fp in complex_file_paths():
        print fp
        complex_id = os.path.basename(fp)
        complex_id = complex_id.split(".")[0]

        g_chain,b_chain = complex_chain_cfg(complex_id,cfg)
        g_lines , b_lines = filter_file(open(fp), g_chain, b_chain)

        fp_dir = os.path.join(data237_root,"splitted_complex" ,complex_id)

        if not os.path.exists(fp_dir):
            os.makedirs(fp_dir)

        write2file(g_lines, b_lines, fp_dir)

def filter_file(file_obj,g_chain, b_chain):
    g_lines,b_lines = [],[]
    for l in file_obj.readlines():
        if not l.startswith("ATOM"):
            continue
        chain_name = l[21]
        if chain_name in g_chain:g_lines.append(l)
        elif chain_name in b_chain:b_lines.append(l)
        else:raise ValueError("Unknown chain %s" %chain_name)
    return g_lines,b_lines

def write2file(g_lines,b_lines,fp_dir):
    g_fp = os.path.join(fp_dir,"antigen.pdb")
    b_fp = os.path.join(fp_dir,"antibody.pdb")
    
    #antigen to file
    with open(g_fp,"w") as f:
        for l in g_lines:
            f.write(l)
    #antibidy to file
    with open(b_fp,"w") as f:
        for l in b_lines:
            f.write(l)
if __name__ == "__main__":
    split()
