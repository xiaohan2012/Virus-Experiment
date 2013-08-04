import os
import re

from glob import glob

from ve.config import *

def read_chains(fp = os.path.join(data237_root,"list_237.txt")):
    cfg = {}
    with open(fp,"r") as f:
        for l in f.readlines():
            name,g_chain,b_chain,_ = re.split(r"\s+",l,3)
            cfg[name] = (g_chain.split(","),b_chain.split(","))
    return cfg

def complex_chain_cfg(complex_id,cfg):
    return  cfg[complex_id]

def complex_file_paths(f_dir = os.path.join(data237_root,"cp_pdb","*")):
    for fp in glob(f_dir):
        yield fp

def split(chain_file_path, complex_file_path, splitted_complex_path):
    cfg = read_chains(chain_file_path)
    for fp in complex_file_paths(complex_file_path):
        print fp
        complex_id = os.path.basename(fp)
        complex_id = complex_id.split(".")[0]

        g_chain,b_chain = cfg[complex_id]
        g_lines , b_lines = separate_lines(open(fp), g_chain, b_chain)

        fp_dir = os.path.join(splitted_complex_path ,complex_id)

        if not os.path.exists(fp_dir):
            os.makedirs(fp_dir)

        write2file(g_lines, b_lines, fp_dir)

def separate_lines(file_obj,g_chain, b_chain):
    g_lines,b_lines = [],[]
    for l in file_obj.readlines():
        if not (l.startswith("ATOM") or l.startswith("TER")):
            continue
        chain_name = l[21].upper()
        if chain_name in g_chain:g_lines.append(l)
        elif chain_name in b_chain:b_lines.append(l)
        else:print("Ignore unknown chain %s" %chain_name)
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
    import sys
    var = sys.argv[1]
    split(chain_file_path="../data/three-groups/chains/%s.txt" %var,
          complex_file_path="../data/three-groups/complexes/%s/*" %var,
          splitted_complex_path="../data/three-groups/split-complexes/%s" %var)
