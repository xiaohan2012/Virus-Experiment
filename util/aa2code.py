import os
import glob
import re

from config import *

def gen_402_codemap():
    return get_codes_from_file(pdb_src)
def gen_402_inv_codemap():
    return get_inv_codes_from_file(pdb_src)

def gen_166_codemap():
    return dict((i,line.strip())
                    for i,line in enumerate(open(os.path.join(chain_list_path )).readlines()))

def gen_166_inv_codemap():
    return dict((line.strip(),i)
                    for i,line in enumerate(open(os.path.join(chain_list_path )).readlines()))

def get_codes_from_file(pdb_src):
    return dict((i,os.path.basename(pdb_path).split('.')[0]) 
                    for i,pdb_path in enumerate(glob.glob(pdb_src)))

def get_inv_codes_from_file(pdb_src):
    return dict((os.path.basename(pdb_path).split('.')[0] , i)
                    for i,pdb_path in enumerate(glob.glob(pdb_src)))

def get_rmsd_codemap():
    f = open(rmsd_row_path,"r")
    return dict((i,re.findall(r"(\w{4}_\w)",line)[0])
                    for i,line in enumerate(f.readlines()))

def get_rmsd_inv_codemap():
    f = open(rmsd_row_path,"r")
    return dict((re.findall(r"(\w{4}_\w)",line)[0],i)
                    for i,line in enumerate(f.readlines()))

