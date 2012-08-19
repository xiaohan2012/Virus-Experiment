import os
import glob

def get_codes_from_file(pdb_src):
    codes = {}
    for i,pdb_path in enumerate(glob.glob(pdb_src)):
        code= os.path.basename(pdb_path).split('.')[0] 
        codes[i] = code
    
    return codes

def get_inv_codes_from_file(pdb_src):
    codes = {}
    for i,pdb_path in enumerate(glob.glob(pdb_src)):
        code= os.path.basename(pdb_path).split('.')[0] 
        codes[code] = i
    
    return codes
