import os
import numpy as np
from schrodinger import structure

from util.load_pdb import *

def get_paratope(src):
    antibody = load_pdb_struct(os.path.join(src,"antibody.pdb"))
    antigen = load_pdb_struct(os.path.join(src,"antigen.pdb"))
    
    para_list = []
    
    count = 0
    for atb_res in antibody.residue:
        count += 1
        print "%d / %d" %(count,len(antibody.residue))
        flag = False
        for atb_atom in atb_res.atom:
            for atg_res in antigen.residue:
                res_center = np.average([np.array(a.xyz) for a in atg_res.atom],axis = 0)
                if np.sqrt(np.sum((np.array(atb_atom.xyz) - np.array(res_center)) ** 2)) <= 4.0:
                    flag = True
                    para_list.append(atb_res)
                    print "bingo"
                    break
            if flag:
                #break
        #if flag:
            #continue
    para_path = os.path.join(src,"paratopy.pdb")
    
    print para_path
    with open(para_path, "w") as f:
        a_i = 0#atom index
        for r_i, r in enumerate(para_list):
            for a in r.atom:
                a_i += 1
                f.write("%-6s%5d %-4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6s%6s      %-4s%2s%-2s\n" 
                %("ATOM",#1-6,ATOM
                  a_i,#7-11, atom serial number
                  a.pdbname,#13-16, atom name
                  '',#17
                  r.pdbres.strip(),#18-20, residue name
                  a.chain,#22, chain name
                  r_i,#23-26, residue index
                  '',#27 code for insertion of residues
                  a.x,
                  a.y,
                  a.z,
                  "",
                  "",
                  "",
                  a.element,
                  "",
                ))
    return None
if __name__ == "__main__":
    from config import *
    import glob

    for fname in glob.glob(os.path.join(data_root ,"complex/*")):
        complex_id = os.path.basename(fname) 
        if os.path.exists(os.path.join(data_root,"complex/%s/paratope.pdb" %complex_id)):
            print "%s already done" %complex_id
            continue
        get_paratope(os.path.join(data_root,"complex/%s" %complex_id))
