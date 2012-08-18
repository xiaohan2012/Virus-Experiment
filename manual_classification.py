import os
from collections import defaultdict
import glob

from customcollections import OrderedDefaultDict
from config import *

def get_166_manual_groups():
    pdb_names = []
    
    for fname in glob.glob(data_src):
        complex_id = os.path.basename(fname).split('.')[0]
        pdb_names.append(complex_id.strip())
    print pdb_names,len(pdb_names)
    
    pdb_fp = os.path.join(boss_root , 'manual_classification_result/166_pdbname.txt')
    type_fp = os.path.join(boss_root , 'manual_classification_result/166_type.txt')
    
    class_d = OrderedDefaultDict(list)
    for name,c_type in zip(open(pdb_fp).readlines(),\
                           open(type_fp).readlines()):
        name = '_'.join(name.strip().split())
        c_type = c_type.strip()
        if name and c_type and name in pdb_names:#not empty line
            class_d[c_type].append(name)
    count = 0                    
    for c,pdbs in class_d.items():
        print c,pdbs
        count += len(pdbs)
    print count        
    return class_d.values()

class PdbGroupRelation(defaultdict):
    def __init__(self,groups):
        defaultdict.__init__(self,dict)
        for g in groups:
            for i in xrange(len(g)):
                for j in xrange( i+1 ,len(g)):
                    pdb_i = g[i]
                    pdb_j = g[j]
                    self[pdb_i][pdb_j] = 1
                    self[pdb_j][pdb_i] = 1
                    
    def is_in_one_group(self,pdb1,pdb2):
        try:
            return self[pdb1][pdb2]
        except KeyError:
            return -1
