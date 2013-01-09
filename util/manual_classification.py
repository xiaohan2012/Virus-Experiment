"""
load how complex files are manually classified into groups 
"""

import os
from collections import defaultdict
import glob

from customcollections import OrderedDefaultDict
from config import *

def get_manual_groups(group_id = "157" ):
    print "generating manual classification and grouping"
    pdb_names = []
    
    for fname in glob.glob(pdb_src):
        complex_id = os.path.basename(fname).split('.')[0]
        pdb_names.append(complex_id.strip())
    #print pdb_names,len(pdb_names)
    
    pdb_fp = os.path.join(data_root , 'manual_classification_result/%s_pdbname.txt' %group_id)
    type_fp = os.path.join(data_root , 'manual_classification_result/%s_type.txt' %group_id)
    
    class_d = OrderedDefaultDict(list)
    for name,c_type in zip(open(pdb_fp).readlines(),\
                           open(type_fp).readlines()):
        name = '_'.join(name.strip().split())
        c_type = c_type.strip()
        if name and c_type and name in pdb_names:#not empty line
            class_d[c_type].append(name)
    count = 0                    
    for c,pdbs in class_d.items():
        #print c,pdbs
        count += len(pdbs)
    #print count        
    return class_d.values()

class PdbGroupRelation(defaultdict):
    """Matrix that indicate whether two pdbs are in one group"""
    def __init__(self,groups):
        defaultdict.__init__(self,dict)
        for g in groups:
            for i in xrange(len(g)):
                for j in xrange( i+1 ,len(g)):
                    pdb_i = g[i]
                    pdb_j = g[j]
                    self[pdb_i][pdb_j] = 1
                    self[pdb_j][pdb_i] = 1
        self.pdbs = [pdb for g in groups for pdb in g]
                    
    def is_in_one_group(self,pdb1,pdb2):
        try:
            return self[pdb1][pdb2]
        except KeyError:
            return -1
if __name__ == "__main__"        :
    get_manual_groups()
