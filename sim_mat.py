"""
load similarity matrix from database with cache functionality
"""
import os
from pickle import dump,load
from numpy import array,zeros

from ve.util.aa2code import *
from ve.config import *

def get_sim_matrix_from_db(code_map , col):
    size = len(code_map.keys())
    mat = zeros((size,size))
    print code_map.items()
    for i1,code1 in code_map.items():
        #print i1
        for i2,code2 in code_map.items():
            res = col.find_one({"complex1":code1,"complex2":code2}) or col.find_one({"complex1":code2,"complex2":code1})
            if res:
                #print i1 , i2 , res["val1"] , res["val2"] , res["val3"]
                mat[i1,i2] = res["val1"] + res["val2"] + res["val3"]
            else: 
                print code1 , code2
                raise ValueError("the value does not exist")                
    return mat

def load_sim_mat(mat_obj_id , force_refresh = False):
    if not os.path.exists(sim_dist_pickle_dir):
        os.mkdir(sim_dist_pickle_dir)

    obj_fp = os.path.join(sim_dist_pickle_dir , "%s.mat" %mat_obj_id)
    
    print obj_fp

    code_map = gen_402_codemap()
    try:
        mat = load(open(obj_fp,'r'))
        print "cached already"
        if force_refresh:#force to refresh
            raise IOError("force to refresh")
    except IOError:
        print "not exist, try create it"

        mat = get_sim_matrix_from_db(code_map , db[mat_obj_id] )
        dump(mat,open(obj_fp,'w'))
    return mat


if __name__ == "__main__":
    load_sim_mat("WILM950103" , force_refresh = True)
