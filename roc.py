"""
provides ROC related functionality
1, test performing
2,roc curve plotting(using yard)
3,auc statistics(using yard)
"""
import os
import glob
from collections import OrderedDict
from customcollections import OrderedDefaultDict
from numpy import linspace,min,max
from pickle import dump,load
import matplotlib.pyplot as plt

from util.aa2code import *
from util.manual_classification import *
from config import *

from dist_mat import DistanceMatrix

def perform_roc_test(groups , mat,roc_step_count = 11):
    test_result = OrderedDict()
    pdbs = groups.keys()
    for cutoff in linspace(0,1,roc_step_count):
        test_result[cutoff] = OrderedDefaultDict(int)
        for i in xrange(len(pdbs)):
            for j in xrange(i,len(pdbs)):
                pdb_i = pdbs[i]
                pdb_j = pdbs[j]
                dist = mat.get_distance_between(pdb_i,pdb_j)
                if dist > cutoff and groups.is_in_one_group(pdb_i,pdb_j) == 1:
                    test_result[cutoff]["tp"] += 1
                elif dist > cutoff and groups.is_in_one_group(pdb_i,pdb_j) == -1:                     
                    test_result[cutoff]["fp"] += 1
                elif dist <= cutoff and groups.is_in_one_group(pdb_i,pdb_j) == 1:
                    test_result[cutoff]["fn"] += 1
                elif dist <= cutoff and groups.is_in_one_group(pdb_i,pdb_j) == -1:
                    test_result[cutoff]["tn"] += 1
                else:    
                    raise RuntimeError("other kind of error")
    cutoff_array = []                
    sensitivity_array = []
    specificity_array = []
    for cutoff,stat in test_result.items():
        temp = test_result[cutoff]
        temp["sensitivity"] = temp["tp"] / float(temp["tp"] + temp["fn"])
        temp["specificity"] = temp["fp"] / float(temp["tn"] + temp["fp"])

        sensitivity_array.append(temp["sensitivity"])
        specificity_array.append(temp["specificity"])
        cutoff_array.append(cutoff)

    return  specificity_array , sensitivity_array , cutoff_array

def generate_yard_file(groups , mat , yard_file_name = "yard_output.txt"):
    with open(yard_file_name ,'w') as f:
        f.write("output\tmethod1\n")
        for i in xrange(len(groups.keys())):
            for j in xrange(i + 1,len(groups.keys())):
                pdb_i = pdbs[i]
                pdb_j = pdbs[j]
                f.write("%f\t%f\n" %(1 if groups.is_in_one_group(pdb_i,pdb_j) == 1 else -1,\
                                    mat.get_distance_between(pdb_i,pdb_j)))

def print_auc(i_path):
    cmd = "yard-auc %s" %i_path
    os.system(cmd)

def draw_roc(x_arr  , y_arr ):
    plt.plot(x_arr,y_arr)
    plt.show()        
    
if __name__ == "__main__":
    groups = get_166_manual_groups()
    group_rel = PdbGroupRelation(groups)

    pdbs = [pdb for g in groups for pdb in g]

    inv_code_map = get_inv_codes_from_file(data_src)
    mat  = DistanceMatrix(mat_id  = "ARGP820101" , code_map = inv_code_map)

    roc_data = perform_roc_test(pdbs, group_rel , mat,roc_step_count = 21)
    #for cutoff,stat in roc_data.items():
        #print "cutoff:%f\tsensitivity:%f\tspecificity:%f" %(cutoff,stat["sensitivity"],stat["specificity"])

    #generate_yard_file(pdbs, group_rel , mat)
