import os
import glob
from collections import OrderedDict
from customcollections import OrderedDefaultDict
from numpy import linspace,min,max
from pickle import dump,load
import matplotlib.pyplot as plt

from util import *
from sim_mat import  load_sim_mat
from manual_classification import *
from config import *



class DistanceMatrix(object):
    def __init__(self,mat_id = "dist_mat_402",code_map = {}):
        self.data = load_sim_mat(mat_id)
        self.data = self.data / self.data.diagonal()
        self.cm = code_map

    def get_distance_between(self,pdb1,pdb2):
        ind1 = self.cm[pdb1]
        ind2 = self.cm[pdb2]
        return self.data[ind1][ind2]


def get_roc_data(pdbs, groups , mat,roc_step_count = 11):
    roc_data = OrderedDict()
    for cutoff in linspace(0,1,roc_step_count):
        roc_data[cutoff] = OrderedDefaultdict(int)
        for i in xrange(len(pdbs)):
            for j in xrange(i,len(pdbs)):
                pdb_i = pdbs[i]
                pdb_j = pdbs[j]
                dist = mat.get_distance_between(pdb_i,pdb_j)
                if dist >= cutoff and groups.is_in_one_group(pdb_i,pdb_j) == 1:
                    roc_data[cutoff]["tp"] += 1
                elif dist >= cutoff and groups.is_in_one_group(pdb_i,pdb_j) == -1:                     
                    roc_data[cutoff]["fp"] += 1
                elif dist < cutoff and groups.is_in_one_group(pdb_i,pdb_j) == 1:
                    roc_data[cutoff]["fn"] += 1
                elif dist < cutoff and groups.is_in_one_group(pdb_i,pdb_j) == -1:
                    roc_data[cutoff]["tn"] += 1
                else:    
                    raise RuntimeError("other kind of error")
    cutoff_array = []                
    sensitivity_array = []
    specificity_array = []
    for cutoff,stat in roc_data.items():
        temp = roc_data[cutoff]
        temp["sensitivity"] = temp["tp"] / float(temp["tp"] + temp["fn"])
        temp["specificity"] = temp["fp"] / float(temp["tn"] + temp["fp"])

        sensitivity_array.append(temp["sensitivity"])
        specificity_array.append(temp["specificity"])
        cutoff_array.append(cutoff)

    print sensitivity_array
    print specificity_array
    print cutoff_array
    return roc_data                    

def generate_yard_file(pdbs, groups , mat , yard_file_name = "yard_output.txt"):
    with open(yard_file_name ,'w') as f:
        f.write("output\tmethod1\n")
        for i in xrange(len(pdbs)):
            for j in xrange(i + 1,len(pdbs)):
                pdb_i = pdbs[i]
                pdb_j = pdbs[j]
                f.write("%f\t%f\n" %(1 if groups.is_in_one_group(pdb_i,pdb_j) == 1 else -1,\
                                    mat.get_distance_between(pdb_i,pdb_j)))

def draw_roc():
    y_arr = [1.0, 0.9852941176470589, 0.9411764705882353, 0.9264705882352942, 0.8823529411764706, 0.8823529411764706, 0.8529411764705882, 0.8529411764705882, 0.8235294117647058, 0.8088235294117647, 0.75, 0.6470588235294118, 0.6176470588235294, 0.6029411764705882, 0.4852941176470588, 0.45588235294117646, 0.35294117647058826, 0.2647058823529412, 0.19117647058823528, 0.17647058823529413, 0.0]
    x_arr = [0.9402903527942554, 0.7204183577895723, 0.4120355916328442, 0.19247580393381206, 0.0758663752731814, 0.03332812987823915, 0.018342179206993443, 0.014205432407118327, 0.013034655010927257, 0.012722447705276304, 0.012566344052450827, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089, 0.012488292226038089]
    cutoff = [0.0, 0.050000000000000003, 0.10000000000000001, 0.15000000000000002, 0.20000000000000001, 0.25, 0.30000000000000004, 0.35000000000000003, 0.40000000000000002, 0.45000000000000001, 0.5, 0.55000000000000004, 0.60000000000000009, 0.65000000000000002, 0.70000000000000007, 0.75, 0.80000000000000004, 0.85000000000000009, 0.90000000000000002, 0.95000000000000007, 1.0]
    plt.plot(x_arr,y_arr , "o-")
    for c , x , y in zip(cutoff,x_arr,y_arr):
        print c,x,y
        plt.annotate("%.2f" %c,\
                    xy = (x,y)\
                    )
    plt.show()        
    


if __name__ == "__main__":
    pdb_src = "epi_166/pdb_file/*"

    groups = get_166_manual_groups()
    group_rel = PdbGroupRelation(groups)

    pdbs = []
    for g in groups:
        for i in g:
            pdbs.append(i)
    #print set(wanted).difference(set(pdbs))
    #print len(pdbs),len(set(pdbs))
    inv_code_map = get_inv_codes_from_file(pdb_src)
    mat  = DistanceMatrix(mat_id  = "dist_mat_epi_166.dat" , code_map = inv_code_map)

    roc_data = get_roc_data(pdbs, group_rel , mat,roc_step_count = 21)
    for cutoff,stat in roc_data.items():
        print "cutoff:%f\tsensitivity:%f\tspecificity:%f" %(cutoff,stat["sensitivity"],stat["specificity"])

    generate_yard_file(pdbs, group_rel , mat)
