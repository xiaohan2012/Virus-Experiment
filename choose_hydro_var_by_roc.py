"""
compare different hydrogen contexts by drawing the corresponding roc curves
the roc curve source files are in the format of `yard` file
"""
import os

from draw_roc import generate_yard_file , DistanceMatrix
from sim_mat import *
from config import *
from util.manual_classification import get_166_manual_groups , PdbGroupRelation
from util.aa2code import  get_inv_codes_from_file
from sim_mat import load_sim_mat
from cal_diff_hydro_vars import load_hydro_var

def init_groups_and_pdbs():
    groups = get_166_manual_groups()
    group_rel = PdbGroupRelation(groups)

    pdbs = []
    for g in groups:
        for i in g:
            pdbs.append(i)
    return pdbs , group_rel
                
def plot_roc_curve_for_hydro_vars(mat_ids):
    pdbs , group_rel = init_groups_and_pdbs()
    for mat_id in mat_ids:
        inv_code_map = get_inv_codes_from_file(data_src)
        try:
            mat  = DistanceMatrix(mat_id  = mat_id , code_map = inv_code_map)
        except EOFError:
            print "unexpected"            

        print "roc curve for %s generated" %mat_id
        generate_yard_file(pdbs, group_rel , mat , os.path.join(hydro_yard_dir,"%s.txt" %mat_id))

        yard_path = os.path.join(hydro_yard_dir,"%s.txt" %mat_id)
        gen_img_from_yard_file( yard_path ,\
                                 os.path.join(hydro_roc_curve_dir ,"%s.png" %mat_id) )

def gen_img_from_yard_file(i_path,o_path):
    cmd = "yard-plot --show-auc -o %s %s" %(o_path,i_path)
    print cmd
    os.system(cmd)

def get_auc(i_path):
    cmd = "yard-auc %s" %i_path
    os.system(cmd)

def get_all_auc():
    hydros = map(lambda a:a + "_dist_mat" , load_hydro_var()) #append `_dist_mat` to hydro names
    for mat_id in hydros:
        yard_path = os.path.join(hydro_yard_dir,"%s.txt" %mat_id)
        print mat_id
        get_auc(yard_path)

if __name__ == "__main__":
    hydros = map(lambda a:a + "_dist_mat" , load_hydro_var()) #append `_dist_mat` to hydro names
    plot_roc_curve_for_hydro_vars(hydros)

    #get_all_auc()
