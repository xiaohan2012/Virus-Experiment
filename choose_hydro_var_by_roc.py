"""
compare different hydrogen contexts by drawing the corresponding roc curves
the roc curve source files are in the format of `yard` file
"""
import os , re, sys

from roc import generate_yard_file , print_auc , perform_roc_test
from dist_mat import DistanceMatrix
from sim_mat import *
from config import *
from util.manual_classification import get_166_manual_groups , PdbGroupRelation
from sim_mat import load_sim_mat
from cal_diff_hydro_vars import load_hydro_var

def init_pdb_groups_relation():
    groups = get_166_manual_groups()
    group_rel = PdbGroupRelation(groups)

    return group_rel

####using yard                
def plot_roc_for_hydro_vars(mat_ids):
    group_rel = init_pdb_groups_relation()
    for mat_id in mat_ids:
        mat  = DistanceMatrix(mat_id  = mat_id )

        print "roc curve for %s generated" %mat_id
        generate_yard_file(group_rel , mat , os.path.join(hydro_yard_dir,"%s.txt" %mat_id))

        yard_path = os.path.join(hydro_yard_dir,"%s.txt" %mat_id)
        gen_img_from_yard_file( yard_path ,\
                                 os.path.join(hydro_roc_curve_dir ,"%s.png" %mat_id) )

def gen_img_from_yard_file(i_path,o_path):
    cmd = "yard-plot --show-auc -o %s %s" %(o_path,i_path)
    print cmd
    os.system(cmd)


def get_all_auc():
    hydros = map(lambda a:a + "_dist_mat" , load_hydro_var()) #append `_dist_mat` to hydro names
    for mat_id in hydros:
        yard_path = os.path.join(hydro_yard_dir,"%s.txt" %mat_id)
        print mat_id
        print_auc(yard_path)

####using matplotlib to plot###
import matplotlib.pyplot as plt

def batch_plot(mat_ids , test_iter_counts = 22):
    """ plot the given ids in one plot"""

    group_rel = init_pdb_groups_relation()
    plt.hold(True)
    color_scheme = []
    for mat_id in mat_ids:
        mat  = DistanceMatrix(mat_id  = mat_id )
        x , y , _ = perform_roc_test(group_rel  , mat , test_iter_counts )
        plt.plot(x , y )
    plt.savefig("tmp.jpg")

def print_js_style_roc_test_result(mat_id):
    group_rel = init_pdb_groups_relation()
    mat  = DistanceMatrix(mat_id  = mat_id )
    x_arr , y_arr , _ = perform_roc_test(group_rel  , mat  ,101)
    xy_arr = []
    x_arr.reverse();y_arr.reverse();#reverse the order so that jqplot can draw properly
    for x,y in zip(x_arr , y_arr):
        xy_arr.append("[%f , %f]" %(x , y))
    sys.stderr.write("var $%s = [ %s ];\n" %(re.findall(r"(\w+)_dist_mat" , mat_id )[0], ','.join(xy_arr)))        

def batch_plotting(hydros):
    hydros = map(lambda a:a + "_dist_mat" , hydros) #append `_dist_mat` to hydro names
    for mat_id in hydros:
        #var_name = re.findall(r"(\w+)_dist_mat" , mat_id )[0]
        #sys.stderr.write("$%s ," %var_name );
        print_js_style_roc_test_result(mat_id);
    
    #write other js snippets
    sys.stderr.write("var hydro_names = [%s];\n" %(",".join(["\"%s\"" %hydro[:-9] for hydro in hydros])))
    sys.stderr.write("var data = [%s];\n" %(",".join(["$%s" %hydro[:-9] for hydro in hydros])))
    if len(hydros) <= 6:
        sys.stderr.write("var groups = [{data:data,\nnames:hydro_names}];\n")
    else:
        sys.stderr.write("""
var groups = []; 
//divide the data into 6 groups, each containing 6 cases.
for(var i = 0; i < 6 ;i++){
    groups.push({
    data:data.slice(i * 6 , i * 6 + 6), 
    names:hydro_names.slice(i * 6 , i * 6 + 6)
    }); 
}
""")
    
if __name__ == "__main__":
    #load_hydro_var()

    batch_plotting(["ARGP820101"])
    #plot_roc_for_hydro_vars(hydros)
    #batch_plot(hydros[:2])

    #get_all_auc()
