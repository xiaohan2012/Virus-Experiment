from draw_roc import generate_yard_file , DistanceMatrix
from sim_mat import *
from config import *
from manual_classification import get_166_manual_groups , PdbGroupRelation
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
        generate_yard_file(pdbs, group_rel , mat , os.path.join(hydro_roc_curve_dir,"%s.txt" %mat_id))

if __name__ == "__main__":
    hydros = load_hydro_var()
    plot_roc_curve_for_hydro_vars(hydros)
