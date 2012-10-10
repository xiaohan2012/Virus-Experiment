"""

Configuration files for data file path, database connection,etc
"""
import os

from machine_setting import base
from dbconfig import *

data_root = os.path.join(base,"data")
fp_root = os.path.join(base,"data/402/fp_result")
pdb_src = os.path.join(base,"data/402/pdb_file/*")
output_dir = os.path.join(base,"data/402/fp_result")
sim_dist_pickle_dir = os.path.join(base,"data/402/sim_dist_pickle")
hydro_yard_dir = os.path.join(base,"data/402/yard_file")
hydro_roc_curve_dir = os.path.join(base,"data/402/roc_curve")
qsub_scripts_path = os.path.join(base,"qsubscripts")
mat_csv_path = os.path.join(base,"data/402/mat_csv")

rmsd_root = os.path.join(data_root,"rmsd")
rmsd_row_path = os.path.join(rmsd_root,"rows.txt")
rmsd_matrix_path = os.path.join(rmsd_root,"matrix.txt")

epi166_root = os.path.join(data_root, "epi_166")
pdb_166_src = os.path.join(epi166_root ,"pdb_file/*")
chain_list_path = os.path.join(epi166_root ,"chains.txt")
dist_mat_path = os.path.join(epi166_root, "dist_mat_csv/dist_mat.csv")
