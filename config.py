"""

Configuration files for data file path, database connection,etc
"""
import pymongo
import os


from machine_setting import base

boss_root = os.path.join(base,"data")
fp_root = os.path.join(base,"data/402/fp_result")
data_src = os.path.join(base,"data/402/pdb_file/*")
output_dir = os.path.join(base,"data/402/fp_result")
sim_dist_pickle_dir = os.path.join(base,"data/402/sim_dist_pickle")
hydro_yard_dir = os.path.join(base,"data/402/yard_file")
hydro_roc_curve_dir = os.path.join(base,"data/402/roc_curve")
qsub_scripts_path = os.path.join(base,"qsubscripts")
mat_csv_path = os.path.join(base,"data/402/mat_csv")

