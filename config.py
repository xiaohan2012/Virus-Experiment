"""

Configuration files for data file path, database connection,etc
"""
import pymongo

host = "anode4"
dbname = "virus_cluster"

boss_root = "/home/rxzhu/code/Virus-Experiment/data"
fp_root = "/home/rxzhu/code/Virus-Experiment/data/402/fp_result"
data_src = "/home/rxzhu/code/Virus-Experiment/data/402/pdb_file/*"
output_dir = "/home/rxzhu/code/Virus-Experiment/data/402/fp_result"
sim_dist_pickle_dir = "/home/rxzhu/code/Virus-Experiment/data/402/sim_dist_pickle"
hydro_yard_dir = "/home/rxzhu/code/Virus-Experiment/data/402/yard_file"
hydro_roc_curve_dir = "/home/rxzhu/code/Virus-Experiment/data/402/roc_curve"
qsub_scripts_path = "/home/rxzhu/code/Virus-Experiment/qsubscripts"

db = pymongo.database.Database(pymongo.Connection(host), dbname)

