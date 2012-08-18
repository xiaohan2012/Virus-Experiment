import pymongo

host = "anode4"
dbname = "virus_cluster"

boss_root = "/home/rxzhu/code/Virus-Experiment/"
fp_root = "/home/rxzhu/code/Virus-Experiment/402/fp_result"
data_src = "/home/rxzhu/code/Virus-Experiment/402/pdb_file/*"
output_dir = "/home/rxzhu/code/Virus-Experiment/402/fp_result"
sim_dist_pickle_dir = "/home/rxzhu/code/Virus-Experiment/402/sim_dist_pickle"
hydro_roc_curve_dir = "/home/rxzhu/code/Virus-Experiment/402/roc_curve"

db = pymongo.database.Database(pymongo.Connection(host), dbname)

