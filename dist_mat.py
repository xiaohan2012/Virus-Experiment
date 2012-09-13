from sim_mat import  load_sim_mat
from util.aa2code import  get_inv_codes_from_file
from config import *

class DistanceMatrix(object):
    def __init__(self,mat_id = "dist_mat_402" , code_map = {}):
        self.data = load_sim_mat(mat_id)
        self.data = self.data / self.data.diagonal()
        if code_map:
            self.cm = code_map
        else:            
            self.cm = get_inv_codes_from_file(data_src)

    def get_distance_between(self,pdb1,pdb2):
        ind1 = self.cm[pdb1]
        ind2 = self.cm[pdb2]
        return self.data[ind1][ind2]
