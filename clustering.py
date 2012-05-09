import Levenshtein
from schrodinger.structure import StructureReader
from schrodinger.structutils import  measure
from collections import OrderedDict,defaultdict
from UserDict import UserDict
import math


class residue_fp(object):
    def __init__(self,fp_str,comp,residue_id):
        self.fp_str = fp_str
        self.complex = comp
        self.residue_id = residue_id
        for res in self.complex.residue:
            if self.residue_id == res.resnum:
                self.res = res
                break

    def get_edit_dist_to(self,residue_fp1):
        return Levenshtein.distance(self.fp_str,residue_fp1.fp_str)

    def within_range_of(self,o_res):
        range_dist = 10
        def distance(xyz1,xyz2):
            #it can be more elegant
            x1,y1,z1 = xyz1
            x2,y2,z2 = xyz2
            return math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        outer_c = 0
        o_res_total = len([r for r in o_res.atom])
        res_total = len([r for r in self.res.atom])
        for a1 in self.res.atom:
            inner_c = 0
            for a2 in o_res.atom:
                if distance(a1.xyz,a2.xyz) < range_dist:
                    inner_c += 1
            if inner_c >= o_res_total / 2.:
                outer_c += 1
        if outer_c >= res_total / 2.:
            print "within range"
            return True
        print "out of range"
        return False
    
    def get_string(self):
        return self.fp_str

class residue_fp_list(UserDict):
    def __init__(self,fp_path = '',comp = None):
        UserDict.__init__(self)
        """construct fp list from file"""
        self.comp = comp
        self.data = OrderedDict()
        with open(fp_path,'r') as f:
            for line in f.readlines()[2:]:
                residue_id = int(line.split()[0])
                fp_str= ''.join(["%d" %float(item) for item in line.split()[1:]])
                self.data[residue_id] = residue_fp(fp_str,self.comp,residue_id)


class dist_mat(UserDict):
    def __init__(self,c1_path = '',c2_path  = '',fp1_path = '',fp2_path = ''):
        self.comp1 = StructureReader(c1_path).next()
        self.comp2 = StructureReader(c2_path).next()
        self.fp1 = residue_fp_list(fp1_path,self.comp1)
        self.fp2 = residue_fp_list(fp2_path,self.comp2)

        self.data = defaultdict(dict)

        for res1, fp1 in self.fp1.items():
            for res2, fp2 in self.fp2.items():
                self.data[res1][res2] = Levenshtein.distance(fp1.get_string(),fp2.get_string())
        self.extended_set = set()
        self.discussed_set = set()

    def find_extending_tuple(self):
        """
        get the minimum distance in distance matrix
        and return the corresponding residue pair
        """
        min_dist = float("inf")
        min_tuple = None
        for res1,_ in self.data.items():
            for res2,dist in _.items():
                if dist < min_dist and (res1,res2,dist) not in self.extended_set:
                    min_dist = dist
                    min_tuple  = (res1,res2,dist)
        return min_tuple                    

    def find_candidate_tuple(self):
        #find residues in range for each complex
        min_dist = float("inf")
        min_tuple = None
        for res1,_ in self.data.items():
            for res2,dist in _.items():
                center1 = self.fp1[self.extending_tuple[0]].res
                center2 = self.fp2[self.extending_tuple[1]].res
                if dist < min_dist and \
                (res1,res2,dist) != self.extending_tuple and (res1,res2,dist) not in self.discussed_set and \
                self.fp1[res1].within_range_of(center1) and self.fp2[res2].within_range_of(center2):
                    min_dist = dist
                    min_tuple = (res1,res2,dist)
        return min_tuple                    
    
    def find_cluster(self):
        self.extending_tuple = self.find_extending_tuple()
        self.extended_set.add(self.extending_tuple)
        print "extending tuple" , self.extending_tuple
        while True:
            cur_tuple = self.find_candidate_tuple()
            if not cur_tuple:
                break
            else:                
                self.discussed_set.add(cur_tuple)
                print "current tuple" , cur_tuple
if __name__ == "__main__":
    data_prefix = "clustering_data"
    c1_path = '%s/1CE1_antibody.pdb' %data_prefix 
    c2_path = '%s/1JTO_antibody.pdb' %data_prefix 
    fp1_path = '%s/fp1.dat' %data_prefix 
    fp2_path = '%s/fp2.dat' %data_prefix 
    pair = dist_mat(c1_path, c2_path ,fp1_path ,fp2_path)
    pair.find_cluster()
    print pair.extended_set,pair.discussed_set

