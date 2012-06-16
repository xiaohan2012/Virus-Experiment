import Levenshtein
from schrodinger.structure import StructureReader
from collections import OrderedDict,defaultdict
from UserDict import UserDict
import math

from time_measure import measure_time

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

    #@measure_time
    def within_range_of(self,o_res,range_dist = 20):
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
            return True
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

        self.clustered_fp1_res = set()
        self.clustered_fp2_res = set()
        self.clusters = []

    def find_closest_tuple(self):
        """
        get the minimum distance in distance matrix
        and return the corresponding residue pair
        """
        min_dist = float("inf")
        min_tuple = None
        for res1,_ in self.data.items():
            for res2,dist in _.items():
                if dist < min_dist and \
                   res1 not in self.clustered_fp1_res and \
                   res2 not in self.clustered_fp2_res and \
                   (res1,res2,dist) not in self.not_suitable_tuple:
                    min_dist = dist
                    min_tuple  = (res1,res2,dist)
        return min_tuple                    

    def find_cluster(self):
        cluster = set()
        self.not_suitable_tuple = set()

        self.extending_tuple = self.find_closest_tuple()
        if self.extending_tuple is None:
            return
        ext_res1 , ext_res2 , dist = self.extending_tuple
        cluster.add(self.extending_tuple)
        self.clustered_fp1_res.add(ext_res1)
        self.clustered_fp2_res.add(ext_res2)

        #print "extending tuple" , self.extending_tuple

        while True:
            t = self.find_closest_tuple()#get the next closest(edit distance) tuple
            if t is None:#cannot find any appropriate tuple
                #print "cluster full, step out"
                break
            center_res1 = self.fp1[ext_res1].res
            center_res2 = self.fp2[ext_res2].res
            res1 , res2 , dist = t
            if self.fp1[res1].within_range_of(center_res1) and self.fp2[res2].within_range_of(center_res2):
                #check if it is within the range of extending tuple,if so, add it to the cluster
                cluster.add(t)

                self.clustered_fp1_res.add(res1)
                self.clustered_fp2_res.add(res2)
                #print "cluster size: %d,total residue number: f1 = %d, f2 = %d" %(len(cluster),len(self.fp1),len(self.fp2))
            else:
                ##print "out of range"
                self.not_suitable_tuple.add(t)
                            
        self.clusters.append(cluster)
        self.extending_tuple = None
        return cluster

res_sim_mat = {"AA" : 1,"AC" : 0,"AD" : 0,"AE" : 0,"AF" : 0,"AG" : 0,"AH" : 0,"AI" : 0,"AK" : 0,"AL" : 0,"AM" : 0,"AN" : -2,"AP" : -1,"AQ" : -1,"AR" : -1,"AS" : 1,"AT" : 0,"AV" : 0,"AW" : -3,"AY" : -2,"CC" : 9,"CD" : -3,"CE" : -4,"CF" : -2,"CG" : -3,"CH" : -3,"CI" : -1,"CK" : -3,"CL" : -1,"CM" : -1,"CN" : -3,"CP" : -3,"CQ" : -3,"CR" : -3,"CS" : -1,"CT" : -1,"CV" : -1,"CW" : -2,"CY" : -2,"DD" : 6,"DE" : 2,"DF" : -3,"DG" : -1,"DH" : -1,"DI" : -3,"DK" : -1,"DL" : -4,"DM" : -3,"DN" : 1,"DP" : -1,"DQ" : 0,"DR" : -2,"DS" : 0,"DT" : -1,"DV" : -3,"DW" : -4,"DY" : -3,"EE" : 5,"EF" : -3,"EG" : -2,"EH" : 0,"EI" : -3,"EK" : 1,"EL" : -3,"EM" : -2,"EN" : 0,"EP" : -1,"EQ" : 2,"ER" : 0,"ES" : 0,"ET" : -1,"EV" : -2,"EW" : -3,"EY" : -2,"FF" : 6,"FG" : -3,"FH" : -1,"FI" : 0,"FK" : -3,"FL" : 0,"FM" : 0,"FN" : -3,"FP" : -4,"FQ" : -3,"FR" : -3,"FS" : -2,"FT" : -2,"FV" : -1,"FW" : 1,"FY" : 3,"GG" : 6,"GH" : -2,"GI" : -4,"GK" : -2,"GL" : -4,"GM" : -3,"GN" : 0,"GP" : -2,"GQ" : -2,"GR" : -2,"GS" : 0,"GT" : -2,"GV" : -3,"GW" : -2,"GY" : -3,"HH" : 8,"HI" : -3,"HK" : -1,"HL" : -3,"HM" : -2,"HN" : 1,"HP" : -2,"HQ" : 0,"HR" : 0,"HS" : -1,"HT" : -2,"HV" : -3,"HW" : -2,"HY" : 2,"II" : 4,"IK" : -3,"IL" : 2,"IM" : 1,"IN" : -3,"IP" : -3,"IQ" : -3,"IR" : -3,"IS" : -2,"IT" : -1,"IV" : 3,"IW" : -3,"IY" : -1,"KK" : 5,"KL" : -2,"KM" : -1,"KN" : 0,"KP" : -1,"KQ" : 1,"KR" : 2,"KS" : 0,"KT" : -1,"KV" : -2,"KW" : -3,"KY" : -2,"LL" : 4,"LM" : 2,"LN" : -3,"LP" : -3,"LQ" : -2,"LR" : -2,"LS" : -2,"LT" : -1,"LV" : 1,"LW" : -2,"LY" : -1,"MM" : 5,"MN" : -2,"MP" : -2,"MQ" : 0,"MR" : -1,"MS" : -1,"MT" : -1,"MV" : 1,"MW" : -1,"MY" : -1,"NN" : 6,"NP" : -2,"NQ" : 0,"NR" : 0,"NS" : 1,"NT" : 0,"NV" : -3,"NW" : -4,"NY" : -2,"PP" : 7,"PQ" : -1,"PR" : -2,"PS" : -1,"PT" : -1,"PV" : -2,"PW" : -4,"PY" : -3,"QQ" : 5,"QR" : 1,"QS" : 0,"QT" : -1,"QV" : -2,"QW" : -2,"QY" : -1,"RR" : 5,"RS" : -1,"RT" : -1,"RV" : -3,"RW" : -3,"RY" : -2,"SS" : 4,"ST" : 1,"SV" : -2,"SW" : -3,"SY" : -2,"TT" : 5,"TV" : 0,"TW" : -2,"TY" : -2,"VV" : 4,"VW" : -3,"VY" : -1,"WW" : 11,"WY" : 2,"YY" : 7}

def get_similarity(clusters,pair):
    val1,val2,val3 = 0,0,0

    #calculate value1
    for c in clusters:
        for t in c:
            val1 -= t[2]#the edit distance
    #calculate value2
    for c in clusters:
        for t in c:
            res1,res2,dist = t
            res1_code = pair.fp1[res1].res.getCode()
            res2_code = pair.fp2[res2].res.getCode()
            try:
                val2 += res_sim_mat[res1_code + res2_code]
            except KeyError:
                val2 += res_sim_mat[res2_code + res1_code ]

    #calculate value3
    for c in clusters:
        val3 += len(c)#pair count
    #print "value1:%d,value2:%d,value3:%d" %(val1,val2,val3)
    return val1+val2+val3


    #data_prefix = "clustering_data"
    #c1_path = '%s/1MEL&&1MLC/receptor1.pdb' %data_prefix 
    #c2_path = '%s/1MEL&&1MLC/receptor2.pdb' %data_prefix 
    #fp1_path = '%s/1MEL&&1MLC/avg_sift1.out' %data_prefix 
    #fp2_path = '%s/1MEL&&1MLC/avg_sift2.out' %data_prefix 

    #pair = dist_mat(c1_path, c2_path ,fp1_path ,fp2_path)
    ##print pair.find_cluster()
    #while pair.find_cluster():
        #pass
    ##print pair.clusters
    #get_similarity(pair.clusters,pair)
