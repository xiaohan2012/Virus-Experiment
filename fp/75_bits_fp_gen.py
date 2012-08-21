"""
75 bits version of finger print, in which
the first 15 bits are from the original SIFT algorithm
the next 30 bits describes the spacial features in respective of the antigen
the last 30 bits describes the same thing inrespective of the antiboby
"""

from collections import defaultdict
import numpy as np
from UserList import UserList

def get_dist_group(dist , bound_list = [4. , 8. , 12. , 16. , 20.]):
    """
    get the group index of the distance it should belong to
    for example:
        distance groups are :
            0 -> 4 , 4 -> 8 ,  8 -> 12
        then disttance like :
            1.2 belong to dist group 1 (0 -> 4)
            4.8 belong to dist group 2 (4 -> 8)
            etc
    """
    for level,upper_bound in enumerate(bound_list):
        if dist <= upper_bound:
            return level
    #not in the surrounding
    return -1            

def get_dist_between_res(res1,res2):
    """
    the spatial distance(defined) of two residues
    """
    raise NotImplementedError("Please, oh, no")

def get_dist_group_relation(complex1 , complexa2):
    """
    get the distance group relations of the residues in the given complexes
    """
    rlt = defaultdict(dict)
    for res1 in complex1.residue:
        for res2 in complex2.residue:
            dist_group = get_dist_group(get_dist_between_res(res1 , res2))
            rlt[res1][res2] = dist_group 
            rlt[res2][res1] = dist_group 
    return rlt


def get_res_groups_by_dist(complex1 , complex2):
    """
    assign residues to their belonging distance groups,
    return two lists, in which complex1 and complex2 take turns as the host body
    """
    rlt = get_dist_group_relation(complex1 , complex2)
    host1_list = [] , host2_list = []
    for res1 , res2_list in rlt.items():
        for res2 in res2_list:
            group = rlt[res1][res2]
            host1_list[group].append(res2)
            host2_list[group].append(res1)
    return host1_list , host2_list

class fingerprint_60_bits(UserList):
    """
    the 60 bit finger print class
    """
    def __init__(self,length = 60):
        UserList.__init__(self,[0] * length)

    def inc_bit(self, dist_group , inner_bit_ind , is_first_complex):
        if first_complex:
            self[dist_group * 6 + inner_bit_ind] += 1#the former 30 bits
        else:         
            self[len(self) / 2 + dist_group * 6 + inner_bit_ind] += 1#the latter 30 bits
        
def gen_the_60_bits_fp(complex1,complex2):
    """
    physical properties:

    Polar :TYR,ASN,GLU,SER,CYS,THR,GLY
    Hydrop: PHE, LEU, ILE, TRP, VAL, MET, PRO, ALA
    Charged: ARG, ASP, GLU, LYS, HIS

    chemical properties:

    Lipids: ALA, VAL, LEU, ILE, MET, ASN, GLU, LYS, ARG, GLY, SER, THR, CYS, ASP, PHE
    Aromatic: PHE, TYR, TRP
    Heterocyclic: PRO, HIS
    """
    fp = fingerprint_60_bits(60) # 30 bits in this part of finger print
    fp_info = [
     ['TYR', 'ASN', 'GLU', 'SER', 'CYS', 'THR', 'GLY'],#polar 
     ['PHE', 'LEU', 'ILE', 'TRP', 'VAL', 'MET', 'PRO', 'ALA'],#hydrop
     ['ARG', 'ASP', 'GLU', 'LYS', 'HIS'],#charged
     ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'ASN', 'GLU', 'LYS', 'ARG', 'GLY', 'SER', 'THR', 'CYS', 'ASP', 'PHE'],#lipids
     ['PHE', 'TYR', 'TRP'],#aromatic
     ['PRO','HIS'],#heterocyclic
    ]
    #convert it to inverted index
    inv_fp_info = defaultdict(list)
    for group_index , group_items in enumerate(fp_info):
        for item in group_items:
            inv_fp_info[item].append(item)

    host1_list , host2_list = get_res_groups_by_dist(complex1,complex2)
    for dist_group , residues in enumerate(host1_list):
        for res in residues:
            res_code = res.pdbres.strip()
            for inner_pos in inv_fp_info[res_code]:
                fp.inc_bit(dist_group , inner_pos , is_first_complex = True)

    for dist_group , residues in enumerate(host2_list):
        for res in residues:
            res_code = res.pdbres.strip()
            for inner_pos in inv_fp_info[res_code]:
                fp.inc_bit(dist_group , inner_pos , is_first_complex = False)
            
    return fp

    



