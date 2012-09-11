"""
75 bits version of finger print, in which
the first 15 bits are from the original SIFT algorithm
the next 30 bits describes the spacial features in respective of the antigen
the last 30 bits describes the same thing inrespective of the antiboby
"""

from collections import defaultdict
import numpy as np
from UserList import UserList
from schrodinger import structure, structureutil

import os

from config import *
from util import load_pdb_struct


from fp_gen import distance_tree , distance_data , sift_gen , sift

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
        #print upper_bound,dist
        if dist <= upper_bound:
            return level
    #not in the surrounding
    return -1            

def get_dist_between_res(res1,res2):
    """
    the spatial distance(defined) of two residues
    """
    center1_xyz = np.average([atom.xyz for atom in res1.atom],axis = 0)
    center2_xyz = np.average([atom.xyz for atom in res2.atom],axis = 0)
    diff = np.matrix( center1_xyz - center2_xyz )
    return np.sqrt(( diff * diff.T ).sum())


def get_dist_group_relation(complex1 , complex2):
    """
    get the distance group relations of the residues in the given complexes
    """
    rlt = defaultdict(dict)
    for res1 in complex1.residue:
        for res2 in complex2.residue:
            dist_group = get_dist_group(get_dist_between_res(res1 , res2)) #get the group number id
            rlt[res1][res2] = dist_group 
            rlt[res2][res1] = dist_group 
    return rlt


def get_res_groups_by_dist(complex1 , complex2):
    """
    assign residues to their belonging distance groups,
    return two lists, in which complex1 and complex2 take turns as the host body
    """
    rlt = get_dist_group_relation(complex1 , complex2)
    #print rlt
    host1_list = defaultdict(list)
    host2_list = defaultdict(list)
    for res1 , res2_list in rlt.items():
        for res2 in res2_list:
            group = rlt[res1][res2]
            if group == -1:continue#ignore residues too far away
            host1_list[group].append(res2)
            host2_list[group].append(res1)
    return host1_list , host2_list

class fingerprint_60_bits(UserList):
    """
    the 60 bit finger print class
    """
    def __init__(self):
        UserList.__init__(self,[0] * 60)

    def inc_bit(self, dist_group , inner_bit_ind , is_first_complex):
        if is_first_complex:
            self[dist_group * 6 + inner_bit_ind] += 1#the former 30 bits
        else:         
            self[30 + dist_group * 6 + inner_bit_ind] += 1#the latter 30 bits
        
def get_the_60_bits_fp(complex1,complex2):
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
    fp = fingerprint_60_bits() # 30 bits in this part of finger print
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
            inv_fp_info[item].append(group_index)

    host1_list , host2_list = get_res_groups_by_dist(complex1,complex2)
    #print host1_list
    for dist_group , residues in host1_list.items():
        for res in residues:
            res_code = res.pdbres.strip()
            for inner_pos in inv_fp_info[res_code]:
                fp.inc_bit(dist_group , inner_pos , is_first_complex = True)

    for dist_group , residues in host2_list.items():
        for res in residues:
            res_code = res.pdbres.strip()
            for inner_pos in inv_fp_info[res_code]:
                fp.inc_bit(dist_group , inner_pos , is_first_complex = False)
            
    return fp

#########################
######the 15 bit section
########################

class distance_tree_15(distance_tree):
    
    def __init__(self,  pv_file = None):
        distance_tree.__init__(self , pv_file)
        self.fingerprints = sift_gen_15(self.receptor)

    def find_close_residues(self,  ligand,  cutoff = 4.0):
        lig_name = ligand.title
        d_ = defaultdict(list)
        for atom in ligand.atom:
            #the distance to original point?
            lig_dist = self.measure_dist([atom.xyz])

            #range should refer to the node level in the red black tree
            close_atoms = self.distances.find_nodes_in_range(lig_dist - cutoff, lig_dist + cutoff)

            
            for close in close_atoms:
                atom_atom_dist = self.measure_dist([atom.xyz],  close.data.coords)
                if atom_atom_dist <= cutoff:
                    self.update_min_max(close.data.res_num)
                    d_[close.data.atom].append(atom)

                    #Just take into account the special case, H_DONOR and H_RECEPTOR
                    self.fingerprints.add_sift_chunk_special_case(close.data.atom,  atom,  ligand.title,  atom_atom_dist)
                    
        max_atom_nums = max([len(atom_list) for atom_list in d_.values()])
        #modified definition of interaction
        #print "max in-range atoms count ",max_atom_nums 
        for close_atom, atoms in d_.items():
            #print lig_names
            self.fingerprints.add_sift_chunk_ordinary(close_atom,ligand.title , atoms)#we have only on ligand per complex

class sift_bitset_15(object):
    """ one single sift vector"""
    def __init__(self,  res_num = 0,  bit_set = None):
        
        self.res_num = res_num
        #all bits are of length 2, except position 1,5,6
        self.bit_len_config = [2] * 9
        self.bit_len_config[0] = self.bit_len_config[5] = self.bit_len_config[6] = 1
        self.bit_set = ["0" * length for length in self.bit_len_config]
        print "bit set" , self.bit_set

    def turn_sift_bit_on(self,  sift_bit, strength):
        if strength == "strong":
            self.bit_set[sift_bit] = "11"
        elif strength == "exists":
            if len(self.bit_set[sift_bit]) == 1:#those can never be strong
                self.bit_set[sift_bit] = "1"
            elif len(self.bit_set[sift_bit]) == 2:#those with strong possibilities 
                self.bit_set[sift_bit] = "10"

    def turn_sift_bit_off(self,  sift_bit):
        try:
            self.bit_set[sift_bit] = self.bit_len_config[ sift_bit ] * '0'
        except:
            print "Cannot turn bit " + str(sift_bit) + " off. Too bad..."


class sift_15(sift):
    def __init__(self,  ligand_name , bitset_class):
        sift.__init__(self,  ligand_name , bitset_class)

    def turn_sift_bit_on(self,  sift_bit,  res_num , strength):
        self.sift[res_num].turn_sift_bit_on(sift_bit, strength)
            

class sift_gen_15(sift_gen):
    def __init__(self , receptor):
        sift_gen.__init__(self , receptor)
        self.stren_threshold = {'BACKBONE' : 3,  'SIDECHAIN' : 3,  'POLAR' : 3,  'HYDROPHOBIC' : 3,  'H_ACCEPTOR' : 3,  'H_DONOR' : 3,  'AROMATIC' : 3,  'CHARGED' : 1}

    def get_cur_sift(self , lig_name , rec_atom):
        print "resnum:%4d" %rec_atom.resnum
        if self.sifts.has_key(lig_name):
            if rec_atom.resnum in self.sifts[lig_name].sift.keys():
                cur_sift = self.sifts[lig_name]
            else:
                self.sifts[lig_name].add_sift_chunk(rec_atom.resnum)
                cur_sift = self.sifts[lig_name]
        else:
            new_sift = sift_15(lig_name, bitset_class = sift_bitset_15)
            new_sift.add_sift_chunk(rec_atom.resnum)
            
            self.sifts[lig_name] = new_sift
            cur_sift = new_sift
        return cur_sift                        
        
    def add_sift_chunk_special_case(self,  rec_atom,  lig_atom,  lig_name,  dist):
        """
        the case without strength difference
        """
        cur_sift = self.get_cur_sift(lig_name , rec_atom)

        if structureutil.match_hbond(lig_atom,rec_atom,distance= dist):
            if rec_atom.atomic_number == 1:
                if 'H_DONOR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['H_DONOR'],  rec_atom.resnum, "exists")
            else:
                if 'H_ACCEPTOR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['H_ACCEPTOR'],  rec_atom.resnum, "exists")

    def add_sift_chunk_ordinary(self,  rec_atom,lig_name , atoms):
        """
        the case with strength difference
        """
        cur_sift = self.get_cur_sift(lig_name , rec_atom)
                
        if 'CONTACT' in self.active_bits:
            cur_sift.turn_sift_bit_on(self.bit_pos['CONTACT'],  rec_atom.resnum, "exists")
        if int(rec_atom) in self.backbone_set:
            if 'BACKBONE' in self.active_bits:
                strength = "strong" if len(atoms) >= self.stren_threshold['BACKBONE'] else "exists"
                cur_sift.turn_sift_bit_on(self.bit_pos['BACKBONE'],  rec_atom.resnum, strength)
        else:
            if 'SIDECHAIN' in self.active_bits:
                strength = "strong" if len(atoms) >= self.stren_threshold['SIDECHAIN'] else "exists"
                cur_sift.turn_sift_bit_on(self.bit_pos['SIDECHAIN'],  rec_atom.resnum, strength)
            
            if int(rec_atom) in self.polar_set:
                if 'POLAR' in self.active_bits:
                    strength = "strong" if len(atoms) >= self.stren_threshold['POLAR'] else "exists"
                    cur_sift.turn_sift_bit_on(self.bit_pos['POLAR'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.hydrophobic_set:
                if 'HYDROPHOBIC' in self.active_bits:
                    strength = "strong" if len(atoms) >= self.stren_threshold['HYDROPHOBIC'] else "exists"
                    cur_sift.turn_sift_bit_on(self.bit_pos['HYDROPHOBIC'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.aromatic_set:
                if 'AROMATIC' in self.active_bits:
                    strength = "strong" if len(atoms) >= self.stren_threshold['AROMATIC'] else "exists"
                    cur_sift.turn_sift_bit_on(self.bit_pos['AROMATIC'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.charged_set:
                if 'CHARGED' in self.active_bits:
                    strength = "strong" if len(atoms) >= self.stren_threshold['CHARGED'] else "exists"
                    cur_sift.turn_sift_bit_on(self.bit_pos['CHARGED'],  rec_atom.resnum, strength)

def gen_fp_to_file(receptor=None,binder=None,fp_path=''):
    rec_tree = distance_tree_15()

    rec_tree.set_receptor_structure(receptor)
    rec_tree.parse_receptor()

    rec_tree.find_close_residues(binder , 10.0)

    with open(fp_path,  'w') as out_fp:
        lig_name = rec_tree.fingerprints.sifts.keys()[0]
        for res_num , sift_bitset in rec_tree.fingerprints.sifts[lig_name].sift.items():
            print res_num , sift_bitset
            out_fp.write('%d %s\n' %(res_num , ' '.join([bit for bit in sift_bitset.bit_set])))

    print 'saved to',fp_path

    return fp_path

if __name__ == "__main__":
    comp1 = load_pdb_struct( os.path.join(os.path.dirname( data_src) ,\
                                                     "1DEE_G.pdb") )
    comp2 = load_pdb_struct( os.path.join(os.path.dirname( data_src) ,\
                                                     "1DEE_H.pdb") )

    fp = get_the_60_bits_fp(comp1 , comp2 )
    print fp
    gen_fp_to_file(comp1 , comp2 , "tmp.txt")
