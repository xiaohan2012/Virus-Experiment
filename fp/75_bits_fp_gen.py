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
from customcollections import OrderedDefaultDict

import os

from config import *
from util.load_pdb import load_pdb_struct


from fp_gen import distance_tree , distance_data , sift_gen , sift

class FingerPrint_60(OrderedDefaultDict):
    def __init__(self , antigen , antibody):
        OrderedDefaultDict.__init__(self , lambda :list([0] * 60))

        self.antigen = antigen
        self.antibody = antibody
        self.nearby_relation = defaultdict(dict)#for cache

        self.nearby_reses_in_antigen = OrderedDefaultDict(lambda : defaultdict(list))#surrouding residues in antigen for each residue 
        self.nearby_reses_in_antibody = OrderedDefaultDict(lambda : defaultdict(list))#surrouding residues in antibody for each residue

        self.fp_rule = {#property id and the corresponding residue code
             0 : ['TYR', 'ASN', 'GLU', 'SER', 'CYS', 'THR', 'GLY'],         #polar 
             1 : ['PHE', 'LEU', 'ILE', 'TRP', 'VAL', 'MET', 'PRO', 'ALA'],  #hydrop
             2 : ['ARG', 'ASP', 'GLU', 'LYS', 'HIS'],                       #charged
             3 : ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'ASN', 'GLU', 'LYS',\
              'ARG', 'GLY', 'SER', 'THR', 'CYS', 'ASP', 'PHE'],             #lipids
             4 : ['PHE', 'TYR', 'TRP'],                                     #aromatic
             5 : ['PRO','HIS'],                                             #heterocyclic
        }#the key represents the group index, value for the residue code

        self.res_prop_ids = defaultdict(list)#the property ids that a given residue has
        #we need to do some conversion for fp_rule for better performance

        for prop_id , residues in self.fp_rule.items():
            for res_code in residues:
                self.res_prop_ids[res_code].append(prop_id)
        #print "res_prop_ids",self.res_prop_ids                

        self.atom_dist_cutoff = 4.0

    def residue_nearby_enough(self,res1 , res2):
        """
        determine whether two atoms are nearby enough given the `atom_dist_cutoff`
        """
        def atom_distance(atom1 , atom2):#the distance between two atoms
            diff = np.matrix( np.array(atom1.xyz) - np.array(atom2.xyz))
            return np.sqrt(( diff * diff.T ).sum())

        if self.nearby_relation[res1].has_key(res2):#if it has been computed
            return  self.nearby_relation[res1][res2]

        for atom1 in res1.atom:
            for atom2 in res2.atom:
                if atom_distance(atom1 , atom2) <= self.atom_dist_cutoff:
                    self.nearby_relation[res1][res2] = True#cache the result
                    self.nearby_relation[res2][res1] = True#the symetrical case
                    return True
        self.nearby_relation[res1][res2] = False#cache the result
        self.nearby_relation[res2][res1] = False#the symetrical case
        return False                

    def grouping_residue_by_distance(self):
        """iterate every residue in the complex and group their surrounding residues by distance"""
        def res_distance(res1,res2):
            """residues distance """
            diff = np.matrix( np.average([atom.xyz for atom in res1.atom],axis = 0) - \
                              np.average([atom.xyz for atom in res2.atom],axis = 0) )
            return np.sqrt(( diff * diff.T ).sum())

        def get_dist_group(dist , bound_list = [4. , 8. , 12. , 16. , 20.]):
            """get the group index it should belong to according to the distance """
            for level,upper_bound in enumerate(bound_list):
                #print upper_bound,dist
                if dist <= upper_bound:
                    return level
            #not in the surrounding
            return -1            
        
        #grouping the residues in antigen
        for res1 in self.antigen.residue:
            for res2 in self.antigen.residue:
                if res1 is res2:continue
                if self.residue_nearby_enough(res1 , res2):
                    dist = res_distance(res1 , res2 )#get the distance between res1 and res2
                    dist_group = get_dist_group(dist)#fit it into a group 
                    self.nearby_reses_in_antigen[res1][dist_group].append(res2)#updating the group list

        #grouping the residues in antibody
        for res1 in self.antigen.residue:
            for res2 in self.antibody.residue:
                if res1 is res2:continue
                if self.residue_nearby_enough(res1 , res2):
                    dist = res_distance(res1 , res2 )#get the distance between res1 and res2
                    dist_group = get_dist_group(dist)#fit it into a group 
                    self.nearby_reses_in_antibody[res1][dist_group].append(res2)#updating the group list

        #print "nearby_reses_in_antigen" , self.nearby_reses_in_antigen

    def get_fingerprint(self):
        if not self:#not computed
            self.grouping_residue_by_distance()#first group those residues
            
            #the first 30 bits
            for res , groups in self.nearby_reses_in_antigen.items():
                for group_index , residues in groups.items():
                    for residue in residues:
                        for prop_id in self.res_prop_ids[residue.pdbres.strip().upper()]:
                            #increment the count of property at given position
                            self[res][group_index * 6 + prop_id] += 1
                            #print "#####for residue %d" %res.resnum
                            #print "%s" %(" ".join("%dp%d" %(g,p) for g in xrange(5) for p in xrange(6)))
                            #print ' '.join("%2d " %count for count in self[res])
                            #print res.resnum , group_index , prop_id
            #the 30 ~ 60 bits                            
            for res , groups in self.nearby_reses_in_antibody.items():
                for group_index , residues in groups.items():
                    for residue in residues:
                        for prop_id in self.res_prop_ids[residue.pdbres.strip().upper()]:
                            #increment the count of property at given position, offset by 30
                            self[res][30 + group_index * 6 + prop_id] += 1
            return self

    def display_fingerprint(self , start = None , end = None):
        print "%s%s" %(' ' * 11 , " ".join("%dp%d" %(g,p) for g in xrange(5) for p in xrange(6)))
        for residue , fp in self.items():
            if start and end:
                fp = fp[start:end]
            elif not start and end:
                fp = fp[:end]
            elif start and not end:
                fp = fp[start:]

            print "%8d : %s" %(residue.resnum , ' '.join("%2d " %count for count in fp))


        
    def display_group_info(self):
        def _display_group_info(nearby_reses):
            for residue , groups in nearby_reses :
                print "%8d" %residue.resnum , 
                for group_index , residues in groups.items():
                    #twisted statement, hehe!
                    print "%d:%d(%s)" %( group_index , len(residues) ,\
                             ' '.join("%s(%s)" %(res.pdbres.strip().upper(),\
                                                ','.join('%d'%prop_id  for prop_id in self.res_prop_ids[res.pdbres.strip().upper()]))\
                                                     for res in  residues)),
                print                
            return#for clearity

        print "antigen part(first 30 bit)"
        _display_group_info(self.nearby_reses_in_antigen.items())
        print "antibody part(30 ~ 60 bit)"
        _display_group_info(self.nearby_reses_in_antibody.items())

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
    antigen = load_pdb_struct( os.path.join(os.path.dirname( data_src) ,\
                                                     "1DEE_G.pdb") )
    antibody  = load_pdb_struct( os.path.join(os.path.dirname( data_src) ,\
                                                     "1DEE_H.pdb") )
    fp = FingerPrint_60(antigen , antibody)

    fp.get_fingerprint()

    #print ' '.join("%d" %res.resnum for res in antigen.residue)#check it now

    fp.display_fingerprint(start = 30)
    fp.display_group_info()

    #gen_fp_to_file(antigen , antibody , "tmp.txt")
