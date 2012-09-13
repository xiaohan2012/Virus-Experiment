#!/usr/bin/env python

import os,  sys
from schrodinger import structure, structureutil
import schrodinger.utils.fileutils as fileutils
import numpy as np
from collections import defaultdict

from fp_gen import distance_tree , distance_data , sift , gen_fp_to_file


BLACK = 0
RED = 1

class distance_tree_22(distance_tree):
    
    def __init__(self,  pv_file = None):
        distance_tree.__init__(self , pv_file)

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
                    d_[close.data.atom].append((ligand.title , atom))

                    #Just take into account the special case, H_DONOR and H_RECEPTOR
                    self.fingerprints.add_sift_chunk_special_case(close.data.atom,  atom,  ligand.title,  atom_atom_dist)
                    
        max_atom_nums = max([len(atom_list) for atom_list in d_.values()])
        #modified definition of interaction
        #print "max in-range atoms count ",max_atom_nums 
        strong_c = 0;
        exists_c = 0;
        for close_atom, lig_names in d_.items():
            #print lig_names
            if len(lig_names) >= 3:
                #print "mark as strong"
                self.fingerprints.add_sift_chunk_ordinary(close_atom,lig_names[0][0] ,"strong")#we have only on ligand per complex
                strong_c += 1
            else:
                #print "mark as exists"
                self.fingerprints.add_sift_chunk_ordinary(close_atom,lig_names[0][0] ,"exists")
                exists_c += 1
        print sorted([len(lig_names) for lig_names in d_.values()])
        
        ####### start of environmental descriptors ########
        #print "strong count",strong_c ,"exists count",exists_c
        #addition of surrounding environment                
        POLAR_RESIDUES = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","ARN","ASH","GLH","HID","HIE","LYN"]
        HYDROPHOBIC_RESIDUES =["PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","CYX"]
        AROMATIC_RESIDUES = ["PHE","TYR","TRP","TYO"]
        CHARGED_RESIDUES = ["ARG","ASP","GLU","LYS","HIP","CYT","SRO","TYO","THO"]

        polar_list = structureutil.evaluate_asl( ligand,
                                                    "res. " + ', '.join(POLAR_RESIDUES ))
        polar_set = set( polar_list )
        hydrophobic_list = structureutil.evaluate_asl( ligand,
                                                    "res. " + ', '.join(HYDROPHOBIC_RESIDUES))
        hydrophobic_set = set( hydrophobic_list )
        aromatic_list = structureutil.evaluate_asl( ligand,
                                                    "res. " + ', '.join(AROMATIC_RESIDUES))
        aromatic_set = set(aromatic_list)
        charged_list = structureutil.evaluate_asl ( ligand,
                                                    "res. " + ', '.join(CHARGED_RESIDUES))
        charged_set = set( charged_list )

        #print polar_set, hydrophobic_set ,aromatic_set , charged_set 
        env_d_ = dict()
        for close_atom, lig_names in d_.items():
            for lig_atom in [l[1] for l in lig_names]:
                lig_atom = int(lig_atom)
                if close_atom not in env_d_.keys():
                    env_d_[close_atom] = defaultdict(int)
                if lig_atom in polar_set:
                    #print "it is polar"
                    env_d_[close_atom]['ENV_POLAR'] += 1
                elif lig_atom in aromatic_set:
                    #print "it is aromatic"
                    env_d_[close_atom]['ENV_AROMATIC'] += 1
                elif lig_atom in hydrophobic_set:
                    #print "it is h"
                    env_d_[close_atom]['ENV_HYDROPHOBIC'] += 1
                elif lig_atom in charged_set:
                    #print "it is c"
                    env_d_[close_atom]['ENV_CHARGED'] += 1
                    
        atom_nums = list()
        for a,fp in env_d_.items():
            atom_nums += fp.values()
        max_atom_nums = max(atom_nums)
        env_stat = {"p_list":[0],"a_list":[0],"h_list":[0],"c_list":[0]}
        env_stat["p_list"] += [fp["ENV_POLAR"] for a,fp in env_d_.items() if fp["ENV_POLAR"] != 0]
        env_stat["a_list"] += [fp["ENV_AROMATIC"] for a,fp in env_d_.items() if fp["ENV_AROMATIC"] != 0]
        env_stat["h_list"] += [fp["ENV_HYDROPHOBIC"] for a,fp in env_d_.items() if fp["ENV_HYDROPHOBIC"] != 0]
        env_stat["c_list"] += [fp["ENV_CHARGED"] for a,fp in env_d_.items() if fp["ENV_CHARGED"] != 0]
        print env_stat
        print max(env_stat["p_list"]),max(env_stat["a_list"]),max(env_stat["h_list"]),max(env_stat["c_list"])
        strong_c = 0;
        exists_c = 0;
        for a,fp in env_d_.items():
            for t,c in fp.items():
                if c >= 3:
                    self.fingerprints.add_sift_chunk_env(a,lig_name,t,"strong")
                    strong_c += 1
                else:                    
                    self.fingerprints.add_sift_chunk_env(a,lig_name,t,"exists")
                    exists_c += 1

class sift_bitset(object):
    """ one single sift vector"""
    def __init__(self,  res_num = 0,  bit_set = None):
        
        self.res_num = res_num
        
        if bit_set == None:
            self.bit_set = ["00" for x in range(13)]
        else:
            self.set_bitset(bit_set)


    def set_bitset(self,  sift_string):
        
        if len(sift_string) == 13:
            self.bit_set = sift_string
        else:
            print "set_bitset(): invalid sift length, skipping..."

    def turn_sift_bit_on(self,  sift_bit, strength):
        
        try:
            if strength == "strong":
                self.bit_set[sift_bit] = "11"
            elif strength == "exists":
                self.bit_set[sift_bit] = "10"
        except:
            print "Cannot turn bit " + str(sift_bit) + " on. Too bad..."

    def turn_sift_bit_off(self,  sift_bit):
        try:
            self.bit_set[sift_bit] = "00"
        except:
            print "Cannot turn bit " + str(sift_bit) + " off. Too bad..."

class sift_gen(object):
    
    def __init__(self,  receptor):
        
        self.POLAR_RESIDUES = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","ARN","ASH","GLH","HID","HIE","LYN"]
        self.HYDROPHOBIC_RESIDUES =["PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","CYX"]
        self.AROMATIC_RESIDUES = ["PHE","TYR","TRP","TYO"]
        self.CHARGED_RESIDUES = ["ARG","ASP","GLU","LYS","HIP","CYT","SRO","TYO","THO"]

        # Bit positions for each interaction:
        self.bit_pos = {'CONTACT' : 0,  'BACKBONE' : 1,  'SIDECHAIN' : 2,  'POLAR' : 3,  'HYDROPHOBIC' : 4,  'H_ACCEPTOR' : 5,  'H_DONOR' : 6,  'AROMATIC' : 7,  'CHARGED' : 8, 'ENV_POLAR' : 9,  'ENV_HYDROPHOBIC' : 10,'ENV_AROMATIC' : 11,  'ENV_CHARGED' : 12}

        self.active_bits = ['CONTACT',  'BACKBONE',  'SIDECHAIN',  'POLAR', 'HYDROPHOBIC',  'H_ACCEPTOR',  'H_DONOR',  'AROMATIC',  'CHARGED','ENV_POLAR' ,  'ENV_HYDROPHOBIC' ,'ENV_AROMATIC',  'ENV_CHARGED']
        self.max_bits_per_residue =  13 
        self.used_bits_number = len(self.active_bits)
        
        self.sifts = {}

        self.receptor_st = receptor

        backbone_list = structureutil.evaluate_asl( self.receptor_st,
                                                    "backbone")
        self.backbone_set = set(backbone_list)

        """
        structureutil.evaluate_asl
        
        Search for substructures matching the ASL (Atom Specification Language)

        the atom index is returned
        """
        
        polar_list = structureutil.evaluate_asl( self.receptor_st,
                                                    "res. " + ', '.join(self.POLAR_RESIDUES ))
        self.polar_set = set( polar_list )
        hydrophobic_list = structureutil.evaluate_asl( self.receptor_st,
                                                    "res. " + ', '.join(self.HYDROPHOBIC_RESIDUES))
        self.hydrophobic_set = set( hydrophobic_list )
        aromatic_list = structureutil.evaluate_asl( self.receptor_st,
                                                    "res. " + ', '.join(self.AROMATIC_RESIDUES))
        self.aromatic_set = set(aromatic_list)
        charged_list = structureutil.evaluate_asl ( self.receptor_st,
                                                    "res. " + ', '.join(self.CHARGED_RESIDUES))
        self.charged_set = set( charged_list )

    def set_active_bits(self,  bit_list):
    
        self.active_bits = bit_list
    
    
    def get_active_bit_list(self):
        
        return self.active_bits
    
    
    def print_active_bit_list(self):
        
        for bit in self.active_bits:
            print bit + "\t"
        print


    def fill_missing_zeros(self, start_res,  end_res,  lig_name):
        if lig_name not in self.sifts.keys():
            new_sift = sift(lig_name)
        
        for res in range(start_res,  end_res):
            if res not in self.sifts[lig_name].sift.keys():
                self.sifts[lig_name].add_sift_chunk(res)
                
    def add_sift_chunk_special_case(self,  rec_atom,  lig_atom,  lig_name,  dist):
        """
        the case without strength difference
        """
        if self.sifts.has_key(lig_name):
            if rec_atom.resnum in self.sifts[lig_name].sift.keys():
                cur_sift = self.sifts[lig_name]
            else:
                self.sifts[lig_name].add_sift_chunk(rec_atom.resnum)
                cur_sift = self.sifts[lig_name]
        else:
            new_sift = sift(lig_name)
            new_sift.add_sift_chunk(rec_atom.resnum)
            
            self.sifts[lig_name] = new_sift
            cur_sift = new_sift

        if structureutil.match_hbond(lig_atom,rec_atom,distance= dist):
            if rec_atom.atomic_number == 1:
                if 'H_DONOR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['H_DONOR'],  rec_atom.resnum, "exists")
            else:
                if 'H_ACCEPTOR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['H_ACCEPTOR'],  rec_atom.resnum, "exists")
    

    def add_sift_chunk_env(self,  rec_atom,lig_name,i_type, strength):
        if self.sifts.has_key(lig_name):
            if rec_atom.resnum in self.sifts[lig_name].sift.keys():
                cur_sift = self.sifts[lig_name]
            else:
                self.sifts[lig_name].add_sift_chunk(rec_atom.resnum)
                cur_sift = self.sifts[lig_name]
        else:
            new_sift = sift(lig_name)
            new_sift.add_sift_chunk(rec_atom.resnum)
            
            self.sifts[lig_name] = new_sift
            cur_sift = new_sift

        cur_sift.turn_sift_bit_on(self.bit_pos[i_type],  rec_atom.resnum, strength)

    def add_sift_chunk_ordinary(self,  rec_atom,lig_name, strength):
        if self.sifts.has_key(lig_name):
            if rec_atom.resnum in self.sifts[lig_name].sift.keys():
                cur_sift = self.sifts[lig_name]
            else:
                self.sifts[lig_name].add_sift_chunk(rec_atom.resnum)
                cur_sift = self.sifts[lig_name]
        else:
            new_sift = sift(lig_name)
            new_sift.add_sift_chunk(rec_atom.resnum)
            
            self.sifts[lig_name] = new_sift
            cur_sift = new_sift
        
        if 'CONTACT' in self.active_bits:
            cur_sift.turn_sift_bit_on(self.bit_pos['CONTACT'],  rec_atom.resnum, "exists")
        if int(rec_atom) in self.backbone_set:
            if 'BACKBONE' in self.active_bits:
                cur_sift.turn_sift_bit_on(self.bit_pos['BACKBONE'],  rec_atom.resnum, strength)
        else:
            if 'SIDECHAIN' in self.active_bits:
                cur_sift.turn_sift_bit_on(self.bit_pos['SIDECHAIN'],  rec_atom.resnum, strength)
            
            if int(rec_atom) in self.polar_set:
                if 'POLAR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['POLAR'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.hydrophobic_set:
                if 'HYDROPHOBIC' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['HYDROPHOBIC'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.aromatic_set:
                if 'AROMATIC' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['AROMATIC'],  rec_atom.resnum, strength)
            if int(rec_atom) in self.charged_set:
                if 'CHARGED' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['CHARGED'],  rec_atom.resnum, strength)


    def get_sift_string(self, lig_name):
        if self.used_bits_number != 13:
            out_sift = ''
            for chunk in self.sifts[lig_name].sift:
                for bit in self.active_bits:
                    out_sift.append(chunk.bit_set[bit])
            
            return ''.join(out_sift)
        else:
            """
            concatenate the bit string
            """
            out_sift = []
            for fp_chunk in self.sifts[lig_name].sift.keys():
                out_sift.append("".join(self.sifts[lig_name].sift[fp_chunk].bit_set))
                #for bit in self.sifts[lig_name].sift[fp_chunk].bit_set:
                    #out_sift += str(bit)
            return "".join(out_sift)

if __name__ == "__main__":
    rec_file = '/home/xiaohan/Downloads/protein/1CE1/1CE1_antibody.pdb'
    bind_file = '/home/xiaohan/Downloads/protein/1CE1/1CE1_antigen.pdb'

    from avg_sift import gen_avg_sift
    cutoff = 4.0
    fp_file = '/home/xiaohan/Desktop/1CE1_fp_%d.dat' %(cutoff)
    pat_file = '/home/xiaohan/Desktop/1CE1_pat_%d.dat' %(cutoff)
    gen_fp(rec_file,bind_file,fp_file,cutoff = cutoff)
    gen_avg_sift(fp_file,pat_file)


