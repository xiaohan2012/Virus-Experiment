#!/usr/bin/env python

import os,  sys
from schrodinger import structure, structureutil
import schrodinger.utils.fileutils as fileutils
import numpy as np

from rb_tree import *



class distance_data(object):
    def __init__(self,  res_num = 0,  res_name = None,  coords = [0, 0, 0],  atom = None):
        
        self.res_num = res_num
        self.res_name = res_name
        self.coords = coords
        self.atom = atom

class distance_tree(object):
    
    def __init__(self,  pv_file = None):
        
        self.distances = red_black_tree()
        self.min_res = None
        self.max_res = None
        
        self.pv_file = pv_file
        
        if self.pv_file is not None:
            self.get_receptor_structure()
        else:
            self.receptor = structure.create_new_structure()
        
        self.fingerprints = sift_gen(self.receptor)


    def get_receptor_structure(self):
        
        if os.path.exists(self.pv_file):
            self.receptor = structure.StructureReader(self.pv_file).next()
        else:
            self.receptor = structure.create_new_structure()


    def set_receptor_structure(self,  struct):
        self.receptor = struct

    def set_pv_file(self,  filename):
        
        self.pv_file = filename


    def measure_dist(self,  end_point,  start_point = [0, 0, 0,]):
        """
        the distance between start point and end point
        """
        u = end_point
        v = start_point
        u = np.asarray(u, order='c')
        v = np.asarray(v, order='c')
        q=np.matrix(u-v)
        return np.sqrt((q*q.T).sum())


    def prepare_distance_data(self,  atom):
        
        res_name = atom.pdbres
        res_num = atom.resnum
        coords = [atom.x,  atom.y,  atom.z]
        
        data = distance_data(res_num,  res_name,  coords,  atom)
        
        return data


    def insert_atom(self,  distance,  data):
        
        self.distances.insert_node(distance,  data)


    def update_min_max(self,  res_num):
        
        #Update min
        if self.min_res is None:
            self.min_res = res_num
        elif self.min_res > res_num:
            self.min_res = res_num
        
        #Update max
        if self.max_res is None:
            self.max_res = res_num
        elif self.max_res < res_num:
            self.max_res = res_num


    def parse_receptor(self):
        
        for atom in self.receptor.atom:
            data = self.prepare_distance_data(atom)
            dist = self.measure_dist(data.coords)
            self.insert_atom(dist, data)
            

    def find_close_residues(self,  ligand,  cutoff = 4.0):
        
        for atom in ligand.atom:
            #the distance to original point?
            lig_dist = self.measure_dist([atom.xyz])

            #range should refer to the node level in the red black tree
            close_atoms = self.distances.find_nodes_in_range(lig_dist - cutoff, lig_dist + cutoff)

            
            for close in close_atoms:
                atom_atom_dist = self.measure_dist([atom.xyz],  close.data.coords)
                if atom_atom_dist <= cutoff:
                    self.update_min_max(close.data.res_num)
                    self.fingerprints.add_sift_chunk(close.data.atom,  atom,  ligand.title,  atom_atom_dist)

class sift_bitset(object):
    """ one single sift vector"""
    def __init__(self,  res_num = 0,  bit_set = None):
        
        self.res_num = res_num
        
        if bit_set == None:
            self.bit_set = [0 for x in range(9)]
        else:
            self.set_bitset(bit_set)


    def set_bitset(self,  sift_string):
        
        if len(sift_string) == 9:
            self.bit_set = sift_string
        else:
            print "set_bitset(): invalid sift length, skipping..."

    def turn_sift_bit_on(self,  sift_bit):
        
        try:
            self.bit_set[sift_bit] = 1
        except:
            print "Cannot turn bit " + str(sift_bit) + " on. Too bad..."


    def turn_sift_bit_off(self,  sift_bit):
        try:
            self.bit_set[sift_bit] = 0
        except:
            print "Cannot turn bit " + str(sift_bit) + " off. Too bad..."


class sift(object):
    """a set of sift_bitset vectors"""
    def __init__(self,  ligand_name = None , bitset_class = sift_bitset):
        self. sift = {}
        self.ligand_name = ligand_name
        self.bitset_class = bitset_class


    def turn_sift_bit_on(self,  sift_bit,  res_num):
        try:
            self.sift[res_num].turn_sift_bit_on(sift_bit)
        except:
            print "Cannot turn bit " + str(sift_bit) + " on. Too bad..."


    def turn_sift_bit_off(self,  sift_bit,  res_num):
        try:
            self.sift[res_num].turn_sift_bit_off(sift_bit)
        except:
            print "Cannot turn bit " + str(sift_bit) + " off. Too bad..."


    def add_sift_chunk(self,  res_num = 0,  bit_set = None):
        new_chunk = self.bitset_class(res_num,  bit_set)
        self.sift[res_num] = new_chunk
        


class sift_gen(object):
    
    def __init__(self,  receptor):
        
        self.POLAR_RESIDUES = ["ARG","ASP","GLU","HIS","ASN","GLN","LYS","SER","THR","ARN","ASH","GLH","HID","HIE","LYN"]
        self.HYDROPHOBIC_RESIDUES =["PHE","LEU","ILE","TYR","TRP","VAL","MET","PRO","CYS","ALA","CYX"]
        self.AROMATIC_RESIDUES = ["PHE","TYR","TRP","TYO"]
        self.CHARGED_RESIDUES = ["ARG","ASP","GLU","LYS","HIP","CYT","SRO","TYO","THO"]

        # Bit positions for each interaction:
        self.bit_pos = {'CONTACT' : 0,  'BACKBONE' : 1,  'SIDECHAIN' : 2,  'POLAR' : 3,  'HYDROPHOBIC' : 4,  'H_ACCEPTOR' : 5,  'H_DONOR' : 6,  'AROMATIC' : 7,  'CHARGED' : 8}
        self.active_bits = ['CONTACT',  'BACKBONE',  'SIDECHAIN',  'POLAR', 'HYDROPHOBIC',  'H_ACCEPTOR',  'H_DONOR',  'AROMATIC',  'CHARGED']
        self.max_bits_per_residue =  9 
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

    def add_sift_chunk(self,  rec_atom,  lig_atom,  lig_name,  dist):
        
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
            cur_sift.turn_sift_bit_on(self.bit_pos['CONTACT'],  rec_atom.resnum)
        if int(rec_atom) in self.backbone_set:
            if 'BACKBONE' in self.active_bits:
                cur_sift.turn_sift_bit_on(self.bit_pos['BACKBONE'],  rec_atom.resnum)
        else:
            if 'SIDECHAIN' in self.active_bits:
                cur_sift.turn_sift_bit_on(self.bit_pos['SIDECHAIN'],  rec_atom.resnum)
            
            if int(rec_atom) in self.polar_set:
                if 'POLAR' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['POLAR'],  rec_atom.resnum)
            if int(rec_atom) in self.hydrophobic_set:
                if 'HYDROPHOBIC' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['HYDROPHOBIC'],  rec_atom.resnum)
            if int(rec_atom) in self.aromatic_set:
                if 'AROMATIC' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['AROMATIC'],  rec_atom.resnum)
            if int(rec_atom) in self.charged_set:
                if 'CHARGED' in self.active_bits:
                    cur_sift.turn_sift_bit_on(self.bit_pos['CHARGED'],  rec_atom.resnum)
            if structureutil.match_hbond(lig_atom,rec_atom,distance= dist):
                if rec_atom.atomic_number == 1:
                    if 'H_DONOR' in self.active_bits:
                        cur_sift.turn_sift_bit_on(self.bit_pos['H_DONOR'],  rec_atom.resnum)
                else:
                    if 'H_ACCEPTOR' in self.active_bits:
                        cur_sift.turn_sift_bit_on(self.bit_pos['H_ACCEPTOR'],  rec_atom.resnum)


    def get_sift_string(self, lig_name):
        out_sift = ''
        if self.used_bits_number != 9:
            for chunk in self.sifts[lig_name].sift:
                for bit in self.active_bits:
                    out_sift.append(chunk.bit_set[bit])
            
            return ''.join(out_sift)
        else:
            """
            concatenate the bit string
            """
            for fp_chunk in self.sifts[lig_name].sift.keys():
                for bit in self.sifts[lig_name].sift[fp_chunk].bit_set:
                    out_sift += str(bit)
            return out_sift



def gen_fp_to_file(receptor=None,binder=None,fp_path=''):
    rec_tree = distance_tree()

    rec_tree.set_receptor_structure(receptor)
    rec_tree.parse_receptor()

    rec_tree.find_close_residues(binder)

    with open(fp_path,  'w') as out_fp:
        for key in rec_tree.fingerprints.sifts.keys():
            rec_tree.fingerprints.fill_missing_zeros(rec_tree.min_res,  rec_tree.max_res,  key)
            fp_string = rec_tree.fingerprints.get_sift_string(key)
            out_fp.write(rec_tree.receptor.title + ':' + key + ':' + str(rec_tree.min_res) + ':' + fp_string + '\n')

    print 'saved to',fp_path

    return fp_path

if __name__ == "__main__":
    rec_file = 'protein2protein/1A2Y_antigen.pdb'
    bind_file = 'protein2protein/1A2Y_antibody.pdb'

    basename,  ext = fileutils.splitext(os.path.basename(rec_file))

    rec_tree = distance_tree()

    #receptor data
    rec = structure.StructureReader(rec_file).next()

    if rec.title == "":
        rec.title = basename

    rec_tree.set_receptor_structure(rec)
    rec_tree.parse_receptor()

    #bind data
    basename,  ext = fileutils.splitext(os.path.basename(rec_file))
    bind = structure.StructureReader(bind_file).next()

    if bind .title == "":
        bind .title = basename

    
    rec_tree.find_close_residues(bind)


    out_fp = open(rec_tree.receptor.title + '_fp.dat',  'w')


	  
    for key in rec_tree.fingerprints.sifts.keys():
        rec_tree.fingerprints.fill_missing_zeros(rec_tree.min_res,  rec_tree.max_res,  key)
        fp_string = rec_tree.fingerprints.get_sift_string(key)
        out_fp.write(rec_tree.receptor.title + ':' + key + ':' + str(rec_tree.min_res) + ':' + fp_string + '\n')

    out_fp.close()
