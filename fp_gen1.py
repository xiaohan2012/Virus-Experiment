#!/usr/bin/env python

import os,  sys
from schrodinger import structure, structureutil
import schrodinger.utils.fileutils as fileutils
import numpy as np
from collections import defaultdict

BLACK = 0
RED = 1


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



                    


class tree_node(object):
    def __init__(self,  key = None,  data = None,  color = RED):
        self.key = key
        self.data = data
        self.parent = self.left_child = self.right_child = None
        self.color = color
        self.nonzero = 1
        self.count = 1

    def __nonzero__(self):
        return self.nonzero

    def is_red(self):
        if self.color == RED:
            return 1
        else:
            return 0

    def is_black(self):
        if self.color == BLACK:
            return 1
        else:
            return 0

    def append_data(self,  additional_data):
        for elem in self.data:
            if additional_data == elem:
                print "Item already in data list... Ignoring"
            else:
                data.append(additional_data)
                self.count += 1

    def delete_data(self,  del_data):
        if del_data in self.data:
            self.data.remove(del_data)
            self.count -= 1



class red_black_tree(object):
    """
    Class for building and modifying binary red black tree
    """
    def __init__(self):
        self.sentinel = tree_node()
        self.sentinel.left = self.sentinel.right = self.sentinel
        self.sentinel.color = BLACK
        self.sentinel.nonzero = 0
        self.root = self.sentinel
        self.elements = 0

    def compare_keys(self,  key_x,  key_y):
        return cmp(key_x,  key_y)

    def rotate_left(self, node_x):

        node_y = node_x.right

        node_x.right = node_y.left
        if node_y.left != self.sentinel:
            node_y.left.parent =node_x

        if node_y != self.sentinel:
            node_y.parent = node_x.parent
        if node_x.parent:
            if node_x == node_x.parent.left:
                node_x.parent.left = node_y
            else:
                node_x.parent.right = node_y
        else:
            self.root = node_y

        node_y.left = node_x
        if node_x != self.sentinel:
            node_x.parent = node_y


    def rotate_right(self, node_x):
        
        node_y = node_x.left

        node_x.left = node_y.right
        if node_y.right != self.sentinel:
            node_y.right.parent = node_x

        if node_y != self.sentinel:
            node_y.parent = node_x.parent
        if node_x.parent:
            if node_x == node_x.parent.right:
                node_x.parent.right = node_y
            else:
                node_x.parent.left =node_y
        else:
            self.root = node_y

        node_y.right = node_x
        if node_x != self.sentinel:
            node_x.parent = node_y


    def insert_rebalance(self, node_x):

        while node_x != self.root and node_x.parent.is_red():
            if node_x.parent == node_x.parent.parent.left:
                node_y = node_x.parent.parent.right

                if node_y.is_red():
                    node_x.parent.color = BLACK
                    node_y.color = BLACK
                    node_x.parent.parent.color = RED
                    node_x = node_x.parent.parent

                else:
                    if node_x == node_x.parent.right:
                        node_x =node_x.parent
                        self.rotate_left(node_x)

                    node_x.parent.color = BLACK
                    node_x.parent.parent.color = RED
                    self.rotate_right(node_x.parent.parent)
            else:
                node_y = node_x.parent.parent.left

                if node_y.is_red():
                    node_x.parent.color = BLACK
                    node_y.color = BLACK
                    node_x.parent.parent.color = RED
                    node_x = node_x.parent.parent

                else:
                    if node_x == node_x.parent.left:
                        node_x = node_x.parent
                        self.rotate_right(node_x)

                    node_x.parent.color = BLACK
                    node_x.parent.parent.color = RED
                    self.rotate_left(node_x.parent.parent)

        self.root.color = BLACK


    def insert_node(self, key, data):
        
        current = self.root
        parent = None
        
        while current != self.sentinel:
            rc = self.compare_keys(key,  current.key)
            if rc == 0:
                current.append_data(data)
                return current
            parent = current
            if rc < 0:
                current = current.left
            else:
                current = current.right

        new_node = tree_node(key, data)
        new_node.left = new_node.right = self.sentinel
        new_node.parent = parent

        self.elements = self.elements + 1

        if parent:
            if self.compare_keys(key, parent.key) < 0:
                parent.left = new_node
            else:
                parent.right = new_node
        else:
            self.root = new_node

        self.insert_rebalance(new_node)
        return new_node


    def delete_rebalance(self, node_x):
        
        while node_x != self.root and node_x.is_black():
            if node_x == node_x.parent.left:
                node_w = node_x.parent.right
                if node_w.is_red():
                    node_w.color = BLACK
                    node_x.parent.color = RED
                    self.rotate_left(node_x.parent)
                    node_w = node_x.parent.right

                if node_w.left.is_black() and node_w.right.is_black():
                    node_w.color = RED
                    node_x = node_x.parent
                else:
                    if node_w.right.is_black():
                        node_w.left.color = BLACK
                        node_w.color = RED
                        self.rotate_right(node_w)
                        node_w = node_x.parent.right

                    node_w.color = node_x.parent.color
                    node_x.parent.color = BLACK
                    node_w.right.color = BLACK
                    self.rotate_left(node_x.parent)
                    node_x = self.root

            else:
                node_w = node_x.parent.left
                if node_w.is_red():
                    node_w.color = BLACK
                    node_x.parent.color = RED
                    self.rotate_right(node_x.parent)
                    node_w = node_x.parent.left

                if node_w.right.is_black() and node_w.left.is_black():
                    node_w.color = RED
                    node_x = node_x.parent
                else:
                    if node_w.left.is_black():
                        node_w.right.color = BLACK
                        node_w.color = RED
                        self.rotate_left(node_w)
                        node_w = node_x.parent.left

                    node_w.color = node_x.parent.color
                    node_x.parent.color = BLACK
                    node_w.left.color = BLACK
                    self.rotate_right(node_x.parent)
                    node_x = self.root

        node_x.color = BLACK

    def delete_node(self,  convict,  data = None):
        
        if not convict or convict == self.sentinel:
            return
            
        if convict.count > 1 and not data: 
            convict.delete_data(data)
            return

        if convict.left == self.sentinel or convict.right == self.sentinel:
            y = convict
        else:
            y = convict.right
            while y.left != self.sentinel:
                y = y.left

        if y.left != self.sentinel:
            x = y.left
        else:
            x = y.right

        x.parent = y.parent
        if y.parent:
            if y == y.parent.left:
                y.parent.left = x
            else:
                y.parent.right = x
        else:
            self.root = x

        if y != convict:
            convict.key = y.key
            convict.data = y.data

        if y.color == BLACK:
            self.delete_rebalance(x)

        del y
        self.elements = self.elements - 1


    def delete_node_by_key(self,  key,  data = None):
        
        del_node = self.find_node(key)
        
        self.delete_node(del_node,  data)


    def find_node(self, key):

        current = self.root

        while current != self.sentinel:
            rc = self.compare_keys(key, current.key)
            if rc == 0:
                return current
            else:
                if rc < 0:
                    current = current.left
                else:
                    current = current.right

        return None


    def find_nodes_in_range(self,  down,  up):
        
        out = []
        current = self.root
        
        cmp_down = self.compare_keys(down, current.key)
        
        #Searching for lowermost item
        while current.left != self.sentinel and self.compare_keys(down, current.left.key) <= 0:
            if self.compare_keys(down, current.key) <= 0:
                current = current.left
            else:
                current = current.right
        #Trawersing right until reaching rightmost item within search criteria
        while current and current != self.sentinel and self.compare_keys(up, current.key) >= 0:
            out.append(current)
            current = self.next_node(current)
        return out



    def list_nodes(self):
        cur = self.first_node()
        result = []
        while cur:
            result.append(cur)
            cur = self.next_node(cur)
        return result


    def first_node(self):
        cur = self.root
        while cur.left:
            cur = cur.left
        return cur


    def last_node(self):
        cur = self.root
        while cur.right:
            cur = cur.right
        return cur


    def next_node(self, prev):
        cur = prev
        if cur.right:
            cur = prev.right
            while cur.left:
                cur = cur.left
            return cur
        while 1:
            cur = cur.parent
            if not cur:
                return None
            if self.compare_keys(cur.key, prev.key)>=0:
                return cur


    def prev_node(self, next):
        cur = next
        if cur.left:
            cur = next.left
            while cur.right:
                cur = cur.right
            return cur
        while 1:
            cur = cur.parent
            if cur is None:
                return None
            if self.__cmp(cur.key, next.key)<0:
                return cur


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


class sift(object):
    """a set of sift_bitset vectors"""
    def __init__(self,  ligand_name = None):
        
        self. sift = {}
        self.ligand_name = ligand_name


    def turn_sift_bit_on(self,  sift_bit,  res_num, strength):
        
        try:
            self.sift[res_num].turn_sift_bit_on(sift_bit, strength)
        except:
            print "Cannot turn bit " + str(sift_bit) + " on. Too bad..."


    def turn_sift_bit_off(self,  sift_bit,  res_num):
        try:
            self.sift[res_num].turn_sift_bit_off(sift_bit)
        except:
            print "Cannot turn bit " + str(sift_bit) + " off. Too bad..."


    def add_sift_chunk(self,  res_num = 0,  bit_set = None):
        
        new_chunk = sift_bitset(res_num,  bit_set)
        self.sift[res_num] = new_chunk
        


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




def gen_fp(receptor_file="",binder_file="",fp_path='',cutoff = 4.0):
    rec_tree = distance_tree(receptor_file)
    rec_tree.parse_receptor()

    #print 'cutoff %f' %cutoff
    binder= structure.StructureReader(binder_file).next() #read the binder
    rec_tree.find_close_residues(binder,cutoff = cutoff)

    with open(fp_path,  'w') as out_fp:
        for key in rec_tree.fingerprints.sifts.keys():
            rec_tree.fingerprints.fill_missing_zeros(rec_tree.min_res,  rec_tree.max_res,  key)
            fp_string = rec_tree.fingerprints.get_sift_string(key)
            out_fp.write(rec_tree.receptor.title + ':' + key + ':' + str(rec_tree.min_res) + ':' + fp_string + '\n')

    #print 'saved to',fp_path
    
    return fp_path
if __name__ == "__main__":
    rec_file = '/home/xiaohan/Downloads/protein/1CE1/1CE1_antibody.pdb'
    bind_file = '/home/xiaohan/Downloads/protein/1CE1/1CE1_antigen.pdb'

    from avg_sift import gen_avg_sift
    cutoff = 4.0
    fp_file = '/home/xiaohan/Desktop/1CE1_fp_%d.dat' %(cutoff)
    pat_file = '/home/xiaohan/Desktop/1CE1_pat_%d.dat' %(cutoff)
    gen_fp(rec_file,bind_file,fp_file,cutoff = cutoff)
    gen_avg_sift(fp_file,pat_file)


