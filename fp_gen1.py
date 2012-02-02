#!/usr/bin/env python

import os,  sys
from schrodinger import structure, structureutil
import schrodinger.utils.fileutils as fileutils
import numpy as np

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
        while current != self.sentinel and self.compare_keys(up, current.key) >= 0:
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
    def __init__(self,  ligand_name = None):
        
        self. sift = {}
        self.ligand_name = ligand_name


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
        
        new_chunk = sift_bitset(res_num,  bit_set)
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



def gen_fp(receptor=None,binder=None,fp_path=''):
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
