#!/usr/bin/env python
import os,  sys
from schrodinger import structure, structureutil
from collections import defaultdict
from cgkit.cgtypes import *
from numpy import *

class Residue(object):
    def __init__(self,res,comp):
        self.c = comp
        self.fp = [0] * 110
        self.ca = None
        self.body = res

        self.y_cords = linspace(-30,10,9)
        self.x_cords = linspace(0,20,11)

        self.x_min , self.x_max = 0.,20.
        self.x_step = 2.
        self.x_seg_cnt = ( self.x_max - self.x_min ) / self.x_step
        self.x_ind_min , self.x_ind_max = self.x_min / self.x_step, self.x_max / self.x_step - 1

        self.y_min , self.y_max = -30.,10.
        self.y_step = 5.
        self.y_seg_cnt = ( self.y_max - self.y_min ) / self.y_step
        self.y_ind_min , self.y_ind_max = self.y_min / self.y_step, self.y_max / self.y_step - 1

        self.orig_point = vec3(0,0,0)

        self.dist_min = 0
        self.dist_max = 20
        self.dist_step = 2
        self.dist_ind_min = self.dist_min / self.dist_step
        self.dist_ind_max = self.dist_max / self.dist_step - 1
        for a in res.atom:
            if a.pdbname.strip() == "CA":
                self.ca = a
                break
        if self.ca is None:
            raise ValueError("self.ca should not be none")
    
    def turn_on_bit(self,bit_num,count):
        self.fp[bit_num] = count

    def get_surrounding_res(self):
        for res in self.c.residues:
            if res.ca == self.ca:continue
            yield res

    def get_region_number(self,res_ca):
        my_point = vec3(self.ca.xyz)
        v1 = my_point - self.orig_point 
        v2 = vec3(res_ca.xyz) - my_point
        ang = v1.angle(v2)

        v2_len = v2.length()
        x = v2_len * sin(ang)
        y = v2_len * cos(ang)
        
        if x < 0:
            raise ValueError("x cannot be negative")

        x_ind  = floor( x / self.x_step)
        y_ind  = floor( y / self.y_step)
        
        if x_ind >= self.x_ind_min and x_ind <= self.x_ind_max and\
            y_ind >= self.y_ind_min and y_ind <= self.y_ind_max:
            #within the range
            return int((y_ind - self.y_ind_min) * self.x_seg_cnt + x_ind)
        else:
            return None


    def get_struct_fp(self):
        #type one fp
        d_ = defaultdict(int)
        for res in self.get_surrounding_res():
            num = self.get_region_number(res.ca)
            if num is not None:
                d_[num] += 1
        
        for bit_num,count in d_.iteritems():
            self.turn_on_bit(bit_num,count)

        return self.fp           

    def get_surrounding_fp(self):
        #type two fp
        d_ = defaultdict(list)
        for res in self.get_surrounding_res():
            ca = res.ca
            dist = sqrt(sum(power(array(ca.xyz) - array(self.ca.xyz),2)))#distance
            dist_ind = floor( dist / self.dist_step )

            if dist_ind > self.dist_ind_min and dist_ind < self.dist_ind_max:
                #within range
                d_[dist_ind].append(res)

        hydro_dict={'A':0.61,'C':1.07,'D':0.46,'E':0.47,'F':2.02,'G':0.07,'H':0.61,'I':2.22,'K':1.15,'L':1.53,'M':1.18,'N':0.06,'P':1.95,'Q':0.0,'R':0.6,'S':0.05,'T':0.05,'V':1.32,'W':2.65,'Y':1.88}
        charged_dict={'A':-0.01,'C':0.12,'D':0.15,'E':0.07,'F':0.03,'G':0.0,'H':0.08,'I':-0.01,'K':0.0,'L':-0.01,'M':0.04,'N':0.06,'P':0.0,'Q':0.05,'R':0.04,'S':0.11,'T':0.04,'V':0.01,'W':0.0,'Y':0.03}
        h_bond_dict={'A':0,'C':0,'D':1,'E':1,'F':0,'G':0,'H':1,'I':0,'K':2,'L':0,'M':0,'N':2,'P':0,'Q':2,'R':4,'S':1,'T':1,'V':0,'W':1,'Y':1}
        for i in xrange(self.dist_ind_max - self.dist_ind_min + 1):
            # for layer i
            h_bond, charged , hydro = 0 , 0 , 0
            for res in d_[i]:
                # for res in layer i
                code = res.body.getCode()
                hydro += hydro_dict[code]
                charged += charged_dict[code]
                h_bond += h_bond_dict[code]
            
            #fp for layer i,in the 3 aspects
            self.turn_on_bit(80 + i , hydro)
            self.turn_on_bit(90 + i , charged)
            self.turn_on_bit(100 + i , h_bond)

    def __repr__(self):            
        return "ca atom index:%d" %(self.ca.index)

class Complex(object):
    def __init__(self,pdb_fp):
        self.st = structure.StructureReader(pdb_fp).next()
        self.residues = []
        for res in self.st.residue:
            self.residues.append(Residue(res,self))

    def get_fp(self):
        self.fp_list = []
        for i , res in enumerate(self.residues):
            print "residue %d" %i
            res.get_surrounding_fp()
            res.get_struct_fp()
            print res.fp,len(res.fp)
            self.fp_list.append(res.fp)
            
    def write_fp_to_file(self,path):
        with open(path,'w') as f:
            temp_list = []
            for fp in self.fp_list:
                temp_list.append(",".join(map(str,fp)))
            f.write("\n".join(temp_list))
            
if __name__ == "__main__":
    c = Complex("new_fp_data/data/data HIV/1E6J/1E6J.pdb")
    c.get_fp()
    c.write_fp_to_file("/home/xiaohan/Desktop/temp.dat")
