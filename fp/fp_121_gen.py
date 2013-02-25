# -*- coding: utf-8 -*- 
"""
120 bits finger print

"""
from __future__ import division

import os
import numpy as np

from glob import glob
from collections import defaultdict,OrderedDict
from itertools import izip
from cgkit.cgtypes import *

from ve.util.load_pdb import load_pdb_struct
from ve.util.file import fpstr2file
from ve.config import *

from geom import *

logging.basicConfig(level = logging.DEBUG)

class FP121(object):
    def __init__(self,antigen,antibody):
        self.antigen = antigen
        self.antibody = antibody
        
        self.prange = 6#range on the perpendicular direction
        self.pstep = 1

        self.crange = 20#circular range
        self.cstep = 2

        print("calculating antigen and antibody center")
        self.atg_center, self.atb_center = None, None
        self.atg_center,self.atb_center = self.get_both_centers()
       
        self.center_vec = Vector(self.atg_center - self.atb_center)

        print("calculating closest pair and set their center")
        self.set_closest_residue_center()
       
        print("calculating axial plane")
        self.set_axial_plane()

    def set_closest_residue_center(self):
        self.rc_cache = {}
        self.closest_pair = None;
        self.closest_pair = self.get_closest_residue_pair()

        closest_antibody_coords = self.get_r_center(self.closest_pair[0],"b")
        closest_antigen_coords = self.get_r_center(self.closest_pair[1],"g")

        pair_coords = np.array([closest_antibody_coords, closest_antigen_coords])
        self.cp_center = Point(np.average(np.array(pair_coords),0))
        print("closest pair center %s" %repr(self.cp_center))
 
    def set_axial_plane(self):
        self.complex_center = get_perp_point(self.atb_center, self.atg_center, self.cp_center)
        self.axial_plane = get_perp_plane(self.center_vec, self.cp_center)

    def get_group_from_dist(self,pdist,cdist):
        """parameter,dist to plane and dist to center"""
        if pdist <= self.prange and cdist <= self.crange:
            pidx = int(pdist / self.pstep)
            cidx = int(cdist / self.cstep)
            return pidx,cidx

    def gen_all_fp(self):
        self.antigen_fps = OrderedDict((r,self.gen_r_fp(r,"b")) 
                                        for r in self.antigen.residue)

    def get_r_fp_str(self,r_fp, delimiter = ","):
        return delimiter.join(map(lambda i: "%d" %i , r_fp))

    def get_fp_str(self,delimiter = "\n"):
        return delimiter.join(map(lambda tp: "%d: %s" %(tp[0].resnum,self.get_r_fp_str(tp[1])),
                                    self.antigen_fps.items()))

    def gen_r_fp(self,residue,g_or_b):
        """g_or_b refers to antigen or antibody"""
        target_residue_center = self.get_r_center(residue, g_or_b)
        perp_center = self.axial_plane.get_perp_point(target_residue_center)
        center_line = Line(perp_center, target_residue_center)

        #first 60 for antigen, next 60 for antibody
        antigen_side_fp = defaultdict(int)
        antibody_side_fp = defaultdict(int)

        for r in self.antigen.residue:
            residue_center = self.get_r_center(r,"g")
            dist2plane = self.axial_plane.dist2point(residue_center)
            dist2perp_line = center_line.dist2point(residue_center)
            if dist2plane <= self.prange and dist2perp_line <= self.crange:
                yidx = dist2plane // self.pstep
                xidx = dist2perp_line // self.cstep
                idx = (xidx,yidx)
                #antigen_side_fp[idx].append(r)
                antigen_side_fp[idx] += 1

        for r in self.antibody.residue:
            residue_center = self.get_r_center(r,"b")
            dist2plane = self.axial_plane.dist2point(residue_center)
            dist2perp_line = center_line.dist2point(residue_center)
            if dist2plane <= self.prange and dist2perp_line <= self.crange:
                yidx = dist2plane // self.pstep
                xidx = dist2perp_line // self.cstep
                idx = (xidx,yidx)
                #antigen_side_fp[idx].append(r)
                antibody_side_fp[idx] += 1
        fp = [0] * 120
        for yidx in xrange(self.prange // self.pstep):
            for xidx in xrange(self.crange // self.cstep):
                tile_size = self.crange // self.cstep
                idx = yidx * tile_size + xidx
                fp[idx] = antigen_side_fp[(xidx,yidx)]
                fp[60 + idx] = antibody_side_fp[(xidx,yidx)]
        return fp

    def get_both_centers(self):
        if not self.atg_center and not self.atb_center:
            atg_xyz = np.array([np.array(a.xyz) for res in self.antigen.residue for a in res.atom ])
            atb_xyz = np.array([np.array(a.xyz) for res in self.antibody.residue for a in res.atom])
            return Point(np.average(atg_xyz, axis = 0)), Point(np.average(atb_xyz, axis = 0))
        else:
            return Point(self.atg_center), Point(self.atb_center)

    def get_r_center(self,residue,g_or_b):
        """get residue center"""
        key = "%s%d" %(g_or_b, residue.resnum)
        if self.rc_cache.has_key(key):
            return self.rc_cache[key]
        else:
            cent = np.average([a.xyz for a in residue.atom],axis = 0)
            self.rc_cache[key] = cent
            return cent


    def get_rr_dist(self,r1,r2):
        rc1,rc2 = self.get_r_center(r1,"b"), self.get_r_center(r2,"g")
        return np.sqrt(np.sum(np.power(np.array(rc1) - np.array(rc2),2)))

    def get_closest_residue_pair(self):
        if not self.closest_pair:
            min_dist = float("inf")
            closest_pair = None
            for ab in self.antibody.residue:
                for ag in self.antigen.residue:
                    dist = self.get_rr_dist(ab, ag)
                    if dist < min_dist:
                        closest_pair = ab,ag
                        min_dist = dist
        return closest_pair


####Real life date processing part#####
def data237_fp_gen(refresh = True):
    path = os.path.join(data237_root ,"splitted_complex","*")
    print path
    for fp in glob(path):
        complex_id = os.path.basename(fp)
        
        #prepare neccessary directory
        if not os.path.exists(data237_fp_root): os.makedirs(data237_fp_root)

        fp_path = os.path.join(data237_fp_root,"%s.fp" %complex_id)

        if not refresh and os.path.exists(fp_path):
            print "%s processed" %complex_id
        else:
            print "processing %s" %complex_id
            atb_path = os.path.join(fp,"antibody.pdb")
            atg_path = os.path.join(fp,"antigen.pdb")
            antibody = load_pdb_struct(atb_path)
            antigen = load_pdb_struct(atg_path)
        
            fingerprint = FP121(antigen,antibody)
            fingerprint.gen_all_fp()
            
            fpstr2file(fingerprint.get_fp_str(),fp_path)
            print "fingerprint saved"

if __name__ == "__main__":
    data237_fp_gen(False)
