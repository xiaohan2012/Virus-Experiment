"""
finger print for interative force
"""
from __future__ import division

from common import *
from ve.fp.fp_75_gen import get_15bits 

import sys
import logging

logging.basicConfig( stream=sys.stderr )
logging.getLogger("Interactive Force").setLevel( logging.DEBUG )
logger = logging.getLogger("Interactive Force")

def gen_if_fp_atg(self):
    """
    interactive force finger print, antigen as the receptor
    for each residue
    """
    logger.info("interactive force finger print, antigen as the receptor")
    return get_15bits(receptor= self.atg, binder= self.atb)

def gen_if_fp_atb(self):
    """
    interactive force finger print, antibody as the receptor
    for each residue
    """
    logger.info("interactive force finger print, antibody as the receptor")
    return get_15bits(receptor= self.atb, binder= self.atg)

def gen_if_complex_fp(self, res_fps, residues):
    """the interactive finger print for the *compelx*
    (dict of (Residue -> ResidueFingerprint), list of Residue) -> CompelxFingerprint
    given the fps for residues and residue list, output the complex fp
    """

    from collections import defaultdict

    #list of finger print of length 15
    from ve.fp.fp import HeadlessFingerprint
    fp_lists = map(lambda i: HeadlessFingerprint(res_fps.get_bitlength()), xrange(self.group_count))
    
    #for each residue in this side
    for residue in residues:
        #get its distance to the compelx center
        dist = self.get_geom_center().dist2point(residue.center)

        #within range
        if dist < self.boundary:
            #group it based on the center
            group_id = int(dist // self.step)
            
            #vector addition of the fingerprints
            fp_lists[group_id] += res_fps[residue]
            
    #init the complex fp with length 0
    complex_fp  = HeadlessFingerprint(0)
    
    #concatenate the elements in the fp list
    for fp in fp_lists:
        complex_fp = complex_fp.append(fp)

    #return it
    return complex_fp

def gen_if_complex_fp_atg(self):
    """ complex interactive force fp for antigen side"""
    #get the finger prints for residues 
    res_fps = self.gen_if_residue_fp_atg()
    
    #get the residues
    residues = self.atg.residues
    
    #return it
    return self.gen_if_complex_fp(res_fps, residues)
    
def gen_if_complex_fp_atb(self):
    """ complex interactive force fp for antibody side"""
    #get the finger prints for residues 
    res_fps = self.gen_if_residue_fp_atb()
    
    #get the residues
    residues = self.atb.residues
    
    #return it
    return self.gen_if_complex_fp(res_fps, residues)
    


def init_interative_force_fp_trait(self, dist=20, step=2):
    """
    Dependecy:
    1, atg and atb
    """
    #check atg and atb availablity dependency
    if not is_atg_atb_loaded(self):
        raise NoAtgNorAtbLoadedError
    
    from types import MethodType
    
    #residue case
    self.gen_if_residue_fp_atg = MethodType(gen_if_fp_atg, self)
    self.gen_if_residue_fp_atb = MethodType(gen_if_fp_atb, self)

    #compelx case
    self.gen_if_complex_fp = MethodType(gen_if_complex_fp, self)
    self.gen_if_complex_fp_atg = MethodType(gen_if_complex_fp_atg, self)
    self.gen_if_complex_fp_atb = MethodType(gen_if_complex_fp_atb, self)

    #we need to complex center
    from ve.fp.complex_util.geom import init_geom_trait
    init_geom_trait(self)
    
    #geometry related calculation
    self.boundary = dist
    self.step = step
    self.group_count = self.boundary // self.step
    
