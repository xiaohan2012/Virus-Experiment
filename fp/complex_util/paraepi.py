import os

from pickle import load, dump

from types import MethodType

import sys
import logging

from ve.config import data237_paraepi_cache as cache_dir

from common import make_logger
logger = make_logger("Paratope Epitope Detection")

#cache utility
from ve.fp.complex_util.cache import ParatopeCache as PCache, EpitopeCache as ECache

#solving trait dependency issue
from ve.fp.complex_util.propagate_distcache import DistanceCachePropagationTrait

class ParatopeNotFoundError(Exception):pass
class EpitopeNotFoundError(Exception):pass

class FindParaEpiTrait(DistanceCachePropagationTrait):
    """
    find epitope and paratope trait
    
    Dependency:
    1. dist cache propagation
    """

    def __init__(self, **kwargs):
        logger.debug("initializtin paratop epitope detection")
        super(FindParaEpiTrait,self).__init__(**kwargs)

        #paratope and epitope minimum distance
        self.paratope_threshold = 5
        self.epitope_threshold = 5

    def find_paratope(self, refresh = False):
        """
        (bool: False) => None

        refresh: mandate computation or not

        find paratope(residues in the antibody side that are close enough to the antigen side )
        """
        #computed already
        if hasattr(self, "paratope"):
            return self.paratope

        self.paratope = []

        logger.info("finding paratope for %s" %self.c_id)

        #propagating the distance cache to each reisdue
        self.prop_dist_cache()

        #if there is no need to refresh cache and the paratope result is in the cache
        #then load it from cache
        if not refresh and PCache.has_cache(self.c_id):
            self.paratope = PCache.load(self.c_id, self)
        else:
            #needs to be calculated
            #consider each residue  `b_r` in the antibody side
            for b_r in self.atb.residues:
                #if there exists one residue in the antigen side that is close enough to `b_r`, 
                #then `b_r` should belong to the paratope set
                for g_r in self.atg.residues:
                    dist =b_r.dist_to(g_r)
                    if dist <=self.paratope_threshold:
                        self.paratope.append(b_r)
                        break

            #if the paratope set is empty, there is something wrong
            if not self.paratope:
                raise ParatopeNotFoundError

            #save the result to cache
            PCache.dump(self.c_id, self.paratope)
            
        #be functional! return it
        return self.paratope

    def find_epitope(self, refresh=False):
        """
        find epitope(residues in the antigen side that are close enough to the antibody side )
        """
        if hasattr(self, "epitope"):
            return self.epitope
            
        self.epitope = []

        logger.info("finding epitope for %s" %self.c_id)

        #propagating the distance cache to each reisdue
        self.prop_dist_cache()

        #if there is not need to refresh and cache exists, then load it
        if not refresh and ECache.has_cache(self.c_id):
            self.epitope = ECache.load(self.c_id, self)
        else:
            #consider each residue  `g_r` in the antigen side
            for g_r in self.atg.residues:

                #if there exists one residue in the antibody side that is close enough to `g_r`, 
                #then `g_r` should belong to the epitope set
                for b_r in self.atb.residues:
                    if g_r.dist_to(b_r) <= self.epitope_threshold:
                        self.epitope.append(g_r)
                        break
            if not self.epitope:
                raise EpitopeNotFoundError
            #cache the result
            ECache.dump(self.c_id, self.epitope)

        #return it
        return self.epitope
        
#output the paratope and epitope into pdb format string
from ve.config import data237_complex_root as complex_path, data237_paraepi_root as paraepi_path

class ParaEpiIOTrait(FindParaEpiTrait):
        
    def gen_str(self, source_fp, res_ids):
        #load the complex file 
        #and filter out the unrelated lines
        lines = []
        for l in open(source_fp,"r").readlines():
            resnum = ''.join(l[22:26]).strip()
            if resnum in res_ids:
                lines.append(l)

        #return the concatenated string
        return "".join(lines)

    def epitope_str(self):
        #get the epitope residue ids
        res_ids = map(lambda a:str(a.resnum),self.find_epitope())

        return self.gen_str(self.atg.src, res_ids)

    def paratope_str(self):
        #get the paratope residue ids
        res_ids = map(lambda a:str(a.resnum),self.find_paratope())

        return self.gen_str(self.atb.src, res_ids)

    def write_epitope(self, paraepi_dir):
        #get the epitope string
        string = self.epitope_str()

        #get the epitope path
        output_fp = os.path.join(paraepi_dir, self.c_id, "epitope.pdb")

        #if path does not exist, create it
        if not os.path.exists(os.path.join(paraepi_dir, self.c_id)):
            os.makedirs(os.path.join(paraepi_dir, self.c_id))

        #write it
        f = open(output_fp, "w")
        f.write(string)
        f.close()

    def write_paratope(self, paraepi_dir):
        #get the paratope string
        string = self.paratope_str()

        #get the paratope path
        output_fp = os.path.join(paraepi_dir, self.c_id, "paratope.pdb")

        #if path does not exist, create it
        if not os.path.exists(os.path.join(paraepi_dir, self.c_id)):
            os.makedirs(os.path.join(paraepi_dir, self.c_id))

        #write it
        print output_fp
        f = open(output_fp, "w")
        f.write(string)
        f.close()