import os

from pickle import load, dump

from types import MethodType

import sys
import logging

from ve.config import data237_paraepi_cache as cache_dir

logging.basicConfig( stream=sys.stderr )
logging.getLogger("Find ParaEpi").setLevel( logging.DEBUG )
logger = logging.getLogger("Find ParaEpi")

class ParatopeNotFoundError(Exception):pass
class EpitopeNotFoundError(Exception):pass


def find_paratope(self, refresh = False):
    """
    (bool: False) => None
    
    refresh: mandate computation or not

    find paratope(residues in the antibody side that are close enough to the antigen side )
    """
    #computed already
    if self.paratope is not None:
        return self.paratope
        
    self.paratope = []

    logger.info("finding paratope for %s" %self.c_id)

    #propagating the distance cache to each reisdue
    self.prop_dist_cache()

    attr_name = "paratope"
    
    #if there is no need to refresh cache and the paratope result is in the cache
    #then load it from cache
    if not refresh and self.cache_exists(attr_name):
        self.paratope = self.load_cache(attr_name)
    #needs to be calculated
    else:
        #consider each residue  `b_r` in the antibody side
        for b_r in self.atb.residues:
            #if there exists one residue in the antigen side that is close enough to `b_r`, 
            #then `b_r` should belong to the paratope set
            for g_r in self.atg.residues:
                dist =b_r.dist2residue(g_r)
                if dist <=self.paratope_threshold:
                    self.paratope.append(b_r)
                    break

        #if the paratope set is empty, there is something wrong
        if not self.paratope:
            raise ParatopeNotFoundError

        #save the result to cache
        self.save_cache(attr_name)
        
    #be functional! reture it
    return self.paratope

def find_epitope(self, refresh=False):
    """
    find epitope(residues in the antigen side that are close enough to the antibody side )
    """
    if self.epitope is not None:
        return self.epitope
    
    self.epitope = []
    
    logger.info("finding epitope for %s" %self.c_id)

    #propagating the distance cache to each reisdue
    self.prop_dist_cache()

    attr_name = "epitope"
    
    #if there is not need to refresh and cache exists, then load it
    if not refresh and self.cache_exists(attr_name):
        self.epitope = self.load_cache(attr_name)
    else:
        #consider each residue  `g_r` in the antigen side
        for g_r in self.atg.residues:

            #if there exists one residue in the antibody side that is close enough to `g_r`, 
            #then `g_r` should belong to the epitope set
            for b_r in self.atb.residues:
                if g_r.dist2residue(b_r) <= self.paratope_threshold:
                    self.epitope.append(g_r)
                    break
        if not self.epitope:
            raise EpitopeNotFoundError
        #cache the result
        self.save_cache(attr_name)

    #return it
    return self.epitope



def cache_exists(self,cache_type):
    """
    (string) => bool
    given the `cache_type`(paratope or epitope, currently),
    checks if the corresponding cache exists for the given complex
    """
    
    #the cache holding directory
    par_dir = os.path.join(cache_dir, cache_type)
    if not os.path.exists(par_dir):
        os.makedirs(par_dir)
    
    #the cache object path
    cache_path = os.path.join(cache_dir, cache_type, self.c_id)

    #check if it exists
    mark = os.path.exists(cache_path)
    
    if mark:
        logger.info("%s cache for %s exists" %(cache_type, self.c_id))
    
    #return it
    return mark

def paratope_cl_process(self,content):
    """
    (list of int) => (list of Residue)

    paratope cache load helper function    
    given a list of residue numbers, return the corresponding list of residue objects of the paratope side

    """
    return filter(lambda r: r.resnum in content, 
                  self.atb.residues)

def epitope_cl_process(self,content):
    """
    (list of int) => (list of Residue)

    epitope cache load helper function    
    given a list of residue numbers, return the corresponding list of residue objects of the epitope side
    """
    return filter(lambda r: r.resnum in content, 
                  self.atg.residues)

def epi_para_cd_process(self,obj):
    """
    (list of Residue) => (list of int)

    epitope and paratope cache dump helper function
    Given a list of Residue, return the corresponding list of residue id
    """
    #extract the residue ids in the list of residues
    return map(lambda r: r.resnum, obj)


def load_cache(self,cache_type):
    """
    (string) => list of Residue
    Given the cache type(paratope or epitope), return the cacheed residue list specified by the complex_id of `self`
    """
    #get the cache path
    cache_path = os.path.join(cache_dir, cache_type, self.c_id)
    
    #get the cached residue ids
    cache_content = load(open(cache_path,"r"))
    
    #return the residue list
    return self.cache_load_methods[cache_type](cache_content)

def save_cache(self,cache_type):
    """
    (string) => None
    Given the cache type(paratope or epitope), save the cached residue id list specified by the complex_id of `self`
    """
    #get the cache path
    cache_path = os.path.join(cache_dir, cache_type, self.c_id)
    
    #get the residue ids
    val = self.cache_dump_methods[cache_type](getattr(self,cache_type))
    logger.debug("Cache residue ids: %s" %" ".join(map(lambda i: ":%d" %i, val)))
    #dump it
    dump(val, open(cache_path,"w"))

    
def init_find_epiparatope_trait(self):
    """
    find epitope and paratope util
    
    Dependency:
    1. atg and atb
    """
    #solving trait dependency issue
    from ve.fp.complex_util.propagate_distcache import init_propagate_distcache_trait
    init_propagate_distcache_trait(self)

    #paratope and epitope minimum distance
    self.paratope_threshold = 5
    self.epitope_threshold = 5

    #bind the necessary method
    self.find_paratope = MethodType(find_paratope, self)
    self.find_epitope = MethodType(find_epitope, self)

    self.cache_exists = MethodType(cache_exists, self)
    self.load_cache = MethodType(load_cache, self)
    self.save_cache = MethodType(save_cache, self)

    self.paratope_cl_process = MethodType(paratope_cl_process, self)
    self.epitope_cl_process = MethodType(epitope_cl_process, self)

    self.epi_para_cd_process = MethodType(epi_para_cd_process, self)

    
    self.cache_load_methods = {
        "paratope": self.paratope_cl_process,
        "epitope": self.epitope_cl_process,
        }
    self.cache_dump_methods = {
        "paratope": self.epi_para_cd_process,
        "epitope": self.epi_para_cd_process,
        }

    self.paratope = None
    self.epitope = None
    
#output the paratope and epitope into pdb format string

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
    #the antigen path
    source_fp = os.path.join(self.complex_path, self.c_id, "antigen.pdb")
    
    #get the epitope residue ids
    res_ids = map(lambda a:str(a.resnum),self.epitope)
    
    return self.gen_str(source_fp, res_ids)

def paratope_str(self):
    #the antibody path
    source_fp = os.path.join(self.complex_path, self.c_id, "antibody.pdb")
    
    #get the paratope residue ids
    res_ids = map(lambda a:str(a.resnum),self.paratope)
    
    return self.gen_str(source_fp, res_ids)


def write_epitope(self):
    #get the epitope string
    string = self.epitope_str()
    
    #get the epitope path
    output_fp = os.path.join(self.paraepi_path, self.c_id, "epitope.pdb")
    
    #if path does not exist, create it
    if not os.path.exists(os.path.join(self.paraepi_path, self.c_id)):
        os.mkdir(os.path.join(self.paraepi_path, self.c_id))
    
    #write it
    f = open(output_fp, "w")
    f.write(string)
    f.close()

def write_paratope(self):
    #get the paratope string
    string = self.paratope_str()
    
    #get the paratope path
    output_fp = os.path.join(self.paraepi_path, self.c_id, "paratope.pdb")
    
    #if path does not exist, create it
    if not os.path.exists(os.path.join(self.paraepi_path, self.c_id)):
        os.mkdir(os.path.join(self.paraepi_path, self.c_id))
    
    #write it
    f = open(output_fp, "w")
    f.write(string)
    f.close()

from ve.config import data237_complex_root as complex_path, data237_paraepi_root as paraepi_path

def init_io_trait(self,i_path = complex_path, o_path= paraepi_path):
    self.complex_path = i_path
    self.paraepi_path = o_path

    #solve the dependency
    init_find_epiparatope_trait(self)
    
    self.gen_str = MethodType(gen_str, self)
    self.epitope_str = MethodType(epitope_str, self)
    self.paratope_str = MethodType(paratope_str, self)
    self.write_paratope = MethodType(write_paratope, self)
    self.write_epitope = MethodType(write_epitope, self)
        
    
