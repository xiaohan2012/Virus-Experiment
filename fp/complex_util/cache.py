import os
from cPickle import load, dump
from ve.config import data237_cache as cache_dir

class CacheTrait(object):
    @classmethod
    def set_signature(cls, t):
        cls.cache_type = t
        
        #make necessary directory creation
        path = os.path.join(cache_dir, cls.cache_type)
        if not os.path.exists(path):
            os.makedirs(path)
    
    @classmethod
    def extract(cls, obj):
        """(CacheTrait, VE related obj) -> Python built-in obj"""
        raise NotImplementedError
        
    @classmethod
    def assemble(cls, obj):
        """(CacheTrait, Python built-in obj) -> VE related obj"""
        raise NotImplementedError
        
    @classmethod
    def load(cls, c_id, c, **kwargs):
        """
        (class, str, Complex) -> obj
        
        load cache according to `c_id`
        the interface to user
        """
        if cls.cache_type is None:
            raise ValueError("Please specify cache type")
        
        return cls.assemble(cls.load_pickle(c_id), c, **kwargs)
        
    @classmethod    
    def dump(cls, c_id, obj):
        """
        (type, str, object) -> NoneType
        dump cache
        the interface to user
        """
        if cls.cache_type is None:
            raise ValueError("Please specify cache type")

        cls.dump_pickle(c_id, cls.extract(obj))
        
    @classmethod
    def has_cache(cls, c_id):
        """
        (Cache, str) -> bool
        
        if the cache exists
        """
        return os.path.exists(cls.get_dir(c_id))
                 
    @classmethod
    def get_dir(cls, c_id):
        return os.path.join(cache_dir, cls.cache_type, c_id)

    @classmethod    
    def load_pickle(cls, c_id):
        """load built-in Python obj"""
        return load(open(cls.get_dir(c_id), "r"))

    @classmethod    
    def dump_pickle(cls, c_id, obj):
        """dump built-in Python obj"""
        return dump(obj, open(cls.get_dir(c_id), "w"))


class TriangleCache(CacheTrait):
    cache_type = "triangle"

    @classmethod
    def extract(cls, triangles):
        """
        (TriangleCache, list of ResTriangle)
        -> triangle tuples: list of tuple(str, str, str)
        """
        return [tuple(tri.res_ids) for tri in triangles]

    @classmethod
    def assemble(cls, resids_lst, c):
        """
        (TriangleCache,
        triangle tuples of the complex: list of tuple(str, str, str),
        the host complex: Complex)
        -> list of ResTriangle
        """
        from ve.fp.residue_util.res_triangle import ResTriangle
        return [ResTriangle(c.get_res_from_resids(resids)) for resids in resids_lst]

from ve.fp.fp import BaseComplexFingerprint, BaseResidueFingerprint

class ComplexFingerPrintCache(CacheTrait):
    """Cache for complex finger print"""
    cache_type = None #please specify
    
        
    @classmethod
    def extract(cls, fp):
        """
        (ComplexFingerPrintCache, BaseComplexFingerPrint) ->  ComplexFingerPrintPickable
        """
        return fp.to_pickable()
        
    @classmethod
    def assemble(cls, p, c, complex_fp_cls = BaseComplexFingerprint, res_fp_cls = BaseResidueFingerprint):
        """
        (ComplexFingerPrintCache, ComplexFingerPrintPickable, Complex) -> BaseComplexFingerPrint
        """
        from ve.fp.fp import BaseComplexFingerprint
        return complex_fp_cls.from_pickable(p, c, res_fp_cls)

class ParaepiCache(CacheTrait):
    """Cache for paraepi"""
    cache_type = None

    @classmethod
    def extract(cls, residues):
        """
        (ParaepiCache, list of Residue) -> (list of str)
        """
        #extract the residue ids in the list of residues
        return map(lambda r: r.res_id, residues)

    @classmethod
    def _assemble_helper(cls, res_ids, obj):
        """(ParaepiCache, list of str, Complex.atg or Complex.atb) -> list of Residue"""
        res_list =  filter(lambda r: r.res_id in res_ids, 
                           obj.residues)
        
        #abnomal if empty list of returned
        assert(len(res_list) > 0)
        
        return res_list

class EpitopeCache(ParaepiCache):
    """Cache for epitope"""
    cache_type = "epitope"

    @classmethod
    def assemble(cls, res_ids, c):
        return cls._assemble_helper(res_ids, c.atg)
        
class ParatopeCache(ParaepiCache):
    """Cache for paratope"""
    cache_type = "paratope"

    @classmethod
    def assemble(cls, res_ids, c):
        return cls._assemble_helper(res_ids, c.atb)
    


