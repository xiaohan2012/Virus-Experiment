from cPickle import load, dump
from ve.config import data237_cache as cache_dir

class CacheTrait(object):
    @classmethod
    def set_signature(cls, t):
        cls.cache_type = t
        
        #make necessary directory creation
        import os
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
    def load(cls, c_id, c):
        """
        (class, str, Complex) -> obj
        
        load cache according to `c_id`
        the interface to user
        """
        if cls.cache_type is None:
            raise ValueError("Please specify cache type")
        
        return cls.assemble(cls.load_pickle(c_id), c)
        
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
        import os
        return os.path.exists(cls.get_dir(c_id))
                 
    @classmethod
    def get_dir(cls, c_id):
        import os
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
    def assemble(cls, p, c):
        """
        (ComplexFingerPrintCache, ComplexFingerPrintPickable, Complex) -> BaseComplexFingerPrint
        """
        from ve.fp.fp import BaseComplexFingerprint
        return BaseComplexFingerprint.from_pickable(p, c)
    