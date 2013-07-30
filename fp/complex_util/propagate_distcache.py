"""
Propagate the complex's distance matrix cache to all its residues
"""
from itertools import chain
from types import MethodType

from ve.util.dist import ResDistCache

from common import make_logger
logger = make_logger("Distance Cache Propagation")

class DistanceCachePropagationTrait(object):
    """
    Distance cache needs to be populated to residues
    propagate the distance matrix cache to all the subordinating residues
    """
    def __init__(self, cache_cls = ResDistCache):
        super(DistanceCachePropagationTrait,self).__init__()
        self.distcache = cache_cls()

        self.prop_dist_cache()
        
    def prop_dist_cache(self):
        """
        propagate the distance matrix cache to all the subordinating residues
        """
        logger.debug("propagating distance cache to residues")
        for r in chain(self.atg.residues, self.atb.residues):
            r.distcache = self.distcache
"""
def init_propagate_distcache_trait(self, cache_cls = ResDistCache):

    logger.info("initializing distance cache propagation trait")

    self.prop_dist_cache = MethodType(prop_dist_cache, self)
    
    self.distcache = cache_cls()

def prop_dist_cache(self):

    logger.debug("propagating distance cache to residues")
    for r in chain(self.atg.residues, self.atb.residues):
        r.distcache = self.distcache

"""
