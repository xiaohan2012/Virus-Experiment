"""
Propagate the complex's distance matrix cache to all its residues
"""
from itertools import chain
from types import MethodType

from ve.util.dist import ResDistCache

import sys
import logging

logging.basicConfig( stream=sys.stderr )
logging.getLogger("Distance Cache Propagation").setLevel( logging.DEBUG )
logger = logging.getLogger("Distance Cache Propagation")


def init_propagate_distcache_trait(self, cache_cls = ResDistCache):
    """ inject this util to the complex class"""
    logger.info("initializing distance cache propagation trait")

    self.prop_dist_cache = MethodType(prop_dist_cache, self)
    
    self.distcache = cache_cls()

def prop_dist_cache(self):
    """
    propagate the distance matrix cache to all the subordinating residues
    """
    logger.debug("propagating distance cache to residues")
    for r in chain(self.atg.residues, self.atb.residues):
        r.distcache = self.distcache


