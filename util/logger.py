"""
logger related utility
"""

import sys
import logging

def make_logger(name):
    logging.basicConfig( stream=sys.stderr )
    logging.getLogger(name).setLevel( logging.DEBUG )
    return logging.getLogger(name)
    
