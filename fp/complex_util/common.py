from ve.util.logger import make_logger

class NoAtgNorAtbLoadedError(Exception):pass

def is_atg_atb_loaded(self):
    return hasattr(self, "atg") and hasattr(self,"atb")

