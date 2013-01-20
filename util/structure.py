from schrodinger.structure import Structure

class mystructure(object):
    def __init__(self,structure):
        self.st = structure

    def __getattr__(self,name):
       if hasattr(self.st, name):
            return getattr(self.st,name)
       else:
            raise AttributeError("\"%s\" method not found." %name)
