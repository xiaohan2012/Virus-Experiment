from schrodinger.structure import Structure

class mystructure(object):
    def __init__(self,structure):
        self.st = structure

    def to_json(self):
        from simplejson import dumps
        return dumps({
            "atoms": [{"element": a.element,"location": a.xyz} for a in self.atom],
            "bonds":[{"atoms": (b.atom1.index,b.atom2.index), "order": b.order} for b in self.bond]
        })
        
    def __getattr__(self,name):
        if hasattr(self.st, name):
            return getattr(self.st,name)
        else:
            raise AttributeError("\"%s\" method not found." %name)

if __name__ == '__main__':
    from ve.util.load_pdb import load_pdb_struct
    c = load_pdb_struct("data/data480/ab_pdb/1SLG_B.pdb")
    print c.to_json()