from schrodinger.structure import StructureReader

def load_pdb_struct(path):
    return StructureReader(path).next()
