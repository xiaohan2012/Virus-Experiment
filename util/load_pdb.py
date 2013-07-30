import os
from glob import glob

from schrodinger.structure import StructureReader
from structure import mystructure
from ve.util.residue import BaseResidue

def load_pdb_struct(path,residue_cls  = BaseResidue):
    st = StructureReader(path).next()

    st = mystructure(st)
    if residue_cls is not None:
        residues = map(residue_cls,st.residue)
        #use custom Residue class
        st.residues = residues

    return st

from ve.config import data237_complex_root as complex_dir

def complex_ids(path=complex_dir):
    """
    (str) => list of str

    get complex ids
    
    >>> ids = complex_ids()
    >>> len(ids)
    236
    """
    return sorted(map(lambda s: s.split(".")[0],
               map(os.path.basename, 
                   glob(os.path.join(path ,"*")))))
    
from ve.util.complex import BaseComplex


def load_complexes(complex_ids, directory = complex_dir, complex_cls=BaseComplex, residue_cls=BaseResidue):
    """
    load complexes in a row

    >>> ids = complex_ids()
    >>> cs = load_complexes(ids)
    >>> len(list(cs))
    236
    >>> from ve.config import data480_complex_root
    >>> ids = complex_ids(data480_complex_root)
    >>> cs = load_complexes(ids, directory = data480_complex_root)
    >>> len(list(cs))
    479
    """
    for cid in complex_ids:        
        atg = load_pdb_struct(os.path.join(directory, cid, "antigen.pdb"), residue_cls)

        atb = load_pdb_struct(os.path.join(directory, cid, "antibody.pdb"), residue_cls)
        
        c = complex_cls(complex_id = cid, antigen = atg, antibody = atb)

        yield c

    
def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    test()
