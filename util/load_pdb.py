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
        st.src = path

    return st

from ve.config import data237_complex_root as complex_dir

def complex_ids(path=complex_dir):
    """
    (str) => list of str

    get complex ids
    
    >>> ids = complex_ids()
    >>> len(ids)
    456
    """
    return sorted(map(lambda s: s.split(".")[0],
               map(os.path.basename, 
                   glob(os.path.join(path ,"*")))))
    
from ve.util.complex import BaseComplex

def load_complex(cid, atg_path, atb_path, complex_cls=BaseComplex, residue_cls=BaseResidue):
    atg = load_pdb_struct(atg_path, residue_cls)

    atb = load_pdb_struct(atb_path, residue_cls)
        
    c = complex_cls(complex_id = cid, antigen = atg, antibody = atb)

    return c
    
def load_complexes(complex_ids, directory = complex_dir, complex_cls=BaseComplex, residue_cls=BaseResidue):
    """
    load complexes in a row

    >>> from ve.config import data480_complex_root
    >>> ids = complex_ids(data480_complex_root)
    >>> n=10
    >>> cs = load_complexes(ids[:n], directory = data480_complex_root)
    >>> cs=list(cs)
    >>> print cs[0].c_id == ids[0]
    True
    >>> len(list(cs))
    10
    """
    for cid in complex_ids:
        yield load_complex(cid, os.path.join(directory, cid, "antigen.pdb"), os.path.join(directory, cid, "antibody.pdb"), complex_cls, residue_cls)

    
def test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    test()
