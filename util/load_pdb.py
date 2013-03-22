import os
from glob import glob

from schrodinger.structure import StructureReader
from structure import mystructure

def load_pdb_struct(path,residue_cls  = None):
    st = StructureReader(path).next()
    st = mystructure(st)
    if residue_cls is not None:
        residues = map(residue_cls,st.residue)
        #use custom Residue class
        st.residues = residues
    return st

def complex_ids(path):
    return map(lambda s: s.split(".")[0],
               map(os.path.basename, 
                   glob(os.path.join(path ,"*"))))
        
