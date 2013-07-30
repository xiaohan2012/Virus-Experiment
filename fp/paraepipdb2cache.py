"""
Convert paratope & epitope pdb files to corresponding res_ids-formated cache files
"""

from ve.fp.complex_util.cache import ParatopeCache as PCache, EpitopeCache as ECache

def convert(c_id, st, CacheClass):
    """(str, pdb structure, Cache) -> None"""
    CacheClass.dump(c_id, st.residues)

def main():
    import os
    from ve.util.load_pdb import complex_ids, load_pdb_struct
    from ve.config import data237_root

    paraepi_dir = os.path.join(data237_root, "paraepi")
    #cids = complex_ids(paraepi_dir)
    cids = ["1SLG_D"]
    for cid in cids:
        paratope_st = load_pdb_struct(os.path.join(paraepi_dir, cid, "paratope.pdb"))
        convert(cid, paratope_st, PCache)
        
        epitope_st = load_pdb_struct(os.path.join(paraepi_dir, cid, "epitope.pdb"))
        convert(cid, epitope_st, ECache)

if __name__ == '__main__':
    main()
        
    