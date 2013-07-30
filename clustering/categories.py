from lmatrix import LMatrix

data = {
    "lysozyme": ['1MLC_F', '1YQV_Y', '1P2C_C', '1JHL_A', '1OP9_B', '1FBI_X', '1J1O_Y', '3D9A_C', '1DZB_X', '2ZNX_Y', '1G7J_C', '1FDL_Y', '1A2Y_C', '1BVK_C', '1SQ2_L', '2I25_M', '1RI8_B', '1RJC_B', '1JTP_L'],
    "influenza_A_HA": ['1EO8_A', '3GBN_B', '3LZF_A'],
    "influenza_A_NA": ['1NCA_N', '1NMB_N', '2AEP_A'],
    "HIV_envelope_glycoprotein": ['1YYM_P', '2NY1_A', '2NY7_G', '3NGB_A', '2XQY_A', '2R29_A'],
    "HIV_transmembrane_glycoprotein": ['3MA9_A', '3MAC_A', '2CMR_A', '2XRA_A']
}


def categorize(mat, categories):
    """(labeled matrix, dict(str -> list)) => list of labeled matrix

    >>> from load_data import load_mat
    >>> mats = categorize(load_mat("data/matt.csv"), data)
    >>> len(mats)
    5
    """
    mats = {}
    for cat_name, cid_list in categories.items():
        labels = cid_list

        sub_mat = LMatrix(labels)
        for cid1 in cid_list:
            for cid2 in cid_list:
                sub_mat[cid1, cid2] = mat[cid1, cid2]
        mats[cat_name] = sub_mat
        
    return mats
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()