from ve.util.residue import BaseResidue

class MyResidue(BaseResidue):
    def __init__(self, res):
        BaseResidue.__init__(self, res)
        
        from ve.fp.residue_util.res_geom import init_geom_trait
        init_geom_trait(self)
        

from ve.util.complex import BaseComplex        

class MyComplex(BaseComplex):
    def __init__(self,c_id, atg, atb):
        self.c_id = c_id
        self.atg = atg
        self.atb = atb
        
        #Split cylinde part finger print
        from ve.fp.complex_util.split_cylinder import init_gcb_split_cylinder_trait
        init_gcb_split_cylinder_trait(self)
        
        #physical chemical finger print
        from ve.fp.complex_util.physi_chemi import init_physi_chemi_fp_trait
        init_physi_chemi_fp_trait(self)
        
        #interactive force finger print
        from ve.fp.complex_util.interactive_force_fp import init_interative_force_fp_trait
        init_interative_force_fp_trait(self)
        
        
    def get_fp(self):
        #antigen side split cylinder fingerprint
        fp1 = self.get_atg_fp_by_split_cylinder()

        #antigen side physical chemical fingerprint
        fp2 = self.get_physi_chemi_atg_fp()

        #antibody side split cylinder fingerprint
        fp3 = self.get_atb_fp_by_split_cylinder()

        #antigen side physical chemical fingerprint
        fp4 = self.get_physi_chemi_atb_fp()
        
        #interactive force finger print 
        #antigen as the receptor
        fp5 = self.gen_if_complex_fp_atg()
        
        return fp1.append(fp2).append(fp3).append(fp4).append(fp5)



