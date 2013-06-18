import os
from common import *

class TriangleCacheTest(unittest.TestCase):
    """ Test case for find paratope method"""
    
    def setUp(self):
        from ve.fp.complex_util.cache import TriangleCache as C
        self.C = C
        self.C.set_signature("triangle")

        from ve.fp.complex_util.triangle import ResidueTriangleTrait
        from ve.fp.fp_80 import Residue as Residue80
        ResidueClass = make_residue_class(Residue80)
        Complex = make_complex_class(ResidueTriangleTrait, residue_class = ResidueClass)
        
        self.c = Complex()
        
    def test_dump(self):
        tri = self.c.get_triangles()
        self.C.dump(self.c.c_id, tri)

        path = self.C.get_dir(self.c.c_id)
        self.assertTrue(os.path.exists(path))
        self.assertTrue(self.C.has_cache(self.c.c_id))

    def test_load(self):
        tri = self.c.get_triangles()
        self.C.dump(self.c.c_id, tri)

        new_tri = self.C.load(self.c.c_id, self.c)

        self.assertEqual(tri, new_tri)
        
    def tearDown(self):
        """don't forget the clean the caches'"""
        path = self.C.get_dir(self.c.c_id)
        os.remove(path)
        

class ComplexFingerPrintCacheTest(unittest.TestCase):
    """ Test case for complex fp cache"""
    
    def setUp(self):
        #cache init
        from ve.fp.complex_util.cache import ComplexFingerPrintCache as C
        self.C = C
        self.C.set_signature("complex_fp_test")

        #complex class init
        from ve.fp.fp_80 import Residue as Residue80
        from ve.fp.complex_util.split_cylinder import SplitCylinderViaComplexPlaneTrait
        Residue = make_residue_class(Residue80)
        Complex = make_complex_class(SplitCylinderViaComplexPlaneTrait, Residue)

        self.c = Complex()

    def test_dump(self):
        fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                            targets=[(0,self.c.atb.residues)])
        self.C.dump(self.c.c_id, fps)

        path = self.C.get_dir(self.c.c_id)
        self.assertTrue(os.path.exists(path))
        self.assertTrue(self.C.has_cache(self.c.c_id))

    def test_load(self):
        #first dump it
        fps = self.c.gen_fp_by_splitting_cylinder(bases=self.c.atb.residues,
                                            targets=[(0,self.c.atb.residues)])
        self.C.dump(self.c.c_id, fps)

        #then load it
        new_fps = self.C.load(self.c.c_id, self.c)

        self.assertEqual(fps, new_fps)

    def test_padded_complex_fp(self):
        pass
        
    def tearDown(self):
        """don't forget the clean the caches'"""
        path = self.C.get_dir(self.c.c_id)
        os.remove(path)
        
if __name__ == '__main__':
    unittest.main()