from ve.fp.complex_util.triangle import ResidueTriangleTrait
from ve.fp.fp_80 import Residue as Residue80

from common import *

ResidueClass = make_residue_class(Residue80)

ComplexClass = make_complex_class(ResidueTriangleTrait, residue_class = ResidueClass)

class FindTriangleTest(unittest.TestCase):
    """ Test case for find paratope method"""
    
    def setUp(self):
        self.c = ComplexClass()
        
    def test_get_triangle_with_refresh(self):
        tri = self.c.get_triangles(refresh=True)
        
        expected = 75
        actual = len(tri)
        self.assertEqual(expected, actual)
        
    def test_get_triangle_without_refresh(self):
        """from cache, probably"""
        tri = self.c.get_triangles(refresh=False)
        
        expected = 75
        actual = len(tri)
        self.assertEqual(expected, actual)

        
if __name__ == '__main__':
    unittest.main()