from ve.util.load_pdb import load_pdb_struct

from common import *

from ve.fp.complex_util.geom import GeometryTrait
ComplexClass = make_complex_class(GeometryTrait)

class GeomTestCase(NumericTestCase):
    def setUp(self):
        self.c = ComplexClass()

    def test_geom_center(self):
        """test for the geometric center"""
        from ve.fp.geom import Point
        
        #ensure it is a point
        self.assertTrue(isinstance(self.c.get_geom_center(), Point))
        
        #value equality
        self.assertArrayAlmostEqual(self.c.get_geom_center(), [19.38842578, 5.40316406, 13.46964063])


if __name__ == "__main__":
    unittest.main()
