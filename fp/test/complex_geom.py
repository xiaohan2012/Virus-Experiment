from ve.util.load_pdb import load_pdb_struct

from common import *


class GeomTestCase(NumericTestCase):
    def setUp(self):

        self.atg = load_pdb_struct(os.path.join(test_data_dir, "antigen.pdb"), TestResidue)
        self.atb = load_pdb_struct(os.path.join(test_data_dir, "antibody.pdb"), TestResidue)

        self.c_id = self.c_id = "1SLG_D"


        from ve.fp.complex_util.geom import init_geom_trait
        init_geom_trait(self)

    def test_geom_center(self):
        """test for the geometric center"""
        from ve.fp.geom import Point
        
        #ensure it is a point
        self.assertTrue(isinstance(self.get_geom_center(), Point))
        
        #value equality
        self.assertArrayAlmostEqual(self.get_geom_center(), [19.38842578, 5.40316406, 13.46964063])


if __name__ == "__main__":
    unittest.main()
