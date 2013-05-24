from common import *

from ve.fp.residue_util.geom import GeometryTrait
Residue = make_residue_class(GeometryTrait)

from ve.fp.complex_util.res_spat_dist import ResidueSpatialDistributionTrait
Complex = make_complex_class(ResidueSpatialDistributionTrait, residue_class = Residue)


class ResidueDistributionTestCase(NumericTestCase):
    """Test case for residue spatial distribution within the cylinder"""
    
    def setUp(self):
        self.c = Complex()

    def test_atg_side_count(self):
        """Residue count distribution on the antigen side"""
        dist, _ = self.c.get_atg_res_spat_dist()
        actual = dist[63]
        expected = 4
        self.assertEqual(actual, expected)

    def test_atb_side_count(self):
        """Residue count distribution on the antibody side"""
        dist, _ = self.c.get_atb_res_spat_dist()
        actual = dist[33]
        expected = 1
        self.assertEqual(actual, expected)

    def test_atg_side_residue(self):
        """Residue distribution on the antigen side"""
        _, dist = self.c.get_atg_res_spat_dist()
        actual = len(dist[63])
        expected = 4
        self.assertEqual(actual, expected)

    def test_atb_side_residue(self):
        """Residue distribution on the antibody side"""
        _, dist = self.c.get_atb_res_spat_dist()
        actual = len(dist[33])
        expected = 1
        self.assertEqual(actual, expected)

        
if __name__ == "__main__":
    unittest.main()
        

