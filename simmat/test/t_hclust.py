import unittest

from ve.simmat.source import load_simmat
from ve.simmat.hclust import sim_to_dist

class SimToDistTestCase(unittest.TestCase):
    """Sim mat to  dist mat test case"""
    def test_general_case(self):
        sim_mat = load_simmat("data/fp_370_atg.txt")
        dist_mat = sim_to_dist(sim_mat)
        self.assertEqual(1- sim_mat[1,0], dist_mat[1,0])


if __name__ == "__main__":
    unittest.main()
