import unittest
from ve.fp.complex_util.split_cylinder import *

class GetIndexTestCase(unittest.TestCase):
    """test for the get_idx method"""

    def setUp(self):
        init_split_cylinder_util(self)

    def test_cylinder_height_range(self):
        self.assertEqual( (-10,10), self.cylinder_height_range)

    def test_layer_size(self):
        self.assertEqual(self.layer_size, 5)
    
    def test_layer_count(self):
        self.assertEqual(self.layer_count, 10)
        
    def test_get_idx_general_case1(self):
        actual = self.get_idx(-9.1 ,4.5)
        expected = 2
        self.assertEqual(actual, expected)
        
    def test_get_idx_general_case2(self):
        actual = self.get_idx(1.2 ,2.5)
        expected = 25 + 1
        self.assertEqual(actual, expected)

    def test_get_idx_general_case3(self):
        actual = self.get_idx(8.2 , 6.1)
        expected = 45 + 3
        self.assertEqual(actual, expected)
    
    def test_get_idx_extreme_case1(self):
        actual = self.get_idx(-9.9 ,0.0)
        expected = 0
        self.assertEqual(actual, expected)

    def test_get_idx_extreme_case2(self):
        actual = self.get_idx(9.9 ,9.9)
        expected = 49
        self.assertEqual(actual, expected)

    def test_get_idx_out_of_range1(self):
        actual = self.get_idx(10 ,9.9)
        expected = None
        self.assertEqual(actual, expected)

    def test_get_idx_out_of_range2(self):
        actual = self.get_idx(-10 ,9.9)
        expected = None
        self.assertEqual(actual, expected)

    def test_get_idx_out_of_range3(self):
        actual = self.get_idx(-9.9 ,10)
        expected = None
        self.assertEqual(actual, expected)

    def test_get_idx_out_of_range4(self):
        actual = self.get_idx(11 ,12)
        expected = None
        self.assertEqual(actual, expected)

class SplitCylinderMethodTestCase(unittest.TestCase):
    """test for the split cylinder method"""
    

if __name__ == "__main__":
    unittest.main()
