import unittest

import os
import numpy as np

from ve.simmat.lmatrix import lmatrix

class LMatrixTestCase(unittest.TestCase):
    """Testcase for labeled matrx"""

    def setUp(self):
        self.m = lmatrix(["a","b","c"])        
        self.m[0,:] = [1,2,3]
        self.m[1,:] = [4,5,6]
        self.m[2,:] = [7,8,9]

    def test_one_item_getting(self):
        """getting one item at a item"""
        actual = self.m["a", "b"]
        expected = 2
        self.assertEqual(actual, expected)

    def test_one_item_setting(self):
        """setting one item at a item"""
        self.m["c","c"] = -9
        expected = -9
        actual = self.m[2,2]
        self.assertEqual(expected, actual)

    def test_row_slicing_getting(self):
        """getting a row using slicing"""
        actual = self.m["a",:]
        expected = np.array([1,2,3])
        self.assertTrue( (actual == expected).all())

    def test_row_slicing_setting(self):
        """setting a row using slicing"""
        self.m["b",:] = [6,5,4]
        actual = self.m[1,:]
        expected = np.array([6,5,4])
        self.assertTrue( (actual == expected).all())

    def test_column_slicing_getting(self):
        """getting a column using slicing"""
        actual = self.m[:,"a"]
        expected = np.array([1,4,7])
        self.assertTrue( (actual == expected).all())

    def test_column_slicing_setting(self):
        """setting a column using slicing"""
        self.m[:,"c"] = [9,6,3]
        actual = self.m[:,2]
        expected = np.array([9,6,3])
        self.assertTrue( (actual == expected).all())

class LoadSimmatTestCase(unittest.TestCase):
    """similarity matrix loading test case"""

    def setUp(self):
        """load the matrix"""
        from ve.simmat.io import load_simmat
        self.m = load_simmat("fp_370_atg.txt")

    def test_content_matching_1(self):
        c1 = '2NLJ_C'
        c2 = '1OTT_A'
        actual = self.m[c1,c2]
        expected = 0.888093103028
        self.assertEqual(actual, expected)
    

if __name__ == "__main__":
    unittest.main()
