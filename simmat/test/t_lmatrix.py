import unittest

import numpy as np

from ve.simmat.lmatrix import lmatrix

class LMatrixTestCase1(unittest.TestCase):
    """Testcase for labeled matrx
    in this case, the matrix data is not set
    """

    def setUp(self):
        labels = ["a","b","c"]
        self.m = lmatrix(labels)        
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


    def test_to_table(self):
        """test for table formulation"""
        actual = self.m.to_csv_str()
        expected = """,a,b,c
a,1.000000,2.000000,3.000000
b,4.000000,5.000000,6.000000
c,7.000000,8.000000,9.000000"""
        self.assertEqual(actual, expected)
        

class LMatrixTestCase2(unittest.TestCase):
    """Testcase for labeled matrx
    in this case, the matrix data is set
    """

    def setUp(self):
        labels = ["a","b","c"]
        self.m = lmatrix(labels, data=np.array([[1,2,3],[4,5,6],[7,8,9]]))
        
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

class LoadSimmatTestCase(unittest.TestCase):
    """similarity matrix loading test case"""

    def setUp(self):
        """load the matrix"""
        from ve.simmat.source import load_simmat
        self.m = load_simmat("data/fp_370_atg.txt")

    def test_content_matching_1(self):
        c1 = '2NLJ_C'
        c2 = '1OTT_A'
        actual = self.m[c1,c2]
        expected = 0.888093103028
        self.assertEqual(actual, expected)
    

class FromDBTest(unittest.TestCase):
    def test_dimension_and_val(self):
        from ve.util.load_pdb import complex_ids
        from ve.config import epi166_fp
        from ve.dbconfig import db
        
        cids = complex_ids(epi166_fp)

        mat = lmatrix.from_db(db["epi_166"], cids)

        self.assertEqual((166, 166), mat.shape)
        
        self.assertAlmostEqual(mat["2VYR_A", "2VYR_A"], 26 + 260 +139)


if __name__ == "__main__":
    unittest.main()
