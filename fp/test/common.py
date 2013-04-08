import unittest

class NumericTestCase(unittest.TestCase):
    def assertArrayEqual(self, first, second):
        self.assertTrue((first == second).all())

    def assertArrayAlmostEqual(self, first, second):
        from itertools import izip
        for a,b in izip(first, second):
            self.assertAlmostEqual(a,b)
        
