import unittest

import os
import numpy as np

from ve.simmat.source import make_dataloader, make_single_line_converter

single_line_headless_converter = make_single_line_converter(slice(0,370))

class FPStringConverterTestCase(unittest.TestCase):
    
    def test_single_line_headless_converter(self):
        string = "0.00,0.00,0.00,2.00,2.00,5.00,4.00,5.00,6.00,3.00,7.00,2.00,0.00,2.00,2.00,2.00,3.00,1.00,4.00,0.00,2.00,2.00,0.00,1.00,3.00,1.00,1.00,7.00,2.00,2.00,4.00,2.00,1.00,0.00,2.00,2.00,3.00,4.00,1.00,3.00,6.00,6.00,1.00,1.00,4.00,2.00,2.00,1.00,3.00,2.00,9.00,2.00,0.00,2.00,3.00,4.00,1.00,5.00,6.00,4.00,4.00,3.00,0.00,0.00,2.00,1.00,1.00,2.00,2.00,3.00,3.00,2.00,0.00,1.95,3.12,8.08,9.67,9.10,22.65,17.73,23.42,26.66,0.00,0.03,0.33,0.10,0.47,0.31,0.45,1.00,0.98,1.82,0.00,1.00,7.00,4.00,8.00,12.00,13.00,25.00,23.00,34.00"
        
        actual = single_line_headless_converter(string)
        expected = np.array([0.00,0.00,0.00,2.00,2.00,5.00,4.00,5.00,6.00,3.00,7.00,2.00,0.00,2.00,2.00,2.00,3.00,1.00,4.00,0.00,2.00,2.00,0.00,1.00,3.00,1.00,1.00,7.00,2.00,2.00,4.00,2.00,1.00,0.00,2.00,2.00,3.00,4.00,1.00,3.00,6.00,6.00,1.00,1.00,4.00,2.00,2.00,1.00,3.00,2.00,9.00,2.00,0.00,2.00,3.00,4.00,1.00,5.00,6.00,4.00,4.00,3.00,0.00,0.00,2.00,1.00,1.00,2.00,2.00,3.00,3.00,2.00,0.00,1.95,3.12,8.08,9.67,9.10,22.65,17.73,23.42,26.66,0.00,0.03,0.33,0.10,0.47,0.31,0.45,1.00,0.98,1.82,0.00,1.00,7.00,4.00,8.00,12.00,13.00,25.00,23.00,34.00])

        self.assertTrue((expected == actual).all())

class DataloadTestCase(unittest.TestCase):
    def test_single_line_headless_loader(self):
        """ test on the single line headess 370-bits finger print data loader"""
        from ve.config import data237_root

        fp_dir = os.path.join(data237_root, "fp_370_atg")
        dataloader = make_dataloader(fp_dir, single_line_headless_converter)
        
        actual = dataloader("4DKF_A")
        expected = np.array([0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.00,0.00,1.00,3.00,2.00,1.00,2.00,0.00,1.00,0.00,0.00,1.00,3.00,2.00,3.00,5.00,2.00,1.00,1.00,2.00,3.00,3.00,3.00,5.00,1.00,2.00,2.00,4.00,0.00,0.00,0.00,0.60,1.93,6.16,8.17,19.60,14.20,7.80,0.00,0.00,0.00,0.04,0.00,0.30,-0.01,0.26,0.47,0.62,0.00,0.00,0.00,4.00,0.00,15.00,2.00,13.00,7.00,8.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,2.00,1.00,2.00,4.00,4.00,7.00,4.00,4.00,0.00,1.00,0.00,3.00,3.00,4.00,3.00,6.00,7.00,5.00,3.00,0.00,3.00,1.00,2.00,4.00,3.00,5.00,3.00,9.00,7.00,0.00,2.00,1.00,2.00,4.00,5.00,3.00,2.00,5.00,3.00,0.00,0.00,0.00,1.00,2.00,0.00,2.00,0.00,1.00,1.00,0.00,0.00,0.00,2.00,1.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,2.34,2.02,0.60,9.79,10.82,12.54,22.44,14.71,21.08,0.00,0.18,0.00,0.15,0.52,0.44,0.71,0.31,0.68,0.96,0.00,2.00,0.00,1.00,8.00,8.00,13.00,5.00,8.00,14.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,1.00,0.00,0.00,1.00,1.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,2.00,0.00,0.00,2.00,2.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,8.00,0.00,0.00,8.00,7.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,5.00,0.00,0.00,5.00,4.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,14.00,0.00,0.00,14.00,9.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,9.00,0.00,0.00,9.00,8.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,7.00,0.00,0.00,7.00,4.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00])
        
        self.assertTrue((actual == expected).all())

if __name__ == "__main__":
    unittest.main()
