import unittest

import numpy as np

from ve.simmat.sim_metric import make_single_line_fp_similatirty_calculator, corr_coef

class SingleLineFPCalculatorTestCase(unittest.TestCase):
    def test_vector_add_sum_case(self):
        """vector addition on the two fp arrays and sum up the result array"""

        metric_func = lambda fp1,fp2: np.sum(fp1 + fp2)

        fp1 = np.array([1,0,1,0,1,0])
        fp2 = np.array([0,1,0,1,0,1])
        
        sim_calc = make_single_line_fp_similatirty_calculator(metric_func)

        actual = sim_calc(fp1, fp2)
        expected = 6
        
        self.assertEqual(actual, expected)

    def test_corrcoef_calc(self):
        """test for the correlation coefficient calculator"""
        
        fp1 = np.array([1,2,3])
        fp2 = np.array([3,2,1])
        
        sim_calc = make_single_line_fp_similatirty_calculator(corr_coef)
        
        actual = sim_calc(fp1, fp2)
        expected = -1

        self.assertEqual(expected, actual)

if __name__ == "__main__":
    unittest.main()
        
