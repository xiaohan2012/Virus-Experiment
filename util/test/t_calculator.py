import unittest

from ve.util.calculator import complex_pairwise_calc

class PairWiseCalculatorTestCase(unittest.TestCase):
    """Test case for pair wise calculator"""
    
    def test_silly_asymmetrical_case(self):
        """very silly and asymmetrical test case"""
        complex_id_list = ["a", "b", "c"]
        dataloader = lambda s: s
        calc_func = lambda c1,c2: c1 + c2
        
        actual,_ = complex_pairwise_calc(complex_id_list, dataloader, calc_func, callback = None, symmetry=False)
        
        expected = [("a", "b", "ab"), ('a', 'c', 'ac'), ('b', 'a', 'ba'), ('b', 'c', 'bc'), ('c', 'a', 'ca'), ('c', 'b', 'cb'),
                    ('a','a','aa'), ('b','b','bb'), ('c','c','cc')]
        
        self.assertEqual(actual, expected)

    def test_silly_symmetrical_case(self):
        """very silly and symmetrical test case"""
        complex_id_list = range(1,4) #1~3
        dataloader = lambda d: d
        calc_func = lambda d1,d2: d1 + d2
        
        actual,_ = complex_pairwise_calc(complex_id_list, dataloader, calc_func, callback = None, symmetry=True)
        
        expected = [(1, 2, 3), (1, 3, 4), (2, 3, 5),
                    (1, 1, 2), (2, 2, 4), (3, 3, 6)]
        
        self.assertEqual(actual, expected)

    def test_callback_with_return_values(self):
        """similar to silly and symmetrical test case,
        but with calback
        """
        complex_id_list = range(1,4) #1~3
        dataloader = lambda d: d
        calc_func = lambda d1,d2: d1 + d2
        callback = lambda c1,c2,v: (c1,c2,v+1)

        actual,_ = complex_pairwise_calc(complex_id_list, dataloader, calc_func, callback = callback, symmetry=True)
        
        expected = [(1, 2, 4), (1, 3, 5), (2, 3, 6), 
                    (1, 1, 3), (2, 2, 5), (3, 3, 7)]
        
        self.assertEqual(actual, expected)
        
    def test_exception_list_empty(self):
        """test for the case no error are encountered,that is to say the exception list is empty"""
        complex_id_list = range(1,4) #1~3
        dataloader = lambda d: d
        calc_func = lambda d1,d2: d1 + d2
        
        _, actual = complex_pairwise_calc(complex_id_list, dataloader, calc_func, symmetry=False)
        expected = []

        self.assertEqual(actual, expected)
        
    def test_exception_list_nonempty(self):
        """test for the case errors are encountered,that is to say the exception list is NOT empty"""
        complex_id_list = map(float, range(0,3)) #0.0 ~ 2.0
        dataloader = lambda d: d
        calc_func = lambda d1,d2: d1 / d2 #DivisionByZeroError
        
        actual_result, actual_exception = complex_pairwise_calc(complex_id_list, dataloader, calc_func, symmetry=False)
        expected_exception = [(1,0),(2,0),(0,0)]
        expected_result = [(0.0, 1.0, 0.0), (0.0, 2.0, 0.0), (1.0, 2.0, 0.5), (2.0, 1.0, 2.0), (1.0, 1.0, 1.0), (2.0, 2.0, 1.0)]
        
        self.assertEqual(actual_result, expected_result)
        self.assertEqual(actual_exception, expected_exception)

if __name__ == "__main__":
    unittest.main()
