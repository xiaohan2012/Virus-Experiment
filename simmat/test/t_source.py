import unittest

from ve.simmat import source as s

cid = "0HEZ_E"

class DataLoaderTestCase(unittest.TestCase):
    """Testcase for dataloader"""
    
    def test_fp370_atb_dataloader_fp_length(self):
        """fp370_atb_dataloader, length should be 370"""
        actual = len(s.fp370_atb_dataloader(cid))
        expected = 370
        self.assertEqual(actual, expected)

    def test_fp370_atg_dataloader_fp_length(self):
        """fp370_atg_dataloader, length should be 370"""
        actual = len(s.fp370_atg_dataloader(cid))
        expected = 370
        self.assertEqual(actual, expected)

    def test_first_110_atg_dataloader_fp_length(self):
        """first_110_atg_dataloader, length should be 110"""
        actual = len(s.first_110_atg_dataloader(cid))
        expected = 110
        self.assertEqual(actual, expected)

    def test_second_110_atg_dataloader_fp_length(self):
        """second_110_atg_dataloader, length should be 110"""
        actual = len(s.second_110_atg_dataloader(cid))
        expected = 110
        self.assertEqual(actual, expected)

    def test_last_150_atg_dataloader_fp_length(self):
        """last_150_atg_dataloader, length should be 110"""
        actual = len(s.last_150_atg_dataloader(cid))
        expected = 150
        self.assertEqual(actual, expected)


    def test_first_110_atb_dataloader_fp_length(self):
        """first_110_atb_dataloader, length should be 110"""
        actual = len(s.first_110_atb_dataloader(cid))
        expected = 110
        self.assertEqual(actual, expected)

    def test_second_110_atb_dataloader_fp_length(self):
        """second_110_atb_dataloader, length should be 110"""
        actual = len(s.second_110_atb_dataloader(cid))
        expected = 110
        self.assertEqual(actual, expected)

    def test_last_150_atb_dataloader_fp_length(self):
        """last_150_atb_dataloader, length should be 110"""
        actual = len(s.last_150_atb_dataloader(cid))
        expected = 150
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
