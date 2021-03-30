import unittest

from BEDResult import BEDResult

class TestBEDResult(unittest.TestCase):

    def test_intersect(self):
        bed1 = BEDResult.read('data/1.bed')
        bed2 = BEDResult.read('data/2.bed')
        bed3 = bed1.intersect(bed2)
        self.assertEqual(bed3.data, BEDResult.read('data/3.bed').data)

if __name__ == '__main__':
    unittest.main()
