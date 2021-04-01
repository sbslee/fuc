import unittest
import pathlib

from VCFResult import VCFResult
from BEDResult import BEDResult

FUC_DIR = pathlib.Path(__file__).parent.absolute()

class TestVCFResult(unittest.TestCase):

    def test_shape(self):
        vcf = VCFResult.read(f'{FUC_DIR}/data/1.vcf')
        self.assertEqual(vcf.shape, (4, 5))

class TestBEDResult(unittest.TestCase):

    def test_intersect(self):
        bed1 = BEDResult.read(f'{FUC_DIR}/data/1.bed')
        bed2 = BEDResult.read(f'{FUC_DIR}/data/2.bed')
        bed3 = bed1.intersect(bed2)
        self.assertEqual(bed3.data, BEDResult.read(f'{FUC_DIR}/data/3.bed').data)

if __name__ == '__main__':
    unittest.main()
