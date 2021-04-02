import unittest
import pathlib

from VCFResult import VCFResult
from BEDResult import BEDResult
from DataFrame import DataFrame

FUC_DIR = pathlib.Path(__file__).parent.absolute()

class TestVCFResult(unittest.TestCase):

    def test_shape(self):
        vcf = VCFResult.read(f'{FUC_DIR}/data/1.vcf')
        self.assertEqual(vcf.shape, (5, 4))

    def test_filter_bed(self):
        vcf1 = VCFResult.read(f'{FUC_DIR}/data/1.vcf')
        bed = BEDResult.read(f'{FUC_DIR}/data/1.bed')
        vcf2 = vcf1.filter_bed(bed)
        self.assertEqual(vcf2.shape, (3, 4))

    def test_merge(self):
        vcf1 = VCFResult.read(f'{FUC_DIR}/data/1.vcf')
        vcf2 = VCFResult.read(f'{FUC_DIR}/data/2.vcf')
        vcf3 = vcf1.merge(vcf2)
        self.assertEqual(vcf3.shape, (7, 5))

class TestBEDResult(unittest.TestCase):

    def test_intersect(self):
        bed1 = BEDResult.read(f'{FUC_DIR}/data/1.bed')
        bed2 = BEDResult.read(f'{FUC_DIR}/data/2.bed')
        bed3 = BEDResult.read(f'{FUC_DIR}/data/3.bed')
        bed4 = bed1.intersect(bed2)
        self.assertEqual(bed3.head, bed4.head)
        self.assertEqual(bed3.data, bed4.data)

class TestDataFrame(unittest.TestCase):

    def test_merge(self):
        df1 = DataFrame.read(f'{FUC_DIR}/data/left.txt')
        df2 = DataFrame.read(f'{FUC_DIR}/data/right.txt')
        df3 = DataFrame.read(f'{FUC_DIR}/data/merged.txt')
        df4 = df1.merge(df2, 'Name')
        self.assertEqual(df3.head, df4.head)
        self.assertEqual(df3.data, df4.data)

if __name__ == '__main__':
    unittest.main()
