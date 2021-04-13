import unittest

from fuc.api.common import fuc_dir

from fuc.api.VCFResult import VCFResult
from fuc.api.BEDResult import BEDResult
from fuc.api.FASTQResult import FASTQResult
from fuc.api.DataFrame import DataFrame

class TestVCFResult(unittest.TestCase):

    def test_shape(self):
        vcf = VCFResult.read(f'{fuc_dir()}/data/1.vcf')
        self.assertEqual(vcf.shape, (5, 4))

    def test_filter_bed(self):
        vcf1 = VCFResult.read(f'{fuc_dir()}/data/1.vcf')
        bed = BEDResult.read(f'{fuc_dir()}/data/1.bed')
        vcf2 = vcf1.filter_bed(bed)
        self.assertEqual(vcf2.shape, (3, 4))

    def test_merge(self):
        vcf1 = VCFResult.read(f'{fuc_dir()}/data/1.vcf')
        vcf2 = VCFResult.read(f'{fuc_dir()}/data/2.vcf')
        vcf3 = vcf1.merge(vcf2)
        self.assertEqual(vcf3.shape, (7, 5))

    def test_compare(self):
        vcf = VCFResult.read(f'{fuc_dir()}/data/1.vcf')
        result = vcf.compare('Steven', 'Sarah')
        self.assertEqual(result, (0, 1, 0, 4))

class TestBEDResult(unittest.TestCase):

    def test_intersect(self):
        bed1 = BEDResult.read(f'{fuc_dir()}/data/1.bed')
        bed2 = BEDResult.read(f'{fuc_dir()}/data/2.bed')
        bed3 = BEDResult.read(f'{fuc_dir()}/data/3.bed')
        bed4 = bed1.intersect(bed2)
        self.assertEqual(bed3.head, bed4.head)
        self.assertEqual(bed3.data, bed4.data)

class TestFASTQResult(unittest.TestCase):

    def test_shape(self):
        fastq = FASTQResult.read(f'{fuc_dir()}/data/1.fastq')
        self.assertEqual(fastq.shape, 3)

class TestDataFrame(unittest.TestCase):

    def test_merge(self):
        df1 = DataFrame.read(f'{fuc_dir()}/data/left.txt')
        df2 = DataFrame.read(f'{fuc_dir()}/data/right.txt')
        df3 = DataFrame.read(f'{fuc_dir()}/data/merged.txt')
        df4 = df1.merge(df2, ['Name'])
        self.assertEqual(df3.head, df4.head)
        self.assertEqual(df3.data, df4.data)

    def test_get_col(self):
        df = DataFrame.read(f'{fuc_dir()}/data/left.txt')
        col = df.get_col(1)
        self.assertEqual(col, ['30', '25', '41', '28'])

if __name__ == '__main__':
    unittest.main()
