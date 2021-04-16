import unittest

from fuc.api.common import fuc_dir
from fuc.api.VcfFrame import VcfFrame
from fuc.api.BedFrame import BedFrame
from fuc.api.FastqFrame import FastqFrame
from fuc.api.DataFrame import DataFrame

class TestVcfFrame(unittest.TestCase):

    def test_shape(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        self.assertEqual(vf.shape, (5, 4))

    def test_filter_bed(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        bed = BedFrame.from_file(f'{fuc_dir()}/data/bed/1.bed')
        vf = vf.filter_bed(bed)
        self.assertEqual(vf.shape, (3, 4))

    def test_merge(self):
        vf1 = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf2 = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/2.vcf')
        vf3 = vf1.merge(vf2)
        self.assertEqual(vf3.shape, (7, 5))

    def test_compare(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        self.assertEqual(vf.compare('Steven', 'Sarah'), (0, 1, 1, 3))

    def test_multiallelic_sites(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        self.assertEqual(vf.multiallelic_sites(), [3])

    def test_reset_samples(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf = vf.reset_samples(['Sarah', 'John'])
        self.assertEqual(vf.samples, ['Sarah', 'John'])

class TestBedFrame(unittest.TestCase):

    def test_intersect(self):
        bf1 = BedFrame.from_file(f'{fuc_dir()}/data/bed/1.bed')
        bf2 = BedFrame.from_file(f'{fuc_dir()}/data/bed/2.bed')
        bf3 = BedFrame.from_file(f'{fuc_dir()}/data/bed/3.bed')
        bf4 = bf1.intersect(bf2)
        self.assertEqual(bf3.meta, bf4.meta)
        self.assertEqual(bf3.data, bf4.data)

class TestFastqFrame(unittest.TestCase):

    def test_shape(self):
        qf = FastqFrame.from_file(f'{fuc_dir()}/data/fastq/1.fastq')
        self.assertEqual(qf.shape, 4)

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
