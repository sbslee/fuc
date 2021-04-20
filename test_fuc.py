import unittest

from fuc.api.common import fuc_dir
from fuc.api.VcfFrame import VcfFrame
from fuc.api.BedFrame import BedFrame
from fuc.api.FastqFrame import FastqFrame

class TestVcfFrame(unittest.TestCase):

    def test_filter_bed(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        bf = BedFrame.from_file(f'{fuc_dir()}/data/bed/1.bed')
        vf = vf.filter_bed(bf)
        self.assertEqual(vf.df.shape, (3, 13))

    def test_merge(self):
        vf1 = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf2 = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/2.vcf')
        vf3 = vf1.merge(vf2, how='outer', format='GT:DP')
        self.assertEqual(vf3.df.shape, (7, 14))

    def test_compare(self):
        vf = VcfFrame.from_file(f'{fuc_dir()}/data/vcf/1.vcf')
        self.assertEqual(vf.compare('Steven', 'Sarah'), (0, 1, 1, 3))

class TestBedFrame(unittest.TestCase):

    def test_intersect(self):
        bf1 = BedFrame.from_file(f'{fuc_dir()}/data/bed/1.bed')
        bf2 = BedFrame.from_file(f'{fuc_dir()}/data/bed/2.bed')
        bf3 = BedFrame.from_file(f'{fuc_dir()}/data/bed/3.bed')
        bf4 = bf1.intersect(bf2)
        self.assertEqual(bf3.to_string(), bf4.to_string())

class TestFastqFrame(unittest.TestCase):

    def test_shape(self):
        qf = FastqFrame.from_file(f'{fuc_dir()}/data/fastq/1.fastq')
        self.assertEqual(qf.shape, 4)

if __name__ == '__main__':
    unittest.main()
