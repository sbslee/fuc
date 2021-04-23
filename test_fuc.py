import unittest
import subprocess

from fuc.api.common import fuc_dir
from fuc import pyvcf, pybed, pyfq

class TestPyvcf(unittest.TestCase):

    def test_filter_empty(self):
        vf = pyvcf.read_file(f'{fuc_dir()}/data/vcf/2.vcf')
        vf = vf.filter_empty()
        self.assertEqual(vf.df.shape, (4, 11))

    def test_filter_bed(self):
        vf = pyvcf.read_file(f'{fuc_dir()}/data/vcf/1.vcf')
        bf = pybed.read_file(f'{fuc_dir()}/data/bed/1.bed')
        vf = vf.filter_bed(bf)
        self.assertEqual(vf.df.shape, (3, 13))

    def test_merge(self):
        vf1 = pyvcf.read_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf2 = pyvcf.read_file(f'{fuc_dir()}/data/vcf/2.vcf')
        vf3 = vf1.merge(vf2, how='outer', format='GT:DP')
        self.assertEqual(vf3.df.shape, (8, 15))

    def test_compare(self):
        vf = pyvcf.read_file(f'{fuc_dir()}/data/vcf/1.vcf')
        self.assertEqual(vf.compare('Steven', 'Sarah'), (0, 1, 1, 3))

    def test_filter_multiallelic(self):
        vf = pyvcf.read_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf = vf.filter_multiallelic()

    def test_reset_samples(self):
        vf = pyvcf.read_file(f'{fuc_dir()}/data/vcf/1.vcf')
        vf = vf.reset_samples(['Sarah', 'John'])

class TestPybed(unittest.TestCase):

    def test_intersect(self):
        bf1 = pybed.read_file(f'{fuc_dir()}/data/bed/1.bed')
        bf2 = pybed.read_file(f'{fuc_dir()}/data/bed/2.bed')
        bf3 = pybed.read_file(f'{fuc_dir()}/data/bed/3.bed')
        bf4 = bf1.intersect(bf2)
        self.assertEqual(bf3.to_string(), bf4.to_string())

class TestPyfq(unittest.TestCase):

    def test_shape(self):
        qf = pyfq.read_file(f'{fuc_dir()}/data/fq/1.fastq')
        self.assertEqual(qf.df.shape, (5, 4))

class TestCli(unittest.TestCase):

    def test_qfcount(self):
        result = subprocess.run(['fuc', 'qfcount', f'{fuc_dir()}/data/fq/1.fastq'], capture_output=True, text=True, check=True)
        self.assertEqual(int(result.stdout.strip()), 5)

    def test_vfmerge(self):
        result = subprocess.run(['fuc', 'vfmerge', f'{fuc_dir()}/data/vcf/1.vcf', f'{fuc_dir()}/data/vcf/2.vcf', '--how', 'outer'], capture_output=True, text=True, check=True)
        self.assertEqual(len(result.stdout.strip().split('\n')), 9)

if __name__ == '__main__':
    unittest.main()
