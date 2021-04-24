import os
import unittest
import subprocess

from fuc.api.common import fuc_dir
from fuc import pyvcf, pybed, pyfq

if not os.path.exists('test_fuc'):
    os.mkdir('test_fuc')

vcf_file1 = f'{fuc_dir()}/data/vcf/1.vcf'
vcf_file2 = f'{fuc_dir()}/data/vcf/2.vcf'
vcf_file3 = f'{fuc_dir()}/data/vcf/3.vcf'

class TestPyvcf(unittest.TestCase):

    def test_shape(self):
        vf = pyvcf.read_file(vcf_file1)
        self.assertEqual(vf.shape, (5, 4))

    def test_filter_empty(self):
        vf = pyvcf.read_file(vcf_file2)
        vf = vf.filter_empty()
        vf.to_file('test_fuc/filter_empty.vcf')
        self.assertEqual(vf.df.shape, (4, 11))

    def test_filter_bed(self):
        vf = pyvcf.read_file(vcf_file1)
        bf = pybed.read_file(f'{fuc_dir()}/data/bed/1.bed')
        vf = vf.filter_bed(bf)
        self.assertEqual(vf.df.shape, (3, 13))

    def test_merge(self):
        vf1 = pyvcf.read_file(vcf_file1)
        vf2 = pyvcf.read_file(vcf_file2)
        vf3 = vf1.merge(vf2, how='outer', format='GT:DP')
        vf3.to_file('test_fuc/test_merge.vcf')
        self.assertEqual(vf3.df.shape, (9, 15))

    def test_compare(self):
        vf = pyvcf.read_file(vcf_file1)
        self.assertEqual(vf.compare('Steven', 'Sarah'), (0, 1, 1, 3))

    def test_filter_multiallelic(self):
        vf = pyvcf.read_file(vcf_file1)
        vf = vf.filter_multiallelic()

    def test_reset_samples(self):
        vf = pyvcf.read_file(vcf_file1)
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
        result = subprocess.run(['fuc', 'vfmerge', vcf_file1, vcf_file2, '--how', 'outer'], capture_output=True, text=True, check=True)
        self.assertEqual(len(result.stdout.strip().split('\n')), 10)

if __name__ == '__main__':
    unittest.main()
