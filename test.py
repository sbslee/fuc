import os
import unittest
import subprocess

from fuc.api.common import fuc_dir
from fuc import pyvcf, pybed, pyfq

fuc_dir = fuc_dir()
vcf_file1 = f'{fuc_dir}/data/vcf/1.vcf'
vcf_file2 = f'{fuc_dir}/data/vcf/2.vcf'
vcf_file3 = f'{fuc_dir}/data/vcf/3.vcf'
bed_file1 = f'{fuc_dir}/data/bed/1.bed'
bed_file2 = f'{fuc_dir}/data/bed/2.bed'
fq_file1 = f'{fuc_dir}/data/fq/1.fastq'
text_file1 = f'{fuc_dir}/data/text/1.txt'
text_file2 = f'{fuc_dir}/data/text/2.txt'

class TestPyvcf(unittest.TestCase):

    def test_shape(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file1)
        self.assertEqual(vf.shape, (5, 4))

    def test_filter_empty(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file2)
        vf = vf.filter_empty()
        self.assertEqual(vf.df.shape, (4, 11))

    def test_filter_bed(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file1)
        bf = pybed.read_file(f'{fuc_dir}/data/bed/1.bed')
        vf = vf.filter_bed(bf)
        self.assertEqual(vf.df.shape, (3, 13))

    def test_merge(self):
        vf1 = pyvcf.VcfFrame.from_file(vcf_file1)
        vf2 = pyvcf.VcfFrame.from_file(vcf_file2)
        vf3 = vf1.merge(vf2, how='outer', format='GT:DP')
        self.assertEqual(vf3.df.shape, (9, 15))

    def test_compare(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file1)
        self.assertEqual(vf.compare('Steven', 'Sarah'), (1, 1, 0, 3))

    def test_filter_multialt(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file1)
        vf = vf.filter_multialt()
        self.assertEqual(vf.shape[0], 4)

    def test_subset(self):
        vf = pyvcf.VcfFrame.from_file(vcf_file1)
        vf = vf.subset(['Sarah', 'John'])
        self.assertEqual(len(vf.samples), 2)

class TestPybed(unittest.TestCase):

    def test_intersect(self):
        bf1 = pybed.read_file(f'{fuc_dir}/data/bed/1.bed')
        bf2 = pybed.read_file(f'{fuc_dir}/data/bed/2.bed')
        bf3 = pybed.read_file(f'{fuc_dir}/data/bed/3.bed')
        bf4 = bf1.intersect(bf2)
        self.assertEqual(bf3.to_string(), bf4.to_string())

class TestPyfq(unittest.TestCase):

    def test_shape(self):
        qf = pyfq.read_file(f'{fuc_dir}/data/fq/1.fastq')
        self.assertEqual(qf.shape, (5, 4))

class TestCli(unittest.TestCase):
    def test_bfintxn(self):
        result = subprocess.run(['fuc', 'bed_intxn', bed_file1, bed_file2], capture_output=True, text=True, check=True)
        self.assertEqual(len(result.stdout.split('\n')), 6)

    def test_bfsum(self):
        result = subprocess.run(['fuc', 'bed_sum', bed_file1], capture_output=True, text=True, check=True)
        self.assertTrue('Total' in result.stdout)

    def test_dfmerge(self):
        result = subprocess.run(['fuc', 'tbl_merge', text_file1, text_file2], capture_output=True, text=True, check=True)
        self.assertTrue('Sarah' in result.stdout)

    def test_dfsum(self):
        result = subprocess.run(['fuc', 'tbl_sum', text_file1], capture_output=True, text=True, check=True)
        self.assertTrue('max' in result.stdout)

    def test_fuccompf(self):
        result = subprocess.run(['fuc', 'fuc_compf', vcf_file1, vcf_file1], capture_output=True, text=True, check=True)
        self.assertTrue('True' in result.stdout)

    def test_fucexist(self):
        result = subprocess.run(['fuc', 'fuc_exist', vcf_file1], capture_output=True, text=True, check=True)
        self.assertTrue('True' in result.stdout)

    def test_fucfind(self):
        result = subprocess.run(['fuc', 'fuc_find', fuc_dir, 'vcf'], capture_output=True, text=True, check=True)
        self.assertTrue('1.vcf' in result.stdout)

    def test_qfcount(self):
        result = subprocess.run(['fuc', 'fq_count', fq_file1], capture_output=True, text=True, check=True)
        self.assertEqual(int(result.stdout.strip()), 5)

    def test_qfsum(self):
        result = subprocess.run(['fuc', 'fq_sum', fq_file1], capture_output=True, text=True, check=True)
        self.assertTrue('# Total: 5' in result.stdout)

    def test_vfmerge(self):
        result = subprocess.run(['fuc', 'vcf_merge', vcf_file1, vcf_file2, '--how', 'outer'], capture_output=True, text=True, check=True)
        self.assertEqual(len(result.stdout.strip().split('\n')), 10)

    def test_vfslice(self):
        result = subprocess.run(['fuc', 'vcf_slice', vcf_file1, 'chr1', '--end', '300'], capture_output=True, text=True, check=True)
        self.assertEqual(len(result.stdout.strip().split('\n')), 5)

if __name__ == '__main__':
    unittest.main()
