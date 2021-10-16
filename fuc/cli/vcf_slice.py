import sys

from .. import api

description = f"""
==============================================================================
Slice VCF file for specified regions.

Input VCF must be already BGZF compressed and indexed (.tbi) to allow random
access. Each region to be sliced must have the format chrom:start-end and be
a half-open interval with (start, end]. This means, for example, chr1:100-103
will extract positions 101, 102, and 103. Alternatively, you can provide a
BED file to specify regions.

Specify regions manually:
  $ fuc {api.common._script_name()} in.vcf.gz 1:100-300 2:400-700 > out.vcf

Speicfy regions with a BED file:
  $ fuc {api.common._script_name()} in.vcf.gz regions.bed > out.vcf
==============================================================================
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Slice VCF file for specified regions.',
        description=description,
    )
    parser.add_argument(
        'file',
        help='VCF file.'
    )
    parser.add_argument(
        'regions',
        help='One or more regions. Also accepts a BED file (zipped or \n'
             'unzipped).'
    )

def main(args):
    api.pyvcf.slice(args.file, args.regions, path='-')
