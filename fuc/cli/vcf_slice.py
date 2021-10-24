import sys

from .. import api

description = f"""
Slice a VCF file for specified regions.
"""

epilog = f"""
[Example] Specify regions manually:
$ fuc {api.common._script_name()} in.vcf.gz 1:100-300 2:400-700 > out.vcf

[Example] Speicfy regions with a BED file:
$ fuc {api.common._script_name()} in.vcf.gz regions.bed > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Slice a VCF file for specified regions.',
    )
    parser.add_argument(
        'file',
        help='Input VCF file must be already BGZF compressed (.gz) and \n'
             'indexed (.tbi) to allow random access.'
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help='One or more regions to be sliced. Each region must have the \n'
             'format chrom:start-end and be a half-open interval with \n'
             '(start, end]. This means, for example, chr1:100-103 will \n'
             'extract positions 101, 102, and 103. Alternatively, you can \n'
             'provide a BED file (zipped or unzipped) to specify regions.'
    )

def main(args):
    api.pyvcf.slice(args.file, args.regions, path='-')
