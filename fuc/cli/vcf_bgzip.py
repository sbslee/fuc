import sys

from .. import api
from Bio import bgzf

description = f"""
#####################################
# Compress a VCF file using bgzip.  #
#####################################

Usage examples:
  $ fuc {api.common._script_name()} in.vcf out.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Compress a VCF file using bgzip.',
        description=description,
    )
    parser.add_argument(
        'input',
        help='VCF file.'
    )
    parser.add_argument(
        'output',
        help='Compressed VCF file.'
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.input)
    vf.to_file(args.output, compression=True)
