import sys

from .. import api

import pysam
from fuc import common, pyvcf

description = f"""
#################################
# Slice a GFF/BED/SAM/VCF file. #
#################################

The Tabix program is used to index a TAB-delimited genome position file (GFF/BED/SAM/VCF) and create an index file (.tbi). The input data file must be position sorted and compressed by bgzip.

Usage examples:
  $ fuc {api.common._script_name()} in.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Slice a GFF/BED/SAM/VCF file.',
        description=description,
    )
    parser.add_argument(
        'file',
        help='File to be sliced.'
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help='One or more regions.'
    )

def main(args):
    vf = pyvcf.VcfFrame.from_file(args.file, meta_only=True)

    for line in vf.meta:
        sys.stdout.write(line + '\n')

    sys.stdout.write('\t'.join(vf.df.columns) + '\n')

    tbx = pysam.TabixFile(args.file)

    for region in args.regions:
        for row in tbx.fetch(*common.parse_region(region)):
             sys.stdout.write(row + '\n')
