import sys

from .. import api

import pysam
from fuc import common, pyvcf

description = """
Slice a GFF/BED/SAM/VCF file with Tabix.

After creating an index file (.tbi), the Tabix program is able to quickly
retrieve data lines overlapping regions specified in the format
'chr:start-end'. Coordinates specified in this region format are 1-based and
inclusive.
"""

epilog = f"""
[Example] Slice a VCF file:
  $ fuc {api.common._script_name()} in.vcf.gz chr1:100-200 > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Slice a GFF/BED/SAM/VCF file with Tabix."""
    )
    parser.add_argument(
        'file',
        help=
"""File to be sliced."""
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help=
"""One or more regions."""
    )

def main(args):
    if '.vcf' in args.file:
        vf = pyvcf.VcfFrame.from_file(args.file, meta_only=True)
        for line in vf.meta:
            sys.stdout.write(line + '\n')
        sys.stdout.write('#' + '\t'.join(vf.df.columns) + '\n')
    elif '.sam' in args.file:
        pass
    elif '.gff' in args.file:
        pass
    elif '.bed' in args.file:
        pass
    else:
        raise ValueError('Unsupported file format')

    tbx = pysam.TabixFile(args.file)

    for region in args.regions:
        for row in tbx.fetch(*common.parse_region(region)):
             sys.stdout.write(row + '\n')
