import sys

from .. import api

import pysam

description = """
Index a VCF file.

This command will create an index file (.tbi) for the input VCF.
"""

epilog = f"""
[Example] Index a compressed VCF file:
  $ fuc {api.common._script_name()} in.vcf.gz

[Example] Index an uncompressed VCF file (will create a compressed VCF first):
  $ fuc {api.common._script_name()} in.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Index a VCF file.',
    )
    parser.add_argument(
        'vcf',
        help=
"""Input VCF file to be indexed. When an uncompressed file is
given, the command will automatically create a BGZF
compressed copy of the file (.gz) before indexing."""
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help=
"""Force to overwrite the index file if it is already present."""
    )

def main(args):
    pysam.tabix_index(args.vcf, preset='vcf', force=args.force)
