import sys

from .. import api

import pysam

description = """
Index a GFF/BED/SAM/VCF file with Tabix.

The Tabix program is used to index a TAB-delimited genome position file
(GFF/BED/SAM/VCF) and create an index file (.tbi). The input data file must
be position sorted and compressed by bgzip.
"""

epilog = f"""
[Example] Index a GFF file:
  $ fuc {api.common._script_name()} in.gff.gz

[Example] Index a BED file:
  $ fuc {api.common._script_name()} in.bed.gz

[Example] Index a SAM file:
  $ fuc {api.common._script_name()} in.sam.gz

[Example] Index a VCF file:
  $ fuc {api.common._script_name()} in.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Index a GFF/BED/SAM/VCF file with Tabix."""
    )
    parser.add_argument(
        'file',
        help=
"""File to be indexed."""
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help=
"""Force to overwrite the index file if it is present."""
    )

def main(args):
    if '.vcf' in args.file:
        preset = 'vcf'
    elif '.sam' in args.file:
        preset = 'sam'
    elif '.gff' in args.file:
        preset = 'gff'
    elif '.bed' in args.file:
        preset = 'bed'
    else:
        raise ValueError('Unsupported file format')

    pysam.tabix_index(args.file, preset=preset, force=args.force)
