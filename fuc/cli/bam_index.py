from .. import api

import pysam

description = f"""
This command will index the input SAM/BAM/CRAM file.

Usage examples:
  $ fuc {api.common._script_name()} in.sam
  $ fuc {api.common._script_name()} in.bam
  $ fuc {api.common._script_name()} in.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Index a SAM/BAM/CRAM file.',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file.')

def main(args):
    pysam.index(args.bam)
