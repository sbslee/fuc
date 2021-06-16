from .. import api
import pysam

description = f"""
This command will index a SAM/BAM/CRAM file. It essentially wraps the
'pysam.index' method from the pysam package.

usage examples:
  $ fuc {api.common._script_name()} in.sam
  $ fuc {api.common._script_name()} in.bam
  $ fuc {api.common._script_name()} in.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] index a SAM/BAM/CRAM file',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    pysam.index(args.bam)
