from .. import api
import pysam

DESCRIPTION = f"""
This command will index a SAM/BAM/CRAM file. It essentially wraps the
'pysam.index' method.

usage examples:
  $ fuc {api.common.script()} in.bam
  $ fuc {api.common.script()} in.sam
  $ fuc {api.common.script()} in.cram
"""

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script(),
        help='[BAM] index a SAM/BAM/CRAM file',
        description=DESCRIPTION
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    pysam.index(args.bam)
