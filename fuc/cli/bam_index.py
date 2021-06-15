from .. import api
import pysam

COMMAND = api.common._script_name()

DESCRIPTION = f"""
This command will index a SAM/BAM/CRAM file. It essentially wraps the
'pysam.index' method.

usage examples:
  $ fuc {COMMAND} in.bam
  $ fuc {COMMAND} in.sam
  $ fuc {COMMAND} in.cram
"""

def create_parser(subparsers):
    parser = subparsers.add_parser(
        COMMAND,
        help='[BAM] index a SAM/BAM/CRAM file',
        description=DESCRIPTION
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    pysam.index(args.bam)
