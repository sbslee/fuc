from .. import api
import pysam

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common._script_name(__file__),
        help='[BAM] index a SAM/BAM/CRAM file',
        description='This command will index a SAM/BAM/CRAM file.'
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    pysam.index(args.bam)
