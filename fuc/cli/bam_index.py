from .. import api
import pysam

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[BAM] index a BAM file',
        description='This command will index a BAM file.'
    )
    parser.add_argument('bam_file', help='BAM file')

def main(args):
    pysam.index(args.bam_file)
