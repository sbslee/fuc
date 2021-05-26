from .. import api
import pysam

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[BAM] slice a BAM file',
        description=
            'This command will slice a BAM file. By default, the command '
            'will create a accompanying index file (.bai) for the output '
            'BAM file.'
    )
    parser.add_argument(
        'input_bam',
        help='input BAM file'
    )
    parser.add_argument(
        'region',
        help='target region'
    )
    parser.add_argument(
        'output_bam',
        help='output BAM file'
    )
    parser.add_argument(
        '--no_index',
        action='store_true',
        help='use to this flag to skip indexing'
    )

def main(args):
    b = pysam.view('-b', args.input_bam, args.region)
    with open(args.output_bam, 'wb') as f:
        f.write(b)
    if not args.no_index:
        pysam.index(args.output_bam)
