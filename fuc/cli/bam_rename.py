from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[BAM] add a new sample name to a BAM file',
        description='This command will add a new sample name to a BAM file.'
    )
    parser.add_argument('input_bam', help='input BAM file')
    parser.add_argument('name', help='sample name')
    parser.add_argument('output_bam', help='output BAM file')

def main(args):
    api.pybam.rename(args.input_bam, args.name, args.output_bam)
