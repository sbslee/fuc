from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common._script_name(__file__),
        help='[BAM] rename the sample in a BAM file',
        description='This command will rename the sample in a BAM file.'
    )
    parser.add_argument('input_bam', help='input BAM file')
    parser.add_argument('name', help='sample name')
    parser.add_argument('output_bam', help='output BAM file')

def main(args):
    api.pybam.rename(args.input_bam, args.name, args.output_bam)
