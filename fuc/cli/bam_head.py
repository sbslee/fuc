from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[BAM] print the header of a BAM file',
        description='This command will print the header of a BAM file.'
    )
    parser.add_argument('bam_file', help='BAM file')

def main(args):
    print(api.pybam.header(args.bam_file))
