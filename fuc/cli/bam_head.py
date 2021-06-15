from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common._script_name(),
        help='[BAM] print the header of a SAM/BAM/CRAM file',
        description='This command will print the header of a SAM/BAM/CRAM file.'
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    print(api.pybam.header(args.bam))
