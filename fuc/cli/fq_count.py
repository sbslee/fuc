from .. import api
import sys

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[FASTQ] count sequence reads in FASTQ files',
        description=
            'This command will count sequence reads in FASTQ '
            'files (both zipped and unzipped). It will look for stdin '
            'if there are no arguments (e.g. $ cat files.list | fuc '
            f'{api.common.script_name(__file__)}).'
    )
    parser.add_argument('paths', nargs='*',
        help='FASTQ file paths (default: stdin)')
    return parser

def main(args):
    if args.paths:
        paths = args.paths
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        qf = api.pyfq.FqFrame.from_file(path)
        print(qf.df.shape[0])
