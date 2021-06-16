from .. import api
import sys

description = f"""
This command will count sequence reads in FASTQ files (both zipped and
unzipped). It will look for stdin if there are no arguments.

usage examples:
  $ fuc {api.common._script_name()} in.fastq
  $ cat fastq.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FASTQ] count sequence reads in FASTQ files',
        description=description,
    )
    parser.add_argument('fastq', nargs='*',
        help='FASTQ files (default: stdin)')
    return parser

def main(args):
    if args.fastq:
        paths = args.fastq
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        qf = api.pyfq.FqFrame.from_file(path)
        print(qf.df.shape[0])
