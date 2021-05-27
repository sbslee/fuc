from .. import api
from pathlib import Path
import sys

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[FUC] check whether files/dirs exist',
        description='This command will check whether files/dirs exist. '
                    "It will return 'True' if they exist and 'False' "
                    'otherwise. The command will look for stdin if there '
                    'are no arguments (e.g. $ cat files.list | fuc '
                    f'{api.common.script_name(__file__)}).'
    )
    parser.add_argument('paths', nargs='*',
        help='file/dir paths (default: stdin)')

def main(args):
    if args.paths:
        paths = args.paths
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        print(Path(path).exists(), path)
