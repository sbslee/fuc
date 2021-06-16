from .. import api
from pathlib import Path
import sys

description = f"""
This command will check whether files/dirs exist. It will return 'True' if
they exist and 'False' otherwise. The command will look for stdin if there
are no arguments.

usage examples:
  $ cat files.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] check whether files/dirs exist',
        description=description,
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
