from .. import api
from pathlib import Path
import sys

description = f"""
This command will check whether files/directories exist. It will return
'True' if they exist and 'False' otherwise. The command will look for stdin
if there are no arguments.

usage examples:
  $ fuc {api.common._script_name()} test.txt
  $ fuc {api.common._script_name()} test_dir
  $ cat test.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] check whether files/directories exist',
        description=description,
    )
    parser.add_argument('files', nargs='*',
        help='test files/directories (default: stdin)')

def main(args):
    if args.files:
        paths = args.files
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        print(Path(path).exists(), path)
