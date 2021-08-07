import sys
from pathlib import Path

from .. import api

description = f"""
This command will check whether or not specified files including directoires exist, returning 'True' if they exist and 'False' otherwise.

The command will look for stdin if there are no arguments.

Usage examples:
  $ fuc {api.common._script_name()} test.txt
  $ fuc {api.common._script_name()} test_dir
  $ cat test.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Check whether certain files exist.',
        description=description,
    )
    parser.add_argument('files', nargs='*',
        help='Files and directories to be tested (default: stdin).')

def main(args):
    if args.files:
        paths = args.files
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        print(Path(path).exists(), path)
