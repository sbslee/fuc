import sys
from pathlib import Path

from .. import api

description = """
Check whether certain files exist.

This command will check whether or not specified files including directories
exist, returning 'True' if they exist and 'False' otherwise.
"""

epilog = f"""
[Example] Test a file:
  $ fuc {api.common._script_name()} in.txt

[Example] Test a directory:
  $ fuc {api.common._script_name()} dir

[Example] When the input is stdin:
  $ cat test.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Check whether certain files exist."""
    )
    parser.add_argument(
        'files',
        nargs='*',
        help=
"""Files and directories to be tested (default: stdin)."""
    )

def main(args):
    if args.files:
        paths = args.files
    elif not sys.stdin.isatty():
        paths = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('no input files detected')
    for path in paths:
        print(Path(path).exists(), path)
