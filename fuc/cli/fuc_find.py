import os
from pathlib import Path

from .. import api

description = """
Retrieve absolute paths of files whose name matches a specified pattern,
optionally recursively.
"""

epilog = f"""
[Example] Retrieve VCF files in the current directory only:
  $ fuc {api.common._script_name()} "*.vcf"

[Example] Retrieve VCF files recursively:
  $ fuc {api.common._script_name()} "*.vcf" -r

[Example] Retrieve VCF files in a specific directory:
  $ fuc {api.common._script_name()} "*.vcf" -d /path/to/dir
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Retrieve absolute paths of files whose name matches a
specified pattern, optionally recursively."""
    )
    parser.add_argument(
        'pattern',
        help=
"""Filename pattern."""
    )
    parser.add_argument(
        '-r',
        '--recursive',
        action='store_true',
        help=
"""Turn on recursive retrieving."""
    )
    parser.add_argument(
        '-d',
        '--directory',
        metavar='PATH',
        default=os.getcwd(),
        help=
"""Directory to search in (default: current directory)."""
    )

def main(args):
    if args.recursive:
        for path in Path(args.directory).rglob(args.pattern):
            print(path.absolute())
    else:
        for path in Path(args.directory).glob(args.pattern):
            print(path.absolute())
