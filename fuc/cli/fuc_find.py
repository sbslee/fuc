import os
from pathlib import Path

from .. import api

description = f"""
This command will recursively find all files with a certain extension and then return their absolute paths.

Usage examples:
  $ fuc {api.common._script_name()} .vcf
  $ fuc {api.common._script_name()} .vcf.gz
  $ fuc {api.common._script_name()} .vcf.gz --dir ~/test_dir
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] Find all files with a certain extension recursively.',
        description=description,
    )
    parser.add_argument(
        'ext',
        help='File extension.'
    )
    parser.add_argument(
        '--dir',
        metavar='PATH',
        default=os.getcwd(),
        help='Directory to search in (default: current directory).'
    )

def main(args):
    for path in Path(args.dir).rglob(f'*{args.ext}'):
        print(path.absolute())
