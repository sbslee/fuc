from .. import api
from pathlib import Path

description = f"""
This command will recursively find files with a certain extension, such as
'.txt' and '.vcf', within the given directory and return their absolute paths.

usage examples:
  $ fuc {api.common._script_name()} path extension
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] find files with certain extension recursively',
        description=description,
    )
    parser.add_argument('path', help='directory path')
    parser.add_argument('extension', help='extension')

def main(args):
    for path in Path(args.path).rglob(f'*.{args.extension}'):
        print(path.absolute())
