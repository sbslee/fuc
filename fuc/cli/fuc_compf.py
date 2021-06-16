from .. import api
import filecmp

description = f"""
This command will compare the contents of two files. It will return 'True'
if they are identical and 'False' otherwise. It essentially wraps the
'filecmp.cmp' method.

usage examples:
  $ fuc {api.common._script_name()} 1.txt 2.txt
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] compare contents of two files',
        description=description,
    )
    parser.add_argument('file1', help='first file')
    parser.add_argument('file2', help='second file')

def main(args):
    print(filecmp.cmp(args.file1, args.file2))
