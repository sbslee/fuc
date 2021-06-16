from .. import api
import filecmp

description = f"""
This command will compare the contents of two files. It will return 'True'
if they are identical and 'False' otherwise. It essentially wraps the
'filecmp.cmp' method from Python.

usage examples:
  $ fuc {api.common._script_name()} left.txt right.txt
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[FUC] compare contents of two files',
        description=description,
    )
    parser.add_argument('left', help='left file')
    parser.add_argument('right', help='right file')

def main(args):
    print(filecmp.cmp(args.left, args.right))
