import filecmp

from .. import api

description = f"""
This command will compare the contents of two files, returning 'True' if they are identical and 'False' otherwise.

Usage examples:
  $ fuc {api.common._script_name()} left.txt right.txt
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Compare the contents of two files.',
        description=description,
    )
    parser.add_argument('left', help='Left file.')
    parser.add_argument('right', help='Right file.')

def main(args):
    print(filecmp.cmp(args.left, args.right))
