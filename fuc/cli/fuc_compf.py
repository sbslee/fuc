import filecmp

from .. import api

description = """
Compare the contents of two files.

This command will compare the contents of two files, returning 'True' if they
are identical and 'False' otherwise.
"""

epilog = f"""
[Example] Compare two files:
  $ fuc {api.common._script_name()} left.txt right.txt
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compare the contents of two files."""
    )
    parser.add_argument(
        'left',
        help=
"""Input left file."""
    )
    parser.add_argument(
        'right',
        help=
"""Input right file."""
    )

def main(args):
    print(filecmp.cmp(args.left, args.right))
