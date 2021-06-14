import argparse

from .version import __version__
from .cli import commands

description = """
The main goal of the fuc package is to wrap some of the most frequently used
commands in the field of bioinformatics into one place. You can use fuc for
both command line interface (CLI) and application programming interface
(API) whose documentations are available at Read the Docs
(https://sbslee-fuc.readthedocs.io/en/latest).
"""

def main():
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'%(prog)s {__version__}',
        help='show the version number and exit'
    )
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='COMMAND',
        help='name of the command'
    )
    subparsers.required = True
    for name, command in commands.items():
        command.create_parser(subparsers)
    args = parser.parse_args()
    commands[args.command].main(args)

if __name__ == '__main__':
    main()
