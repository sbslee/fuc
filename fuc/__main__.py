import argparse

from .version import __version__
from .cli import commands

def main():
    parser = argparse.ArgumentParser()
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
        help='name of the command',
        required=True,
    )
    for name, command in commands.items():
        command.create_parser(subparsers)
    args = parser.parse_args()
    commands[args.command].main(args)

if __name__ == '__main__':
    main()
