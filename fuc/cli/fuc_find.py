from .. import api
from pathlib import Path

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[FUC] find files with certain extension recursively',
        description=
            'This command will recursively find files with a certain '
            "extension -- such as '.txt' and '.vcf' -- within"
            'the given directory and return their absolute paths.'
    )
    parser.add_argument('path', help='directory path')
    parser.add_argument('extension', help='extension')

def main(args):
    for path in Path(args.path).rglob(f'*.{args.extension}'):
        print(path.absolute())
