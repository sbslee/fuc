from .. import api
import filecmp

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[FUC] compare contents of two files',
        description=
            'This command will compare the contents of two files. '
            "It will return 'True' if they are identical and 'False' "
            'otherwise.'
    )
    parser.add_argument('file1', help='first file')
    parser.add_argument('file2', help='second file')

def main(args):
    print(filecmp.cmp(args.file1, args.file2))
