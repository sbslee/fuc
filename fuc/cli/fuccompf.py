from fuc.api.common import get_script_name
import filecmp

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[FUC] compare two files',
        description='This command will compare two files.'
    )
    parser.add_argument('file1', help='first file')
    parser.add_argument('file2', help='second file')

def main(args):
    print(filecmp.cmp(args.file1, args.file2))
