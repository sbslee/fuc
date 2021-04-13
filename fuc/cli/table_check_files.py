from fuc.api.common import get_script_name
from pathlib import Path

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='check whether files exist',
        description='This command will checks whether or not files in the '
            'given list exist inthe operating system.'
    )
    parser.add_argument('input_file', help='input file containing '
        'the list of file paths')
    parser.add_argument('column_header', help='column header')
    parser.add_argument('--delimiter', default=',',
        help="column delimiter (default: ',')")

def main(args):
    with open(args.input_file) as f:
        headers = next(f).strip().split(args.delimiter)
        i = headers.index(args.column_header)
        for line in f:
            fields = line.strip().split(args.delimiter)
            file_path = fields[i]
            print(Path(file_path).exists(), file_path)
