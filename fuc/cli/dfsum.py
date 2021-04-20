from fuc.api.common import get_script_name
import pandas as pd

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[TABLE] summarize a text file',
        description='This command will summarize a text file. It '
                    'essentially wraps the `pandas.DataFrame.describe` '
                    'method.'
    )
    parser.add_argument('text_file', help='text file')
    parser.add_argument('--delimiter', metavar='TEXT', default='\t',
        help="delimiter (default: '\\t')")

def main(args):
    df = pd.read_table(args.text_file, delimiter=args.delimiter)
    print(df.describe(include='all'))
