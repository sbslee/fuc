from fuc.api.common import get_script_name
import pandas as pd
from pandas.api.types import is_numeric_dtype

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[TABLE] summarize a text file',
        description='This command will summarize a text file. It '
                    'essentially wraps the `pandas.DataFrame.describe` '
                    'method.'
    )
    parser.add_argument('text_file', help='text file')
    parser.add_argument('output_file', help='output file')
    parser.add_argument('--delimiter', metavar='TEXT', default='\t',
        help="delimiter (default: '\\t')")
    parser.add_argument('--skiprows', metavar='TEXT',
        help='comma-separated line numbers to skip (0-indexed) or '
             'number of lines to skip at the start of the file '
             "(e.g. '--skiprows 1,' will skip the second line, "
             "'--skiprows 2,4' will skip the third and fifth lines, and "
             "'--skiprows 10' will skip the first 10 lines)")
    parser.add_argument('--na_values', metavar='TEXT', nargs='+',
        help='additional strings to recognize as NA/NaN (by default, '
             "the following values are interpreted as NaN: '', "
             "'#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', "
             "'-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', "
             "'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null')")

def main(args):
    if ',' in args.skiprows:
        skiprows = [int(x) for x in args.skiprows.split(',') if x]
    else:
        skiprows = int(args.skiprows)
    df = pd.read_table(args.text_file, delimiter=args.delimiter,
        skiprows=skiprows, na_values=args.na_values)
    for header in df:
        column = df[header]
        if is_numeric_dtype(column):
            print(column.describe())
            print('='*70)
        else:
            print(column.value_counts())
            print('='*70)
