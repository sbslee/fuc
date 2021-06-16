from .. import api
import pandas as pd
from pandas.api.types import is_numeric_dtype

description = f"""
This command will summarize a table file. It essentially wraps the
'pandas.Series.describe' and 'pandas.Series.value_counts' methods from the
pandas pacakge.

usage examples:
  $ fuc {api.common._script_name()} table.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[TABLE] summarize a table file',
        description=description,
    )
    parser.add_argument(
        'table_file',
        help='table file'
    )
    parser.add_argument(
        '--sep',
        metavar='TEXT',
        default='\t',
        help="delimiter to use (default: '\\t')"
    )
    parser.add_argument(
        '--skiprows',
        metavar='TEXT',
        help='comma-separated line numbers to skip (0-indexed) or '
             'number of lines to skip at the start of the file '
             '(e.g. `--skiprows 1,` will skip the second line, '
             '`--skiprows 2,4` will skip the third and fifth lines, and '
             '`--skiprows 10` will skip the first 10 lines)'
    )
    parser.add_argument(
        '--na_values',
        metavar='TEXT',
        nargs='+',
        help='additional strings to recognize as NA/NaN (by default, '
             "the following values are interpreted as NaN: '', "
             "'#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', "
             "'-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', "
             "'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null')"
    )
    parser.add_argument(
        '--keep_default_na',
        action='store_false',
        help='whether or not to include the default NaN values when '
             'parsing the data (see `pandas.read_table` for details)'
    )
    parser.add_argument(
        '--query',
        metavar='TEXT',
        help='query the columns of a pandas.DataFrame with a '
             """boolean expression (e.g. `--query "A == 'yes'"`)"""
    )
    parser.add_argument(
        '--columns',
        metavar='TEXT',
        nargs='+',
        help='columns to be summarized (by default, all columns '
             'will be included)'
    )

def main(args):
    if args.skiprows is None:
        skiprows = None
    elif ',' in args.skiprows:
        skiprows = [int(x) for x in args.skiprows.split(',') if x]
    else:
        skiprows = int(args.skiprows)
    df = pd.read_table(args.table_file,
                       sep=args.sep,
                       skiprows=skiprows,
                       na_values=args.na_values,
                       keep_default_na=args.keep_default_na)
    if args.query:
        df = df.query(args.query)
    if args.columns:
        headers = args.columns
    else:
        headers = df.columns
    print('='*70)
    for header in headers:
        column = df[header]
        if is_numeric_dtype(column):
            print(column.describe())
            print('='*70)
        else:
            print(column.value_counts())
            print('='*70)
