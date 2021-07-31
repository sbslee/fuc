from .. import api

import pandas as pd
from pandas.api.types import is_numeric_dtype

description = f"""
This command will summarize a table file. It essentially wraps the
'pandas.Series.describe' and 'pandas.Series.value_counts' methods from the
pandas pacakge.

Usage examples:
  $ fuc {api.common._script_name()} table.tsv
  $ fuc {api.common._script_name()} table.csv --sep ,
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Summarize a table file.',
        description=description,
    )
    parser.add_argument(
        'table_file',
        help='Table file.'
    )
    parser.add_argument(
        '--sep',
        metavar='TEXT',
        default='\t',
        help="Delimiter to use (default: '\\t')."
    )
    parser.add_argument(
        '--skiprows',
        metavar='TEXT',
        help=
            'Comma-separated line numbers to skip (0-indexed) or '
            'number of lines to skip at the start of the file '
            '(e.g. `--skiprows 1,` will skip the second line, '
            '`--skiprows 2,4` will skip the third and fifth lines, and '
            '`--skiprows 10` will skip the first 10 lines).'
    )
    parser.add_argument(
        '--na_values',
        metavar='TEXT',
        nargs='+',
        help=
            'Additional strings to recognize as NA/NaN (by default, '
            "the following values are interpreted as NaN: '', "
            "'#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', "
            "'-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', "
            "'NA', 'NULL', 'NaN', 'n/a', 'nan', 'null')."
    )
    parser.add_argument(
        '--keep_default_na',
        action='store_false',
        help=
            'Wwhether or not to include the default NaN values when '
            "parsing the data (see 'pandas.read_table' for details)."
    )
    parser.add_argument(
        '--expr',
        metavar='TEXT',
        help=
            'Query the columns of a pandas.DataFrame with a '
            """boolean expression (e.g. `--query "A == 'yes'"`)."""
    )
    parser.add_argument(
        '--columns',
        metavar='TEXT',
        nargs='+',
        help=
            'Columns to be summarized (by default, all columns '
            'will be included).'
    )
    parser.add_argument(
        '--dtypes',
        metavar='PATH',
        help=
            "File of column names and their data types (etheir "
            "'categorical' or 'numeric'); one tab-delimited pair of column "
            "name and data type per line."
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

    if args.dtypes is not None:
        with open(args.dtypes) as f:
            for line in f:
                fields = line.strip().split('\t')
                col = fields[0]
                dtype = fields[1]
                if dtype == 'numeric':
                    dtype = float
                elif dtype == 'categorical':
                    dtype = 'category'
                else:
                    raise TypeError("Data type must be either 'numeric' or "
                        "'categorical'.")
                df[col] = df[col].astype(dtype)

    if args.expr:
        df = df.query(args.expr)

    if args.columns:
        headers = args.columns
    else:
        headers = df.columns

    print('='*70)

    for header in headers:
        s = df[header]
        print(f'Name: {header}')
        print(f'Count: {len(s)}')
        print('Summary:')
        if is_numeric_dtype(s):
            print(s.describe().to_string())
        else:
            print(s.value_counts().sort_index().to_string())
        print('='*70)
