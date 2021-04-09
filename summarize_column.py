import argparse

from DataFrame import DataFrame
from common import parse_where

def main():
    parser = argparse.ArgumentParser(description='This command will '
        'output a summary table for the given categorical column in the '
        "input text file. For filtering, you can use '--rows' to express "
        'SQLite WHERE clause which will select rows that meet certain '
        'criteria.')
    parser.add_argument('table_file', help='input table file')
    parser.add_argument('groupby', help='input table file')
    parser.add_argument('column', help='input table file')
    parser.add_argument('--delimiter', default='\t',
        help="delimiter for the table (default: '\\t')")
    parser.add_argument('--skiprows', type=int, nargs='+',
        help='line numbers to skip')
    parser.add_argument('--rows', help='SQLite WHERE clause specifying '
        'which rows to summarize')
    parser.add_argument('--exclude_rows', action='store_true',
        help='use this tag to exclude specified rows')
    args = parser.parse_args()
    df = DataFrame.read(args.table_file, delimiter=args.delimiter,
        skiprows=args.skiprows)
    # Filter the rows.
    if args.rows:
        conditions = parse_where(args.rows)
        for k, v in conditions.items():
            df = df.filter_rows(k, v, args.exclude_rows)
    i = df.get_index(args.groupby)
    j = df.get_index(args.column)
    if df.dtypes[i] != 'categorical' or df.dtypes[j] != 'categorical':
        raise TypeError('must provide categorical columns')
    # Prepare the results.
    groups = df.get_unique(i)
    values = df.get_unique(j)
    results = {}
    for group in groups:
        results[group] = {}
        for value in values:
            results[group][value] = 0
    # Get the results.
    for fields in df.get_data():
        group = fields[i]
        value = fields[j]
        results[group][value] += 1
    # Output the results.
    print('\t'.join([args.groupby] + list(values)))
    for group, fields in results.items():
        print('\t'.join([group] + [str(v) for k, v in fields.items()]))

if __name__ == '__main__':
    main()
