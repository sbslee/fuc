import argparse

from api.DataFrame import DataFrame
from api.common import parse_where

def main():
    parser = argparse.ArgumentParser(description='This command will '
        'output a summary of the input text file. For each column, '
        'it will return the counts of unique records for categorical data '
        'and the summary statistics (minimum, maximum, mean, and median) '
        "for numeric data. You can use '--columns' to specify which "
        "columns should be displayed. For filtering, you can use '--rows' "
        'to express SQLite WHERE clause which will select rows that '
        'meet certain criteria.')
    parser.add_argument('table_file', help='input table file')
    parser.add_argument('--delimiter', default='\t',
        help="delimiter for the table (default: '\\t')")
    parser.add_argument('--columns', nargs='+',
        help='specify which columns to summarize')
    parser.add_argument('--exclude_columns', action='store_true',
        help='use this tag to exclude specified columns')
    parser.add_argument('--rows', help='SQLite WHERE clause specifying '
        'which rows to summarize')
    parser.add_argument('--exclude_rows', action='store_true',
        help='use this tag to exclude specified rows')
    args = parser.parse_args()
    df = DataFrame.read(args.table_file, delimiter=args.delimiter)
    # Filter the rows.
    if args.rows:
        conditions = parse_where(args.rows)
        for k, v in conditions.items():
            df = df.filter_rows(k, v, args.exclude_rows)
    # Filter the columns.
    if args.columns:
        df = df.filter_columns(args.columns, args.exclude_columns)
    # Summarize the columns.
    for i in range(df.shape[1]):
        print('#', df.head[i], f'({df.dtypes[i]})')
        results = df.summarize_col(i)
        keys = list(results)
        if len(results) > 10:
            print(keys[0], results[keys[0]])
            print(keys[1], results[keys[1]])
            print('...')
            print(keys[-2], results[keys[-2]])
            print(keys[-1], results[keys[-1]])
            print(f'[ found {len(results)} unique records ]')
        else:
            for k, v in results.items():
                print(k, v, sep='\t')
        print()

if __name__ == '__main__':
    main()
