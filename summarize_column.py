import argparse

from api.DataFrame import DataFrame
from api.common import parse_where

def main():
    parser = argparse.ArgumentParser(description='This command will '
        'output a summary table for the target column in the input '
        'text file. The target column must be categorical. You can also '
        "use '--group_col' to group the observations by another "
        "categorical column. For filtering, you can use '--rows' to express "
        'SQLite WHERE clause which will select rows that meet certain '
        'criteria.')
    parser.add_argument('table_file', help='input table file')
    parser.add_argument('target_col', help='target column')
    parser.add_argument('--group_col', help='column to group by')
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
    i = df.get_index(args.target_col)
    if df.dtypes[i] != 'categorical':
        raise TypeError('target column must be categorical, not numeric')
    if not args.group_col:
        results = df.summarize_col(i)
        print(args.target_col, 'Count', sep='\t')
        for target, count in results.items():
            print(target, count, sep='\t')
    else:
        j = df.get_index(args.group_col)
        if df.dtypes[j] != 'categorical':
            raise TypeError('group column must be categorical, not numeric')
        # Prepare the results.
        targets = df.get_unique(i)
        groups = df.get_unique(j)
        results = {}
        for target in targets:
            results[target] = {}
            for group in groups:
                results[target][group] = 0
        # Get the results.
        for fields in df.get_data():
            target = fields[i]
            group = fields[j]
            results[target][group] += 1
        # Output the results.
        print('\t'.join([args.target_col] + list(groups)))
        for target, fields in results.items():
            print('\t'.join([target] + [str(v) for k, v in fields.items()]))

if __name__ == '__main__':
    main()
