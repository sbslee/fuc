import argparse
from DataFrame import DataFrame

def main():
    parser = argparse.ArgumentParser(description='This command will '
        'output a summary of the input text file.')
    parser.add_argument('table_file', help='input table file')
    parser.add_argument('--delimiter', default='\t',
        help="delimiter for the table (default: '\\t')")
    parser.add_argument('--columns', nargs='+',
        help='specify which columns to summarize')
    parser.add_argument('--exclude_columns', action='store_true',
        help='use this tag to exclude specified columns')
    args = parser.parse_args()
    df = DataFrame.read(args.table_file, delimiter=args.delimiter)

    if args.columns and args.exclude_columns:
        target_cols = [x for x in df.get_head() if x not in args.columns]
    elif args.columns and not args.exclude_columns:
        target_cols = [x for x in df.get_head() if x in args.columns]
    else:
        target_cols = df.get_head()

    for i in [df.get_index(x) for x in target_cols]:
        print('#', df.head[i], f'({df.dtypes[i]})')
        for k, v in df.summarize_col(i).items():
            print(k, v, sep='\t')

if __name__ == '__main__':
    main()
