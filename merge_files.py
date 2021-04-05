import argparse
from DataFrame import DataFrame

def main():
    parser = argparse.ArgumentParser(description='This command '
        'will merge two text files on a shared column.')
    parser.add_argument('left_file', help='left file')
    parser.add_argument('right_file', help='right file')
    parser.add_argument('output_file', help='merged file')
    parser.add_argument('on', help='column name to join on', nargs='+')
    parser.add_argument('--left_delimiter', default='\t',
        help="delimiter for the left file (default: '\\t')")
    parser.add_argument('--right_delimiter', default='\t',
        help="delimiter for the right file (default: '\\t')")
    parser.add_argument('--output_delimiter', default='\t',
        help="delimiter for the output file (default: '\\t')")
    args = parser.parse_args()
    df1 = DataFrame.read(args.left_file, delimiter=args.left_delimiter)
    df2 = DataFrame.read(args.right_file, delimiter=args.right_delimiter)
    df3 = df1.merge(df2, args.on)
    df3.write(args.output_file, delimiter=args.output_delimiter)

if __name__ == '__main__':
    main()
