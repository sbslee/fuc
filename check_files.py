import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description='This command checks whether or not files in the given list exist in the operating system.')
    parser.add_argument('input_file', help='input file containing the list of file paths')
    parser.add_argument('column_header', help='column header')
    parser.add_argument('--delimiter', default=',', 
        help="column delimiter (default: ',')")
    args = parser.parse_args()
    with open(args.input_file) as f:
        headers = next(f).strip().split(args.delimiter)
        i = headers.index(args.column_header)
        for line in f:
            fields = line.strip().split(args.delimiter)
            file_path = fields[i]
            print(Path(file_path).exists(), file_path)

if __name__ == '__main__':
    main()
