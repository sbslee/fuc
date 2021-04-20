import pandas as pd
from fuc.api.common import get_script_name

CHOICES = ['left', 'right', 'outer', 'inner', 'cross']

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[TABLE] merge two text files',
        description='This command will merge two text files using one '
                    'or more shared columns. This essentially wraps the '
                    '`pandas.DataFrame.merge` method.'
    )
    parser.add_argument('left_file', help='left file')
    parser.add_argument('right_file', help='right file')
    parser.add_argument('output_file', help='output file')
    parser.add_argument('--how', metavar='TEXT', choices=CHOICES,
         default='inner', help=f'type of merge to be performed {CHOICES} '
        "(default: 'inner')")
    parser.add_argument('--on', metavar='TEXT', nargs='+',
        help='column names to join on')
    parser.add_argument('--left_delimiter', metavar='TEXT', default='\t',
        help="left delimiter (default: '\\t')")
    parser.add_argument('--right_delimiter', metavar='TEXT', default='\t',
        help="right delimiter (default: '\\t')")
    parser.add_argument('--output_delimiter', metavar='TEXT', default='\t',
        help="output delimiter (default: '\\t')")
    return parser

def main(args):
    df1 = pd.read_table(args.left_file, delimiter=args.left_delimiter)
    df2 = pd.read_table(args.right_file, delimiter=args.right_delimiter)
    df3 = df1.merge(df2, on=args.on, how=args.how)
    df3.to_csv(args.output_file, sep=args.output_delimiter, index=False)
