from .. import api
import pandas as pd

CHOICES = ['left', 'right', 'outer', 'inner', 'cross']

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script(),
        help='[TABLE] merge two table files',
        description=
            'This command will merge two table files using one '
            'or more shared columns. This essentially wraps the '
            '`pandas.DataFrame.merge` method.'
    )
    parser.add_argument('left_file', help='left table file')
    parser.add_argument('right_file', help='right table file')
    parser.add_argument('--how', metavar='TEXT', choices=CHOICES,
         default='inner', help=f'type of merge to be performed {CHOICES} '
        "(default: 'inner')")
    parser.add_argument('--on', metavar='TEXT', nargs='+',
        help='column names to join on')
    parser.add_argument('--left_sep', metavar='TEXT', default='\t',
        help="left delimiter to use (default: '\\t')")
    parser.add_argument('--right_sep', metavar='TEXT', default='\t',
        help="right delimiter to use (default: '\\t')")
    parser.add_argument('--output_sep', metavar='TEXT', default='\t',
        help="output delimiter to use (default: '\\t')")
    return parser

def main(args):
    df1 = pd.read_table(args.left_file, sep=args.left_sep)
    df2 = pd.read_table(args.right_file, sep=args.right_sep)
    df3 = df1.merge(df2, on=args.on, how=args.how)
    print(df3.to_csv(sep=args.output_sep, index=False))
