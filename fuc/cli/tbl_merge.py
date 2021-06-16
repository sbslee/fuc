from .. import api
import pandas as pd

description = f"""
This command will merge two table files using one or more shared columns.
It essentially wraps the 'pandas.DataFrame.merge' method from the pandas
package. For details on the merging algorithms, please visit the method's
documentation page.

usage examples:
  $ fuc {api.common._script_name()} left.tsv right.tsv > merged.tsv
  $ fuc {api.common._script_name()} left.csv right.tsv --lsep , > merged.tsv
  $ fuc {api.common._script_name()} left.tsv right.tsv --how outer > merged.tsv
"""

CHOICES = ['left', 'right', 'outer', 'inner', 'cross']

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[TABLE] merge two table files',
        description=description,
    )
    parser.add_argument('left', help='left file')
    parser.add_argument('right', help='right file')
    parser.add_argument('--how', metavar='TEXT', choices=CHOICES,
         default='inner', help=f'type of merge to be performed {CHOICES} '
        "(default: 'inner')")
    parser.add_argument('--on', metavar='TEXT', nargs='+',
        help='column names to join on')
    parser.add_argument('--lsep', metavar='TEXT', default='\t',
        help="delimiter to use for the left file (default: '\\t')")
    parser.add_argument('--rsep', metavar='TEXT', default='\t',
        help="delimiter to use for the right file (default: '\\t')")
    parser.add_argument('--osep', metavar='TEXT', default='\t',
        help="delimiter to use for the output file (default: '\\t')")
    return parser

def main(args):
    df1 = pd.read_table(args.left, sep=args.lsep)
    df2 = pd.read_table(args.right, sep=args.rsep)
    df3 = df1.merge(df2, on=args.on, how=args.how)
    print(df3.to_csv(sep=args.osep, index=False), end='')
