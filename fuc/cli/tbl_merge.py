import sys

from .. import api

import pandas as pd

description = f"""
This command will merge two table files using one or more shared columns.

The command essentially wraps the 'pandas.DataFrame.merge' method from the pandas package. For details on the merging algorithms, please visit the method's documentation page.

Usage examples:
  $ fuc {api.common._script_name()} left.tsv right.tsv > merged.tsv
  $ fuc {api.common._script_name()} left.csv right.tsv --lsep , > merged.tsv
  $ fuc {api.common._script_name()} left.tsv right.tsv --how outer > merged.tsv
"""

CHOICES = ['left', 'right', 'outer', 'inner', 'cross']

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Merge two table files.',
        description=description,
    )
    parser.add_argument(
        'left',
        help='Left file.'
    )
    parser.add_argument(
        'right',
        help='Right file.'
    )
    parser.add_argument(
        '--how',
        metavar='TEXT',
        choices=CHOICES,
        default='inner',
        help=f"Type of merge to be performed {CHOICES} (default: 'inner')."
    )
    parser.add_argument(
        '--on',
        metavar='TEXT',
        nargs='+',
        help='Column names to join on.'
    )
    parser.add_argument(
        '--lsep',
        metavar='TEXT',
        default='\t',
        help="Delimiter to use for the left file (default: '\\t')."
    )
    parser.add_argument(
        '--rsep',
        metavar='TEXT',
        default='\t',
        help="Delimiter to use for the right file (default: '\\t')."
    )
    parser.add_argument(
        '--osep',
        metavar='TEXT',
        default='\t',
        help="Delimiter to use for the output file (default: '\\t')."
    )
    return parser

def main(args):
    df1 = pd.read_table(args.left, sep=args.lsep)
    df2 = pd.read_table(args.right, sep=args.rsep)
    df3 = df1.merge(df2, on=args.on, how=args.how)
    sys.stdout.write(df3.to_csv(sep=args.osep, index=False))
