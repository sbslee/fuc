import sys

from .. import api

import pandas as pd

description = f"""
###################################################
# Rename the samples in a depth of coverage file. #
###################################################

There are three different renaming modes using the 'names' file:
  - 'MAP': Default mode. Requires two columns, old names in the first and new names in the second.
  - 'INDEX': Requires two columns, new names in the first and 0-based indicies in the second.
  - 'RANGE': Requires only one column of new names but '--range' must be specified.

Usage examples:
  $ fuc {api.common._script_name()} in.tsv old_new.tsv > out.tsv
  $ fuc {api.common._script_name()} in.tsv new_idx.tsv --mode INDEX > out.tsv
  $ fuc {api.common._script_name()} in.tsv new_only.tsv --mode RANGE --range 2 5 > out.tsv
  $ fuc {api.common._script_name()} in.tsv old_new.csv --sep , > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Rename the samples in a depth of coverage file.',
        description=description,
    )
    parser.add_argument(
        'tsv',
        help='TSV file (zipped or unzipped).'
    )
    parser.add_argument(
        'names',
        help='Text file containing information for renaming the samples.'
    )
    parser.add_argument(
        '--mode',
        metavar='TEXT',
        default='MAP',
        choices=['MAP', 'INDEX', 'RANGE'],
        help="Renaming mode (default: 'MAP') (choices: 'MAP', 'INDEX', 'RANGE')."
    )
    parser.add_argument(
        '--range',
        metavar='INT',
        type=int,
        nargs=2,
        help="Index range to use when renaming the samples. Applicable only with the 'RANGE' mode."
    )
    parser.add_argument(
        '--sep',
        metavar='TEXT',
        default='\t',
        help="Delimiter to use for reading the 'names' file (default: '\\t')."
    )

def main(args):
    cf = api.pycov.CovFrame.from_file(args.tsv)
    df = pd.read_table(args.names, sep=args.sep, header=None)
    if args.mode == 'MAP':
        names = dict(zip(df[0], df[1]))
        indicies = None
    elif args.mode == 'INDEX':
        names = df[0].to_list()
        indicies = df[1].to_list()
    else:
        names = df[0].to_list()
        indicies = tuple(args.range)
    cf = cf.rename(names, indicies=indicies)
    sys.stdout.write(cf.to_string())
