import sys

from .. import api

import pysam

description = f"""
############################################################
# Concatenate TSV files containing depth of coverage data. #
############################################################

Usage examples:
  $ fuc {api.common._script_name()} 1.tsv 2.tsv > rows.tsv
  $ fuc {api.common._script_name()} 1.tsv 2.tsv --axis 1 > cols.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Concatenate TSV files containing depth of coverage data.',
        description=description,
    )
    parser.add_argument(
        'tsv',
        metavar='PATH',
        nargs='+',
        help='One or more TSV files.'
    )
    parser.add_argument(
        '--axis',
        metavar='INT',
        default=0,
        type=int,
        help='The axis to concatenate along (default: 0) (chocies: 0, 1 where 0 is index and 1 is columns).'
    )

def main(args):
    cfs = [api.pycov.CovFrame.from_file(x) for x in args.tsv]
    cf = api.pycov.concat(cfs, axis=args.axis)
    sys.stdout.write(cf.to_string())
