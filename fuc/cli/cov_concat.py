import sys

from .. import api

import pysam

description = """
Concatenate depth of coverage files.
"""

epilog = f"""
[Example] Concatenate vertically:
  $ fuc {api.common._script_name()} in1.tsv in2.tsv > out.tsv

[Example] Concatenate horizontally:
  $ fuc {api.common._script_name()} in1.tsv in2.tsv --axis 1 > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Concatenate depth of coverage files."""
    )
    parser.add_argument(
        'tsv',
        nargs='+',
        help=
"""Input TSV files."""
    )
    parser.add_argument(
        '--axis',
        metavar='INT',
        default=0,
        type=int,
        help=
"""The axis to concatenate along (default: 0) (choices:
0, 1 where 0 is index and 1 is columns)."""
    )

def main(args):
    cfs = [api.pycov.CovFrame.from_file(x) for x in args.tsv]
    cf = api.pycov.concat(cfs, axis=args.axis)
    sys.stdout.write(cf.to_string())
