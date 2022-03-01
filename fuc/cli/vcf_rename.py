import sys

from .. import api

import pandas as pd

description = """
Rename the samples in a VCF file.

There are three different renaming modes using the 'names' file:
  - 'MAP': Default mode. Requires two columns, old names in the first
    and new names in the second.
  - 'INDEX': Requires two columns, new names in the first and 0-based
    indicies in the second.
  - 'RANGE': Requires only one column of new names but '--range' must
    be specified.
"""

epilog = f"""
[Example] Using the default 'MAP' mode:
  $ fuc {api.common._script_name()} in.vcf old_new.tsv > out.vcf

[Example] Using the 'INDEX' mode:
  $ fuc {api.common._script_name()} in.vcf new_idx.tsv --mode INDEX > out.vcf

[Example] Using the 'RANGE' mode:
  $ fuc {api.common._script_name()} in.vcf new_only.tsv --mode RANGE --range 2 5 > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Rename the samples in a VCF file."""
    )
    parser.add_argument(
        'vcf',
        help=
"""VCF file (compressed or uncompressed)."""
    )
    parser.add_argument(
        'names',
        help=
"""Text file containing information for renaming the samples."""
    )
    parser.add_argument(
        '--mode',
        metavar='TEXT',
        default='MAP',
        choices=['MAP', 'INDEX', 'RANGE'],
        help=
"""Renaming mode (default: 'MAP') (choices: 'MAP',
'INDEX', 'RANGE')."""
    )
    parser.add_argument(
        '--range',
        metavar='INT',
        type=int,
        nargs=2,
        help=
"""Index range to use when renaming the samples.
Applicable only with the 'RANGE' mode."""
    )
    parser.add_argument(
        '--sep',
        metavar='TEXT',
        default='\t',
        help=
"""Delimiter to use for reading the 'names' file 
(default: '\\t')."""
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
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
    vf = vf.rename(names, indicies=indicies)
    sys.stdout.write(vf.to_string())
