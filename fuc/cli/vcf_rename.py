import sys

from .. import api

import pandas as pd

description = f"""
This command will rename the samples in a VCF file.

There are three renaming modes: 'MAP', 'INDICIES', and 'RANGE'. The default
mode is 'MAP' in which case the 'names' file must contain two columns, one
for the old names and the other for the new names. If the mode is 'INDICIES'
the first column should be the new names and the second column must be
0-based indicies of the samples to be renamed. Lastly, in the 'RANGE' mode
only the first column is required but the 'range' argument must be specified.
For more details on the renaming modes, please visit the
'pyvcf.VcfFrame.rename' method's documentation page.

Usage examples:
  $ fuc {api.common._script_name()} in.vcf old_new.tsv > out.vcf
  $ fuc {api.common._script_name()} in.vcf new_idx.tsv --mode INDICIES > out.vcf
  $ fuc {api.common._script_name()} in.vcf new_only.tsv --mode RANGE --range 2 5 > out.vcf
  $ fuc {api.common._script_name()} in.vcf old_new.csv --sep , > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Rename the samples in a VCF file.',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')
    parser.add_argument('names', help='delimited text file')
    parser.add_argument('--mode', metavar='TEXT', default='MAP', choices=['MAP', 'INDICIES', 'RANGE'], help="renaming mode (default: 'MAP') (choices: 'MAP', 'INDICIES', 'RANGE')")
    parser.add_argument('--range', metavar='INT', type=int, nargs=2, help='specify an index range')
    parser.add_argument('--sep', metavar='TEXT', default='\t', help="delimiter to use (default: '\\t')")

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    df = pd.read_table(args.names, sep=args.sep, header=None)
    if args.mode == 'MAP':
        names = dict(zip(df[0], df[1]))
        indicies = None
    elif args.mode == 'INDICIES':
        names = df[0].to_list()
        indicies = df[1].to_list()
    else:
        names = df[0].to_list()
        indicies = tuple(args.range)
    vf = vf.rename(names, indicies=indicies)
    sys.stdout.write(vf.to_string())
