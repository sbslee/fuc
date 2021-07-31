import sys

from .. import api

import pandas as pd

description = f"""
This command will filter a VCF file (both zipped and unzipped).

Usage examples:
  $ fuc {api.common._script_name()} in.vcf --expr 'GT == "0/0"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'GT != "0/0"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'DP < 30' > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'AD[1] < 10' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'AD[1] < 10 and DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'AD[1] < 10 or DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'AD[1] < 10 or DP < 30' --opposite > out.vcf
  $ fuc {api.common._script_name()} in.vcf --expr 'np.mean(AD) < 10' --greedy --samples sample.list > out.vcf
  $ fuc {api.common._script_name()} in.vcf --drop_duplicates CHROM POS REF ALT > out.vcf
  $ fuc {api.common._script_name()} in.vcf --filter_empty > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Filter a VCF file.',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')
    parser.add_argument(
        '--expr',
        metavar='TEXT',
        help='Expression to evaluate.'
    )
    parser.add_argument(
        '--samples',
        metavar='PATH',
        help='File of sample names to apply the marking (one sample per line).'
    )
    parser.add_argument(
        '--drop_duplicates',
        metavar='TEXT',
        nargs='*',
        help=
            'Only consider certain columns for identifying duplicates, '
            'by default use all of the columns.'
    )
    parser.add_argument(
        '--greedy',
        action='store_true',
        help='Use this flag to mark even ambiguous genotypes as missing.'
    )
    parser.add_argument(
        '--opposite',
        action='store_true',
        help=
            'Use this flag to mark all genotypes that do not satisfy the '
            'query expression as missing and leave those that do intact.'
    )
    parser.add_argument(
        '--filter_empty',
        action='store_true',
        help='Use this flag to remove rows with no genotype calls at all.'
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)

    if args.expr is not None:
        if args.samples is None:
            samples = None
        else:
            samples = pd.read_table(args.samples, header=None)[0].to_list()

        vf = vf.markmiss(args.expr,
                         greedy=args.greedy,
                         opposite=args.opposite,
                         samples=samples)

    if args.drop_duplicates is not None:
        vf = vf.drop_duplicates(args.drop_duplicates)

    if args.filter_empty:
        vf = vf.filter_empty()

    sys.stdout.write(vf.to_string())
