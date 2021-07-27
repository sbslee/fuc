import sys

from .. import api

import pysam

description = f"""
This command will compute read depth from the input SAM/BAM/CRAM files.

Usage examples:
  $ fuc {api.common._script_name()} 1.bam 2.bam --bed in.bed > out.tsv
  $ fuc {api.common._script_name()} in.bam --region chr1:100-200 > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] Compute read depth from SAM/BAM/CRAM files.',
        description=description,
    )
    parser.add_argument(
        'bam',
        nargs='+',
        help='One or more input SAM/BAM/CRAM files.'
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help='BED file.'
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Only report depth in specified region ('chrom:start-end')."
    )

def main(args):
    cf = api.pycov.CovFrame.from_bam(args.bam, bed=args.bed, region=args.region)
    sys.stdout.write(cf.to_string())
