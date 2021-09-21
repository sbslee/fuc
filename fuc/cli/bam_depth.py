import sys

from .. import api

import pysam

description = f"""
###############################################
# Compute read depth from SAM/BAM/CRAM files. #
###############################################

Alignment files must be specified with either '--bam' or '--fn', but it's an error to use both.

By default, the command will count all reads within the alignment files. You can specify target regions with either '--bed' or '--region', but not both. When you do this, pay close attention to the 'chr' string in contig names (e.g. 'chr1' vs. '1'). Note also that '--region' requires the input files be indexed.

Under the hood, the command computes read depth using the 'samtools depth' command.

Usage examples:
  $ fuc {api.common._script_name()} --bam 1.bam 2.bam --bed in.bed > out.tsv
  $ fuc {api.common._script_name()} --fn bam.list --region chr1:100-200 > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Compute read depth from SAM/BAM/CRAM files.',
        description=description,
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more alignment files.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one alignment file per line.'
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help='BED file.'
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Target region ('chrom:start-end')."
    )
    parser.add_argument(
        '--zero',
        action='store_true',
        help='Output all positions including those with zero depth.'
    )

def main(args):
    cf = api.pycov.CovFrame.from_bam(
        bam=args.bam, fn=args.fn, bed=args.bed, region=args.region,
        zero=args.zero
    )
    sys.stdout.write(cf.to_string())
