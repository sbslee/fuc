import sys

from .. import api

import pysam

description = """
Compute read depth from SAM/BAM/CRAM files.

By default, the command will count all reads within the alignment files. You
can specify regions of interest with --bed or --region. When you do this, pay
close attention to the 'chr' string in contig names (e.g. 'chr1' vs. '1').
Note also that --region requires the input files be indexed.
"""

epilog = f"""
[Example] To specify regions with a BED file:
  $ fuc {api.common._script_name()} \\
  --bam 1.bam 2.bam \\
  --bed in.bed > out.tsv

[Example] To specify regions manually:
  $ fuc {api.common._script_name()} \\
  --fn bam.list \\
  --region chr1:100-200 > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Compute read depth from SAM/BAM/CRAM files.',
    )
    parser.add_argument(
        '--bam',
        metavar='PATH',
        nargs='+',
        help='One or more alignment files. Cannot be used with --fn.'
    )
    parser.add_argument(
        '--fn',
        metavar='PATH',
        help='File containing one alignment file per line. Cannot \n'
             'be used with --bam.'
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help='BED file. Cannot be used with --region.'
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Target region ('chrom:start-end'). Cannot be used \n"
             "with --bed."
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
