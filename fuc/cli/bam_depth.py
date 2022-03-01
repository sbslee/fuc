import sys

from .. import api

import pysam

description = """
Compute per-base read depth from BAM files.

Under the hood, the command computes read depth using the 'samtools depth'
command.
"""

epilog = f"""
[Example] Specify regions manually:
  $ fuc {api.common._script_name()} 1.bam 2.bam \\
  -r chr1:100-200 chr2:400-500 > out.tsv

[Example] Specify regions with a BED file:
  $ fuc {api.common._script_name()} bam.list \\
  -r in.bed > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compute per-base read depth from BAM files."""
    )
    parser.add_argument(
        'bams',
        nargs='+',
        help=
"""One or more input BAM files. Alternatively, you can
provide a text file (.txt, .tsv, .csv, or .list)
containing one BAM file per line."""
    )
    parser.add_argument(
        '-r',
        '--regions',
        nargs='+',
        metavar='TEXT',
        help=
"""By default, the command counts all reads in BAM
files, which can be excruciatingly slow for large
files (e.g. whole genome sequencing). Therefore, use
this argument to only output positions in given
regions. Each region must have the format
chrom:start-end and be a half-open interval with
(start, end]. This means, for example, chr1:100-103
will extract positions 101, 102, and 103.
Alternatively, you can provide a BED file (compressed
or uncompressed) to specify regions. Note that the
'chr' prefix in contig names (e.g. 'chr1' vs. '1')
will be automatically added or removed as necessary
to match the input BAM's contig names."""
    )
    parser.add_argument(
        '--zero',
        action='store_true',
        help=
"""Output all positions including those with zero depth."""
    )

def main(args):
    cf = api.pycov.CovFrame.from_bam(
        args.bams, regions=args.regions, zero=args.zero
    )
    sys.stdout.write(cf.to_string())
