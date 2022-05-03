import sys

from .. import api

description = """
Slice a BAM file.
"""

epilog = f"""
[Example] Specify regions manually:
  $ fuc {api.common._script_name()} in.bam 1:100-300 2:400-700 > out.bam

[Example] Speicfy regions with a BED file:
  $ fuc {api.common._script_name()} in.bam regions.bed > out.bam

[Example] Slice a CRAM file:
  $ fuc {api.common._script_name()} in.bam regions.bed --format CRAM --fasta ref.fa > out.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Slice a BAM file."""
    )
    parser.add_argument(
        'bam',
        help=
"""Input BAM file. It must be already indexed to allow random
access. You can index a BAM file with the bam-index command."""
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help=
"""One or more regions to be sliced. Each region must have the
format chrom:start-end and be a half-open interval with
(start, end]. This means, for example, chr1:100-103 will
extract positions 101, 102, and 103. Alternatively, you can
provide a BED file (compressed or uncompressed) to specify
regions. Note that the 'chr' prefix in contig names (e.g.
'chr1' vs. '1') will be automatically added or removed as
necessary to match the input BED's contig names."""
    )
    parser.add_argument(
        '--format',
        metavar='TEXT',
        default='BAM',
        choices=['SAM', 'BAM', 'CRAM'],
        help=
"""Output format (default: 'BAM') (choices: 'SAM', 'BAM',
'CRAM')."""
    )
    parser.add_argument(
        '--fasta',
        metavar='PATH',
        help=
"""FASTA file. Required when --format is 'CRAM'."""
    )

def main(args):
    api.pybam.slice(args.bam, args.regions, format=args.format,
        path='-', fasta=args.fasta)
