import sys

from .. import api

import pysam

description = """
Slice an alignment file (SAM/BAM/CRAM).
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
        help='Slice a SAM/BAM/CRAM file.',
    )
    parser.add_argument(
        'bam',
        help="Input alignment file must be already indexed (.bai) to allow \n"
             "random access. You can index an alignment file with the \n"
             "bam-index command."
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help="One or more regions to be sliced. Each region must have the \n"
             "format chrom:start-end and be a half-open interval with \n"
             "(start, end]. This means, for example, chr1:100-103 will \n"
             "extract positions 101, 102, and 103. Alternatively, you can \n"
             "provide a BED file (compressed or uncompressed) to specify \n"
             "regions. Note that the 'chr' prefix in contig names (e.g. \n"
             "'chr1' vs. '1') will be automatically added or removed as \n"
             "necessary to match the input BED's contig names."
    )
    parser.add_argument(
        '--format',
        metavar='TEXT',
        default='BAM',
        choices=['SAM', 'BAM', 'CRAM'],
        help="Output format (default: 'BAM') (choices: 'SAM', 'BAM', \n"
             "'CRAM')."
    )
    parser.add_argument(
        '--fasta',
        metavar='PATH',
        help="FASTA file. Required when --format is 'CRAM'."
    )

def main(args):
    options = ['-h', '--no-PG']

    # Determine the output format.
    if args.format == 'BAM':
        stdout_method = sys.stdout.buffer.write
        options += ['-b']
    elif args.format == 'CRAM':
        stdout_method = sys.stdout.buffer.write
        if args.fasta is None:
            raise ValueError(
                "A FASTA file must be specified with '--fasta' "
                "when '--format' is 'CRAM'."
            )
        options += ['-C', '-T', args.fasta]
    else:
        stdout_method = sys.stdout.write

    # Parse the regions argument.
    if '.bed' in args.regions[0]:
        args.regions = api.pybed.BedFrame.from_file(args.regions[0]).to_regions()
    else:
        args.regions = api.common.sort_regions(args.regions)

    if api.pybam.has_chr_prefix(args.bam):
        args.regions = api.common.update_chr_prefix(args.regions, mode='add')
    else:
        args.regions = api.common.update_chr_prefix(args.regions, mode='remove')

    alignments = pysam.view(args.bam, *args.regions, *options)

    stdout_method(alignments)
