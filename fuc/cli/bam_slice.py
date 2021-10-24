import sys

from .. import api

import pysam

description = """
Slice a SAM/BAM/CRAM file.
"""

epilog = f"""
[Example] Slice a BAM file:
  $ fuc {api.common._script_name()} in.bam chr1:100-200 chr2:100-200 > out.bam

[Example] Slice a CRAM file:
  $ fuc {api.common._script_name()} in.bam chr1:100-200 --format CRAM --fasta ref.fa > out.cram
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
        help='Alignment file.'
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help="List of regions to be sliced ('chrom:start-end')."
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

    alignments = pysam.view(args.bam, *args.regions, *options)

    stdout_method(alignments)
