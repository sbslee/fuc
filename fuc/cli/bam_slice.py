import sys

from .. import api

import pysam

description = f"""
This command will slice the input SAM/BAM/CRAM file for the specified region(s).

Usage examples:
  $ fuc {api.common._script_name()} in.bam chr1:100-200 > out.bam
  $ fuc {api.common._script_name()} in.bam chr1:100-200 chr2:100-200 > out.bam
  $ fuc {api.common._script_name()} in.bam chr1:100-200 --format SAM > out.sam
  $ fuc {api.common._script_name()} in.bam chr1:100-200 --format CRAM --fasta ref.fa > out.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Slice a SAM/BAM/CRAM file.',
        description=description,
    )
    parser.add_argument(
        'bam',
        help='SAM/BAM/CRAM file.'
    )
    parser.add_argument(
        'region',
        nargs='+',
        help="Space-separated regions ('chrom:start-end')."
    )
    parser.add_argument(
        '--format',
        metavar='TEXT',
        default='BAM',
        choices=['SAM', 'BAM', 'CRAM'],
        help=
            "Output format (default: 'BAM') (choices: 'SAM', 'BAM', 'CRAM'). "
            "A FASTA file must be specified with '--fasta' for 'CRAM'."
    )
    parser.add_argument(
        '--fasta',
        metavar='PATH',
        help="FASTA file. Required when '--format' is 'CRAM'."
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

    alignments = pysam.view(args.bam, *args.region, *options)

    stdout_method(alignments)
