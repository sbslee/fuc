import sys
import tempfile

from .. import api

import pysam

description = f"""
This command will rename the sample(s) in the input SAM/BAM/CRAM file.

Usage examples:
  $ fuc {api.common._script_name()} in.sam NA12878 > out.sam
  $ fuc {api.common._script_name()} in.bam NA12878 > out.bam
  $ fuc {api.common._script_name()} in.cram NA12878 > out.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Rename the samples in a SAM/BAM/CRAM file.',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file.')
    parser.add_argument('name', help='New sample name.')

def main(args):
    # Detect the input format.
    try:
        with open(args.bam, 'r') as f:
            for line in f:
                is_sam = True
                break
    except UnicodeDecodeError:
        is_sam = False

    # Rename the sample(s).
    old_lines = pysam.view('-H', args.bam, '--no-PG').strip().split('\n')
    new_lines = []
    for old_line in old_lines:
        old_fields = old_line.split('\t')
        if old_fields[0] != '@RG':
            new_lines.append(old_line)
            continue
        new_fields = []
        for old_field in old_fields:
            if 'SM:' not in old_field:
                new_fields.append(old_field)
                continue
            new_fields.append(f'SM:{args.name}')
        new_lines.append('\t'.join(new_fields))

    # Write the output file.
    if is_sam:
        for new_line in new_lines:
            sys.stdout.write(new_line + '\n')
        alignments = pysam.view(args.bam, '--no-PG')
        sys.stdout.write(alignments)
    else:
        with tempfile.TemporaryDirectory() as t:
            with open(f'{t}/header.sam', 'w') as f:
                f.write('\n'.join(new_lines))
            alignments = pysam.reheader(f'{t}/header.sam',
                                        args.bam,
                                        '--no-PG')
            sys.stdout.buffer.write(alignments)
