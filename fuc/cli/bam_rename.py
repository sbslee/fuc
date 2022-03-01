import sys
import tempfile

from .. import api

import pysam

description = """
Rename the sample in a BAM file.
"""

epilog = f"""
[Example] Write a new BAM file after renaming:
  $ fuc {api.common._script_name()} in.bam NA12878 > out.bam
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Rename the sample in a BAM file."""
    )
    parser.add_argument(
        'bam',
        help=
"""Input alignment file."""
    )
    parser.add_argument(
        'name',
        help=
"""New sample name."""
    )

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
