import sys
import tempfile

from .. import api

import pysam

description = f"""
This command will rename the sample(s) in the input BAM/CRAM file.

usage examples:
  $ fuc {api.common._script_name()} in.bam NA12878 > out.bam
  $ fuc {api.common._script_name()} in.cram NA12878 > out.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] Rename the samples in a BAM/CRAM file.',
        description=description,
    )
    parser.add_argument('bam', help='BAM/CRAM file.')
    parser.add_argument('name', help='New sample name.')

def main(args):
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

    with tempfile.TemporaryDirectory() as t:
        with open(f'{t}/header.sam', 'w') as f:
            f.write('\n'.join(new_lines))

        alignments = pysam.reheader(f'{t}/header.sam', args.bam, '--no-PG')
        sys.stdout.buffer.write(alignments)
