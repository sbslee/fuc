import sys
import subprocess
import os

from .. import api

from Bio import SeqIO

description = """
Filter sequence records in a FASTA file.
"""

epilog = f"""
[Example] Select certain contigs:
  $ fuc {api.common._script_name()} in.fasta --contigs chr1 chr2 > out.fasta

[Example] Select certain contigs:
  $ fuc {api.common._script_name()} in.fasta --contigs contigs.list --exclude > out.fasta
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Filter sequence records in a FASTA file."""
    )
    parser.add_argument(
        'fasta',
        help=
"""Input FASTA file (compressed or uncompressed)."""
    )
    parser.add_argument(
        '--contigs',
        metavar='TEXT',
        nargs='+',
        help=
"""One or more contigs to be selected. Alternatively, you can
provide a file containing one contig per line."""
    )
    parser.add_argument(
        '--exclude',
        action='store_true',
        help=
"""Exclude specified contigs."""
    )

def main(args):
    if os.path.exists(args.contigs[0]):
        contigs = api.common.convert_file2list(args.contigs[0])
    else:
        contigs = args.contigs

    records = []

    for record in SeqIO.parse(args.fasta, 'fasta'):
        if args.exclude:
            if record.id not in contigs:
                records.append(record)
        else:
            if record.id in contigs:
                records.append(record)

    SeqIO.write(records, sys.stdout, 'fasta')
