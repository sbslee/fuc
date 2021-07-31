import sys
import subprocess

from .. import api

description = f"""
This command will count sequence reads in FASTQ files.

It will look for stdin if there are no arguments.

Usage examples:
  $ fuc {api.common._script_name()} in.fastq
  $ cat fastq.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Count sequence reads in FASTQ files.',
        description=description,
    )
    parser.add_argument(
        'fastq',
        nargs='*',
        help='FASTQ files (zipped or unzipped) (default: stdin).'
    )
    return parser

def main(args):
    if args.fastq:
        fastqs = args.fastq
    elif not sys.stdin.isatty():
        fastqs = sys.stdin.read().rstrip('\n').split('\n')
    else:
        raise ValueError('No input files detected.')

    for fastq in fastqs:
        if fastq.endswith('.gz'):
            cat = 'zcat'
        else:
            cat = 'cat'
        command = f'echo $({cat} < {fastq} | wc -l) / 4 | bc'
        subprocess.run(command, shell=True)
