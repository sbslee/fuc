import sys
import subprocess

from .. import api

description = """
Count sequence reads in FASTQ files.
"""

epilog = f"""
[Example] When the input is a FASTQ file:
  $ fuc {api.common._script_name()} in1.fastq in2.fastq

[Example] When the input is stdin:
  $ cat fastq.list | fuc {api.common._script_name()}
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Count sequence reads in FASTQ files."""
    )
    parser.add_argument(
        'fastq',
        nargs='*',
        help=
"""Input FASTQ files (compressed or uncompressed) (default: stdin)."""
    )

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
