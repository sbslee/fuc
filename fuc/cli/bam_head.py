import sys

from .. import api

import pysam

description = """
Print the header of a BAM file.
"""

epilog = f"""
[Example] Print the header of a BAM file:
  $ fuc {api.common._script_name()} in.bam
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Print the header of a BAM file."""
    )
    parser.add_argument(
        'bam',
        help=
"""Input alignment file."""
    )

def main(args):
    header = pysam.view('-H', args.bam, '--no-PG')
    sys.stdout.write(header)
