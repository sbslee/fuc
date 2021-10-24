from .. import api

import pysam

description = """
Index a SAM/BAM/CRAM file.
"""

epilog = f"""
[Example] Index a BAM file:
  $ fuc {api.common._script_name()} in.bam
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Index a SAM/BAM/CRAM file.',
    )
    parser.add_argument(
        'bam',
        help='Alignment file.'
    )

def main(args):
    pysam.index(args.bam)
