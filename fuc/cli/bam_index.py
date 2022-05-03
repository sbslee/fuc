from .. import api

description = """
Index a BAM file.
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
        help=
"""Index a BAM file."""
    )
    parser.add_argument(
        'bam',
        help=
"""Input alignment file."""
    )

def main(args):
    api.pybam.index(args.bam)
