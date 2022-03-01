import sys

from .. import api

description = """
Find the intersection of BED files.
"""

epilog = f"""
[Example] Find the intersection of three BED files:
  $ fuc {api.common._script_name()} in1.bed in2.bed in3.bed > out.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Find the intersection of BED files."""
    )
    parser.add_argument(
        'bed',
        nargs='+',
        help=
"""Input BED files."""
    )

def main(args):
    bfs = [api.pybed.BedFrame.from_file(x) for x in args.bed]
    final_bf = bfs[0]
    for bf in bfs[1:]:
        final_bf = final_bf.intersect(bf)
    sys.stdout.write(final_bf.to_string())
