from .. import api

description = f"""
This command will compute intersections beween multiple BED files. It
essentially wraps the 'pybed.BedFrame.intersect' method from the fuc API.

usage examples:
  $ fuc {api.common._script_name()} 1.bed 2.bed 3.bed > intersect.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BED] find intersection of two or more BED files',
        description=description,
    )
    parser.add_argument('bed', help='BED files', nargs='+')

def main(args):
    bfs = [api.pybed.BedFrame.from_file(x) for x in args.bed]
    final_bf = bfs[0]
    for bf in bfs[1:]:
        final_bf = final_bf.intersect(bf)
    print(final_bf.to_string(), end='')
