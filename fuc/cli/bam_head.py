from .. import api

description = f"""
This command will print the header of a SAM/BAM/CRAM file. It essentially
wraps the 'pybam.header' method.

usage examples:
  $ fuc {api.common._script_name()} in.sam
  $ fuc {api.common._script_name()} in.bam
  $ fuc {api.common._script_name()} in.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] print the header of a SAM/BAM/CRAM file',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')

def main(args):
    print(api.pybam.header(args.bam))
