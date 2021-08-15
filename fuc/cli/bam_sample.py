from .. import api

description = f"""
This command will extract the sample name of input SAM/BAM/CRAM file.

Usage examples:
  $ fuc {api.common._script_name()} in.sam
  $ fuc {api.common._script_name()} in.bam
  $ fuc {api.common._script_name()} in.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Extract the sample name of SAM/BAM/CRAM file.',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file.')

def main(args):
    print(api.pybam.tag_sm(args.bam)[0])
