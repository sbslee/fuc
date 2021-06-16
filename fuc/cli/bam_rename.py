from .. import api

description = f"""
This command will rename the samples in a SAM/BAM/CRAM file. It essentially
wraps the 'pybam.rename' method from the fuc API.

usage examples:
  $ fuc {api.common._script_name()} in.sam new_name out.sam
  $ fuc {api.common._script_name()} in.bam new_name out.bam
  $ fuc {api.common._script_name()} in.cram new_name out.cram
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] rename the samples in a SAM/BAM/CRAM file',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')
    parser.add_argument('name', help='sample name')
    parser.add_argument('out', help='output file')

def main(args):
    api.pybam.rename(args.bam, args.name, args.out)
