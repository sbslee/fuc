from .. import api
import pysam

description = f"""
This command will slice a SAM/BAM/CRAM file. It essentially wraps the
'pysam.view' method from the pysam package.

By default, the command will index the output file. Use the '--no_index' flag
to skip indexing.

usage examples:
  $ fuc {api.common._script_name()} in.sam 4:300-400 out.sam
  $ fuc {api.common._script_name()} in.bam chr1:100-200 out.bam
  $ fuc {api.common._script_name()} in.cram chr1:100-200 out.cram --no_index
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[BAM] slice a SAM/BAM/CRAM file',
        description=description,
    )
    parser.add_argument('bam', help='SAM/BAM/CRAM file')
    parser.add_argument('region', help="region ('chrom:start-end')")
    parser.add_argument('out', help='output file')
    parser.add_argument('--no_index', action='store_true',
        help='use this flag to skip indexing')

def main(args):
    b = pysam.view('-b', args.bam, args.region)
    with open(args.out, 'wb') as f:
        f.write(b)
    if not args.no_index:
        pysam.index(args.out)
