from .. import api

description = f"""
This command will filter a VCF file (both zipped and unzipped). It essentially
wraps the 'pyvcf.VcfFrame.markmiss' method from the fuc API.

usage examples:
  $ fuc {api.common._script_name()} in.vcf chr1 > sliced.vcf
  $ fuc {api.common._script_name()} in.vcf chr1:100-300 > sliced.vcf
  $ fuc {api.common._script_name()} in.vcf chr1:100 > sliced.vcf
  $ fuc {api.common._script_name()} in.vcf chr1:100- > sliced.vcf
  $ fuc {api.common._script_name()} in.vcf chr1:-300 > sliced.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] filter a VCF file',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')
    parser.add_argument('expr', help='expression to evaluate')
    parser.add_argument('--greedy', action='store_true', help='use this flag to mark even ambiguous genotypes as missing')
    parser.add_argument('--opposite', action='store_true', help='use this flag to mark all genotypes that do not satisfy the query expression as missing and keep those that do')
    parser.add_argument('--samples', metavar='PATH', help='file of sample names to apply the marking (one sample per line)')

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    filtered_vf = vf.markmiss(args.expr, greedy=args.greedy, opposite=args.opposite)
    print(filtered_vf.to_string(), end='')
