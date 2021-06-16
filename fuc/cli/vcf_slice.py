from .. import api

description = f"""
This command will slice a VCF file (both zipped and unzipped). It essentially
wraps the 'pyvcf.VcfFrame.slice' method from the fuc API.

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
        help='[VCF] slice a VCF file',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')
    parser.add_argument('region', help="region ('chrom:start-end')")

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    sliced_vf = vf.slice(args.region)
    print(sliced_vf.to_string(), end='')
