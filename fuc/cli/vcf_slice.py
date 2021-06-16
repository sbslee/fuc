from .. import api

description = f"""
This command will slice a VCF file (both zipped and unzipped).

usage examples:
  $ fuc {api.common._script_name()} in.vcf chr4 --start 300 --end 400 > sliced.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] slice a VCF file',
        description=description,
    )
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('chrom', help='chromosome')
    parser.add_argument('--start', metavar='INT', type=int,
        help='start position')
    parser.add_argument('--end', metavar='INT', type=int,
        help='end position')

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf_file)
    sliced_vf = vf.slice(args.chrom, start=args.start, end=args.end)
    print(sliced_vf.to_string())
