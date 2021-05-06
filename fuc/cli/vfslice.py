from fuc.api.common import get_script_name
from fuc import pyvcf

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[VCF] slice a VCF file',
        description='This command will slice a VCF file (both zipped '
                    'and unzipped).'
    )
    parser.add_argument('vcf_file', help='VCF file')
    parser.add_argument('chrom', help='chromosome')
    parser.add_argument('--start', metavar='INTEGER', type=int,
        help='start position')
    parser.add_argument('--end', metavar='INTEGER', type=int,
        help='end position')

def main(args):
    vf = pyvcf.read_file(args.vcf_file)
    sliced_vf = vf.slice(args.chrom, start=args.start, end=args.end)
    print(sliced_vf.to_string())
