from fuc.api.common import get_script_name
from fuc.api.VCFResult import VCFResult

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='merge multiple VCF files',
        description='This command merges multiple VCF files (both zipped '
                    'and unzipped). By default, only GT subfield of FORMAT '
                    "field is included. Use '--subfield' to include "
                    'additional subfields such as AD and DP.'
    )
    parser.add_argument('input_vcf', help='input VCF files', nargs='+')
    parser.add_argument('output_vcf', help='output VCF file')
    parser.add_argument('--subfield', help='FORMAT subfields', nargs='+')

def main(args):
    vcf_list = [VCFResult.read(x) for x in args.input_vcf]
    merged_vcf = VCFResult()
    for vcf in vcf_list:
        merged_vcf = merged_vcf.merge(vcf, args.subfield)
    merged_vcf.write(args.output_vcf)
