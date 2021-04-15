from fuc.api.common import get_script_name
from fuc.api.VcfFrame import VcfFrame

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[VCF] merge two or more VCF files',
        description='This command will merge multiple VCF files (both zipped '
                    'and unzipped). By default, only the GT subfield of '
                    'the FORMAT field will be included in the merged VCF. '
                    "Use '--format_subfields' to include additional FORMAT "
                    'subfields such as AD and DP.'
    )
    parser.add_argument('input_vcf', help='input VCF files', nargs='+')
    parser.add_argument('--format_subfields', help='FORMAT subfields',
        nargs='+')

def main(args):
    vcf_list = [VcfFrame.from_file(x) for x in args.input_vcf]
    merged_vcf = vcf_list[0]
    for vcf in vcf_list[1:]:
        merged_vcf = merged_vcf.merge(vcf,
            format_subfields=args.format_subfields)
    print(merged_vcf.to_string())
