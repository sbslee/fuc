from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[MAF] convert an annotated VCF file to a MAF file',
        description=
            'This command will convert an annotated VCF '
            'file to a MAF file.'
    )
    parser.add_argument('vcf_file', help='annotated VCF file')

def main(args):
    mf = api.pymaf.MafFrame.from_vcf(args.vcf_file)
    print(mf.to_string())
