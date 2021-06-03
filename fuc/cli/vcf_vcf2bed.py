from .. import api

def create_parser(subparsers):
    parser = subparsers.add_parser(
        api.common.script_name(__file__),
        help='[VCF] convert a VCF file to a BED file',
        description='This command will convert a VCF file to a BED file.'
    )
    parser.add_argument('vcf_file', help='VCF file')

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf_file)
    bf = vf.to_bed()
    print(bf.to_string(), end='')
