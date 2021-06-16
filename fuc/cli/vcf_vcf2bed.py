from .. import api

description = f"""
This command will convert a VCF file to a BED file. It essentially wraps the
'pyvcf.VcfFrame.to_bed' method from the fuc API.

usage examples:
  $ fuc {api.common._script_name()} in.vcf > out.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] convert a VCF file to a BED file',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    bf = vf.to_bed()
    print(bf.to_string(), end='')
