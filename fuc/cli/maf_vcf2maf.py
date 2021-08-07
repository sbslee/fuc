import sys

from .. import api

description = f"""
This command will convert an annotated VCF file to a MAF file.

Usage examples:
  $ fuc {api.common._script_name()} in.vcf > out.maf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Convert a VCF file to a MAF file.',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file.')

def main(args):
    mf = api.pymaf.MafFrame.from_vcf(args.vcf)
    sys.stdout.write(mf.to_string())
