import sys

from .. import api

description = """
Convert a VCF file to a MAF file.
"""

epilog = f"""
[Example] Convert VCF to MAF:
  $ fuc {api.common._script_name()} in.vcf > out.maf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Convert a VCF file to a MAF file."""
    )
    parser.add_argument(
        'vcf',
        help=
"""Annotated VCF file."""
    )

def main(args):
    mf = api.pymaf.MafFrame.from_vcf(args.vcf)
    sys.stdout.write(mf.to_string())
