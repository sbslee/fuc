import sys

from .. import api

description = """
Convert a VCF file to a BED file.
"""

epilog = f"""
[Example] Convert VCF to BED:
  $ fuc {api.common._script_name()} in.vcf > out.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Convert a VCF file to a BED file."""
    )
    parser.add_argument(
        'vcf',
        help=
"""VCF file (compressed or uncompressed)."""
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    bf = vf.to_bed()
    sys.stdout.write(bf.to_string())
