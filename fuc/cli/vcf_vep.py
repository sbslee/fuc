import sys

from .. import api

description = """
Filter a VCF file by annotations from Ensembl VEP.
"""

epilog = f"""
[Example] Select variants in the TP53 gene:
  $ fuc {api.common._script_name()} in.vcf "SYMBOL == 'TP53'" > out.vcf

[Example] Exclude variants from the TP53 gene:
  $ fuc {api.common._script_name()} in.vcf "SYMBOL != 'TP53'" > out.vcf

[Example] Same as above:
  $ fuc {api.common._script_name()} in.vcf "SYMBOL == 'TP53'" --opposite > out.vcf

[Example] Select splice donor or stop-gain variants:
  $ fuc {api.common._script_name()} in.vcf \\
  "Consequence in ['splice_donor_variant', 'stop_gained']" > out.vcf

[Example] Build a complex query to select specific variants:
  $ fuc {api.common._script_name()} in.vcf \\
  "(SYMBOL == 'TP53') and (Consequence.str.contains('stop_gained'))" > out.vcf

[Example] Select variants whose gnomAD AF is less than 0.001:
  $ fuc {api.common._script_name()} in.vcf "gnomAD_AF < 0.001" > out.vcf

[Example] Variants without AF data will be treated as having AF of 0:
  $ fuc {api.common._script_name()} in.vcf "gnomAD_AF < 0.001" --as_zero > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Filter a VCF file by annotations from Ensembl VEP."""
    )
    parser.add_argument(
        'vcf',
        help=
"""VCF file annotated by Ensembl VEP (compressed or uncompressed)."""
    )
    parser.add_argument(
        'expr',
        help=
"""Query expression to evaluate."""
    )
    parser.add_argument(
        '--opposite',
        action='store_true',
        help=
"""Use this flag to return only records that don't
meet the said criteria."""
    )
    parser.add_argument(
        '--as_zero',
        action='store_true',
        help=
"""Use this flag to treat missing values as zero instead of NaN."""
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    filtered_vf = api.pyvep.filter_query(vf, args.expr,
        opposite=args.opposite, as_zero=args.as_zero)
    sys.stdout.write(filtered_vf.to_string())
