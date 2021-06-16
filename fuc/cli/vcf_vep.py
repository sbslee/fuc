from .. import api

description = f"""
This command will filter a VCF file annotated by Ensemble VEP. It
essentially wraps the 'pyvep.filter_query' method from the fuc API. For
details on query expression, please visit the method's documentation page.

usage examples:
  $ fuc {api.common._script_name()} in.vcf 'SYMBOL == "TP53"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'SYMBOL != "TP53"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'SYMBOL == "TP53"' --opposite > out.vcf
  $ fuc {api.common._script_name()} in.vcf \\
      'Consequence in ["splice_donor_variant", "stop_gained"]' > out.vcf
  $ fuc {api.common._script_name()} in.vcf \\
      '(SYMBOL == "TP53") and (Consequence.str.contains("stop_gained"))' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'gnomAD_AF < 0.001' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'gnomAD_AF < 0.001' --as_zero > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] filter a VCF file annotated by Ensemble VEP',
        description=description,
    )
    parser.add_argument('vcf', help='Ensemble VEP-annotated VCF file')
    parser.add_argument('expr', help='query expression to evaluate')
    parser.add_argument(
        '--opposite',
        action='store_true',
        help='use this flag to return records that don’t meet the said criteria'
    )
    parser.add_argument(
        '--as_zero',
        action='store_true',
        help='use this flag to treat missing values as zero instead of NaN'
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    filtered_vf = api.pyvep.filter_query(vf, args.expr,
        opposite=args.opposite, as_zero=args.as_zero)
    print(filtered_vf.to_string(), end='')
