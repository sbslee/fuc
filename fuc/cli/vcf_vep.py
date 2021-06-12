from .. import api
from argparse import RawTextHelpFormatter

command = api.common.script_name(__file__)

description = f"""
This command will filter a VCF file annotated by Ensemble VEP. It essentially wraps the `pandas.DataFrame.query` method. For details on query expression, please visit the method's documentation page (https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html#pandas-dataframe-query).

examples:
  $ fuc {command} in.vcf 'SYMBOL == "TP53"' > out.vcf
  $ fuc {command} in.vcf 'SYMBOL != "TP53"' > out.vcf
  $ fuc {command} in.vcf 'SYMBOL == "TP53"' --opposite > out.vcf
  $ fuc {command} in.vcf 'Consequence in ["splice_donor_variant", "stop_gained"]' > out.vcf
  $ fuc {command} in.vcf '(SYMBOL == "TP53") and (Consequence.str.contains("stop_gained"))' > out.vcf
  $ fuc {command} in.vcf 'gnomAD_AF < 0.001' > out.vcf
  $ fuc {command} in.vcf 'gnomAD_AF < 0.001' --as_zero > out.vcf
"""

def create_parser(subparsers):
    parser = subparsers.add_parser(
        command,
        help='[VCF] filter a VCF file annotated by Ensemble VEP',
        description=description,
        formatter_class=RawTextHelpFormatter,
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
    filtered_vf = api.pyvep.filter_query(vf, args.expr, opposite=args.opposite, as_zero=args.as_zero)
    print(filtered_vf.to_string())
