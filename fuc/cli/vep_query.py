from .. import api
from argparse import RawTextHelpFormatter

command = api.common.script_name(__file__)

description = f"""
This command will filter a VEP-annotated VCF file. It essentially wraps the `pandas.DataFrame.query` method. For details on the method, please visit its documentation page (https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html#pandas-dataframe-query).

examples:
  $ fuc {command} in.vcf 'SYMBOL == "TP53"' > out.vcf
  $ fuc {command} in.vcf 'Consequence in ["splice_donor_variant", "stop_gained"]' > out.vcf
  $ fuc {command} in.vcf '(SYMBOL == "TP53") and (Consequence.str.contains("stop_gained"))' > out.vcf
"""

def create_parser(subparsers):
    parser = subparsers.add_parser(
        command,
        help='[VEP] filter a VEP-annotated VCF file',
        description=description,
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument('vcf', help='VEP-annotated VCF file')
    parser.add_argument('expr', help='query expression to evaluate')

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    df = api.pyvep.to_frame(vf)
    i = df.query(args.expr).index
    filtered_vf = api.pyvcf.VcfFrame(vf.copy_meta(), vf.df.iloc[i])
    print(filtered_vf.to_string())
