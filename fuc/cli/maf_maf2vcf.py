from .. import api
from argparse import RawTextHelpFormatter

command = api.common._script_name(__file__)

description = f"""
This command will convert a MAF file to a VCF file. It essentially wraps the `pymaf.MafFrame.to_vcf` method. For details on the conversion algorithm, please visit the method's documentation page (https://sbslee-fuc.readthedocs.io/en/latest/api.html#fuc.api.pymaf.MafFrame.to_vcf).

examples:
  $ fuc {command} in.maf hs37d5.fa > out.vcf
  $ fuc {command} in.maf --ignore_indels > out.vcf
"""

def create_parser(subparsers):
    parser = subparsers.add_parser(
        command,
        help='[MAF] convert a MAF file to a VCF file',
        description=description,
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument('maf', help='MAF file')
    parser.add_argument('--fasta', metavar='PATH', help='FASTA file (required to include INDELs)')
    parser.add_argument('--ignore_indels', action='store_true', help='use this tag to exclude INDELs from the output')

def main(args):
    mf = api.pymaf.MafFrame.from_file(args.maf)
    vf = mf.to_vcf(fasta=args.fasta, ignore_indels=args.ignore_indels)
    print(vf.to_string())
