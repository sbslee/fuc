from .. import api
import pandas as pd

description = f"""
This command will filter a VCF file (both zipped and unzipped). It essentially
wraps the 'pyvcf.VcfFrame.markmiss' method from the fuc API.

usage examples:
  $ fuc {api.common._script_name()} in.vcf 'GT == "0/0"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'GT != "0/0"' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'DP < 30' > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'AD[1] < 10' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'AD[1] < 10 and DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'AD[1] < 10 or DP < 30' --greedy > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'AD[1] < 10 or DP < 30' --opposite > out.vcf
  $ fuc {api.common._script_name()} in.vcf 'np.mean(AD) < 10' --greedy --samples sample.list > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] filter a VCF file',
        description=description,
    )
    parser.add_argument('vcf', help='VCF file')
    parser.add_argument('expr', help='expression to evaluate')
    parser.add_argument(
        '--greedy',
        action='store_true',
        help='use this flag to mark even ambiguous genotypes as missing'
    )
    parser.add_argument(
        '--opposite',
        action='store_true',
        help='use this flag to mark all genotypes that do not satisfy the '
             'query expression as missing and leave those that do intact'
    )
    parser.add_argument(
        '--samples',
        metavar='PATH',
        help='file of sample names to apply the marking (one sample per line)'
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.vcf)
    if args.samples is None:
        samples = None
    else:
        samples = pd.read_table(args.samples, header=None)[0].to_list()
    filtered_vf = vf.markmiss(args.expr,
                              greedy=args.greedy,
                              opposite=args.opposite,
                              samples=samples)
    print(filtered_vf.to_string(), end='')
