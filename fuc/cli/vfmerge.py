from fuc.api.common import get_script_name
from fuc import pyvcf

def create_parser(subparsers):
    parser = subparsers.add_parser(
        get_script_name(__file__),
        help='[VCF] merge two or more VCF files',
        description='This command will merge multiple VCF files (both zipped '
                    'and unzipped). By default, only the GT subfield of '
                    'the FORMAT field will be included in the merged VCF. '
                    "Use '--format' to include additional FORMAT "
                    'subfields such as AD and DP.'
    )
    parser.add_argument('vcf_files', help='VCF files', nargs='+')
    parser.add_argument('--how', metavar='TEXT', default='inner',
        help='type of merge as defined in `pandas.DataFrame.merge` '
              "(default: 'inner')")
    parser.add_argument('--format', metavar='TEXT', default='GT',
        help="FORMAT subfields to be retained (e.g. 'GT:AD:DP') "
        "(default: 'GT')"
    )

def main(args):
    vfs = [pyvcf.read_file(x) for x in args.vcf_files]
    merged_vf = pyvcf.merge(vfs, format=args.format, how=args.how)
    print(merged_vf.to_string())
