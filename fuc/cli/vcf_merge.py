from .. import api

description = f"""
This command will merge multiple VCF files (both zipped and unzipped). It
essentially wraps the 'pyvcf.merge' method from the fuc API.

By default, only the GT subfield of the FORMAT field will be included in the
merged VCF. Use '--format' to include additional FORMAT subfields such as AD
and DP.

usage examples:
  $ fuc {api.common._script_name()} 1.vcf 2.vcf 3.vcf > merged.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='[VCF] merge two or more VCF files',
        description=description,
    )
    parser.add_argument('vcf_files', help='VCF files', nargs='+')
    parser.add_argument('--how', metavar='TEXT', default='inner',
        help='type of merge as defined in `pandas.DataFrame.merge` '
              "(default: 'inner')")
    parser.add_argument('--format', metavar='TEXT', default='GT',
        help="FORMAT subfields to be retained (e.g. 'GT:AD:DP') "
        "(default: 'GT')"
    )
    parser.add_argument('--sort', action='store_false',
        help='use this flag to turn off sorting of records (default: True)'
    )
    parser.add_argument('--collapse', action='store_true',
        help='use this flag to collapse duplicate records (default: False)'
    )

def main(args):
    vfs = [api.pyvcf.VcfFrame.from_file(x) for x in args.vcf_files]
    merged_vf = api.pyvcf.merge(vfs, format=args.format, how=args.how,
        sort=args.sort, collapse=args.collapse)
    print(merged_vf.to_string(), end='')
