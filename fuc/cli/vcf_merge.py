import sys

from .. import api

description = """
Merge two or more VCF files.
"""

epilog = f"""
[Example] Merge multiple VCF files:
  $ fuc {api.common._script_name()} 1.vcf 2.vcf 3.vcf > merged.vcf

[Example] Keep the GT, AD, DP fields:
  $ fuc {api.common._script_name()} 1.vcf 2.vcf --format GT:AD:DP > merged.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Merge two or more VCF files.',
    )
    parser.add_argument(
        'vcf_files',
        nargs='+',
        help=
"""VCF files (compressed or uncompressed). Note that the 'chr'
prefix in contig names (e.g. 'chr1' vs. '1') will be
automatically added or removed as necessary to match the
contig names of the first VCF."""
    )
    parser.add_argument(
        '--how',
        metavar='TEXT',
        default='inner',
        help=
"""Type of merge as defined in pandas.DataFrame.merge
(default: 'inner')."""
    )
    parser.add_argument(
        '--format',
        metavar='TEXT',
        default='GT',
        help=
"""FORMAT subfields to be retained (e.g. 'GT:AD:DP')
(default: 'GT')."""
    )
    parser.add_argument(
        '--sort',
        action='store_false',
        help=
"""Use this flag to turn off sorting of records
(default: True)."""
    )
    parser.add_argument(
        '--collapse',
        action='store_true',
        help=
"""Use this flag to collapse duplicate records
(default: False)."""
    )

def main(args):
    vfs = [api.pyvcf.VcfFrame.from_file(x) for x in args.vcf_files]
    merged_vf = api.pyvcf.merge(vfs, format=args.format, how=args.how,
        sort=args.sort, collapse=args.collapse)
    sys.stdout.write(merged_vf.to_string())
