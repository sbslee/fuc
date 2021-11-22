import sys

from .. import api

description = f"""
Slice a VCF file for specified regions.
"""

epilog = f"""
[Example] Specify regions manually:
  $ fuc {api.common._script_name()} in.vcf.gz 1:100-300 2:400-700 > out.vcf

[Example] Speicfy regions with a BED file:
  $ fuc {api.common._script_name()} in.vcf.gz regions.bed > out.vcf

[Example] Output a compressed file:
  $ fuc {api.common._script_name()} in.vcf.gz regions.bed | fuc fuc-bgzip > out.vcf.gz
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Slice a VCF file for specified regions.',
    )
    parser.add_argument(
        'vcf',
        help="Input VCF file must be already BGZF compressed (.gz) and \n"
             "indexed (.tbi) to allow random access. A VCF file can be \n"
             "compressed with the fuc-bgzip command and indexed with the \n"
             "vcf-index command."
    )
    parser.add_argument(
        'regions',
        nargs='+',
        help="One or more regions to be sliced. Each region must have the \n"
             "format chrom:start-end and be a half-open interval with \n"
             "(start, end]. This means, for example, chr1:100-103 will \n"
             "extract positions 101, 102, and 103. Alternatively, you can \n"
             "provide a BED file (compressed or uncompressed) to specify \n"
             "regions. Note that the 'chr' prefix in contig names (e.g. \n"
             "'chr1' vs. '1') will be automatically added or removed as \n"
             "necessary to match the input VCF's contig names."
    )

def main(args):
    api.pyvcf.slice(args.vcf, args.regions, path='-')
