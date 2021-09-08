import sys

from .. import api

description = f"""
###########################################
# Slice a VCF file for specified regions. #
###########################################

Target regions can be specified with either '--region', '--bed', or '--vcf'.

Pay attention to the 'chr' string in contig names (e.g. 'chr1' vs. '1').

Usage examples:
  $ fuc {api.common._script_name()} in.vcf --region 1 > out.vcf
  $ fuc {api.common._script_name()} in.vcf --region 1:100-300 > out.vcf
  $ fuc {api.common._script_name()} in.vcf --region 1:100 > out.vcf
  $ fuc {api.common._script_name()} in.vcf --region chr1:100- > out.vcf
  $ fuc {api.common._script_name()} in.vcf --region chr1:-300 > out.vcf
  $ fuc {api.common._script_name()} in.vcf --bed targets.bed > out.vcf
  $ fuc {api.common._script_name()} in.vcf --vcf targets.vcf > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Slice a VCF file for specified regions.',
        description=description,
    )
    parser.add_argument(
        'input',
        help='Input VCF file (zipped or unzipped).'
    )
    parser.add_argument(
        '--region',
        metavar='TEXT',
        help="Target region to use for slicing ('chrom:start-end')."
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        help="BED file to use for slicing (zipped or unzipped)."
    )
    parser.add_argument(
        '--vcf',
        metavar='PATH',
        help="VCF file to use for slicing (zipped or unzipped)."
    )

def main(args):
    vf = api.pyvcf.VcfFrame.from_file(args.input)

    if len([x for x in [args.region, args.bed, args.vcf]
        if x is not None]) > 1:
        raise ValueError('Too many arguments')

    if args.region is not None:
        sliced_vf = vf.slice(args.region)
    elif args.bed is not None:
        sliced_vf = vf.filter_bed(args.bed)
    elif args.vcf is not None:
        sliced_vf = vf.filter_vcf(args.vcf)
    else:
        raise ValueError('Missing target regions')

    sys.stdout.write(sliced_vf.to_string())
