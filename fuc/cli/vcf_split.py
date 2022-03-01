import sys
import os
import shutil

from .. import api

description = """
Split a VCF file by individual.
"""

epilog = f"""
[Example] Split a VCF file by individual:
  $ fuc {api.common._script_name()} in.vcf output_dir
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Split a VCF file by individual."""
    )
    parser.add_argument(
        'vcf',
        help=
"""VCF file to be split."""
    )
    parser.add_argument(
        'output',
        type=os.path.abspath,
        help=
"""Output directory."""
    )
    parser.add_argument(
        '--clean',
        action='store_false',
        help=
"""By default, the command will only return variants present in
each individual. Use the tag to stop this behavior and make
sure that all individuals have the same number of variants."""
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help=
"""Overwrite the output directory if it already exists."""
    )

def main(args):
    if os.path.exists(args.output) and args.force:
        shutil.rmtree(args.output)

    os.mkdir(args.output)

    vfs = api.pyvcf.split(args.vcf, clean=args.clean)

    for vf in vfs:
        vf.to_file(f'{args.output}/{vf.samples[0]}.vcf')
