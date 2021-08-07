import sys

from .. import api

description = f"""
This command will convert a MAF file to a sorted VCF file.

In order to handle INDELs the command makes use of a reference assembly (i.e. FASTA file). If SNVs are your only concern, then you do not need a FASTA file and can just use the '--ignore_indels' flag.

If you are going to provide a FASTA file, please make sure to select the appropriate one (e.g. one that matches the genome assembly).

In addition to basic genotype calls (e.g. '0/1'), you can extract more information from the MAF file by specifying the column(s) that contain additional genotype data of interest with the '--cols' argument. If provided, this argument will append the requested data to individual sample genotypes (e.g. '0/1:0.23').

You can also control how these additional genotype information appear in the FORMAT field (e.g. AF) with the '--names' argument. If this argument is not provided, the original column name(s) will be displayed.

Usage examples:
  $ fuc {api.common._script_name()} in.maf --fasta hs37d5.fa > out.vcf
  $ fuc {api.common._script_name()} in.maf --ignore_indels > out.vcf
  $ fuc {api.common._script_name()} in.maf --fasta hs37d5.fa --cols i_TumorVAF_WU --names AF > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Convert a MAF file to a VCF file.',
        description=description,
    )
    parser.add_argument(
        'maf',
        help='MAF file (zipped or unzipped).'
    )
    parser.add_argument(
        '--fasta',
        metavar='PATH',
        help='FASTA file (required to include INDELs in the output).'
    )
    parser.add_argument(
        '--ignore_indels',
        action='store_true',
        help='Use this flag to exclude INDELs from the output.'
    )
    parser.add_argument(
        '--cols',
        metavar='TEXT',
        nargs='+',
        help='Column(s) in the MAF file.'
    )
    parser.add_argument(
        '--names',
        metavar='TEXT',
        nargs='+',
        help='Name(s) to be displayed in the FORMAT field.'
    )

def main(args):
    mf = api.pymaf.MafFrame.from_file(args.maf)
    vf = mf.to_vcf(
        fasta=args.fasta, ignore_indels=args.ignore_indels,
        cols=args.cols, names=args.names
    )
    sys.stdout.write(vf.to_string())
