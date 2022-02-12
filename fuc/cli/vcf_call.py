import sys

from .. import api

description = f"""
Call variants (SNVs/indels) from BAM files.

This command will run a fully customizable, bcftools-based pipeline for
calling and filtering variants.
"""

epilog = f"""
[Example] Specify regions manually:
  $ fuc {api.common._script_name()} ref.fa in1.bam in2.bam -r chr1:100-200 chr2:300-400 > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Call variants (SNVs/indels) from BAM files.',
    )
    parser.add_argument(
        'fasta',
        help="Reference FASTA file."
    )
    parser.add_argument(
        'bams',
        nargs='+',
        help="One or more input BAM files. Alternatively, you can \n"
             "provide a text file (.txt, .tsv, .csv, or .list) \n"
             "containing one BAM file per line."
    )
    parser.add_argument(
        '--regions',
        nargs='+',
        metavar='TEXT',
        help="Only call variants in given regions. Each region must \n"
             "have the format chrom:start-end and be a half-open \n"
             "interval with (start, end]. This means, for example, \n"
             "chr1:100-103 will extract positions 101, 102, and \n"
             "103. Alternatively, you can provide a BED file \n"
             "(compressed or uncompressed) to specify regions. Note \n"
             "that the 'chr' prefix in contig names (e.g. 'chr1' \n"
             "vs. '1') will be automatically added or removed as \n"
             "necessary to match the input VCF's contig names."
    )
    parser.add_argument(
        '--min-mq',
        metavar='INT',
        type=int,
        default=1,
        help="Minimum mapping quality for an alignment to be used \n"
             "(default: 1)."
    )
    parser.add_argument(
        '--max-depth',
        metavar='INT',
        type=int,
        default=250,
        help="At a position, read maximally this number of reads \n"
             "per input file (default: 250)."
    )

def main(args):
    api.pyvcf.call(
        args.fasta, args.bams, path='-', regions=args.regions,
        min_mq=args.min_mq, max_depth=args.max_depth
    )
