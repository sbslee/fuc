import sys

from .. import api

description = """
Call SNVs and indels from BAM files.

Under the hood, the command utilizes the bcftool program to call variants.
"""

epilog = f"""
[Example] Specify regions manually:
  $ fuc {api.common._script_name()} ref.fa 1.bam 2.bam \\
  -r chr1:100-200 chr2:400-500 > out.vcf

[Example] Specify regions with a BED file:
  $ fuc {api.common._script_name()} ref.fa bam.list \\
  -r in.bed > out.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Call SNVs and indels from BAM files."""
    )
    parser.add_argument(
        'fasta',
        help=
"""Reference FASTA file."""
    )
    parser.add_argument(
        'bams',
        nargs='+',
        help=
"""One or more input BAM files. Alternatively, you can
provide a text file (.txt, .tsv, .csv, or .list)
containing one BAM file per line."""
    )
    parser.add_argument(
        '-r',
        '--regions',
        nargs='+',
        metavar='TEXT',
        help=
"""By default, the command looks at each genomic
position with coverage in BAM files, which can be
excruciatingly slow for large files (e.g. whole
genome sequencing). Therefore, use this argument to
only call variants in given regions. Each region must
have the format chrom:start-end and be a half-open
interval with (start, end]. This means, for example,
chr1:100-103 will extract positions 101, 102, and
103. Alternatively, you can provide a BED file
(compressed or uncompressed) to specify regions. Note
that the 'chr' prefix in contig names (e.g. 'chr1'
vs. '1') will be automatically added or removed as
necessary to match the input BAM's contig names."""
    )
    parser.add_argument(
        '--min-mq',
        metavar='INT',
        type=int,
        default=1,
        help=
"""Minimum mapping quality for an alignment to be used
(default: 1)."""
    )
    parser.add_argument(
        '--max-depth',
        metavar='INT',
        type=int,
        default=250,
        help=
"""At a position, read maximally this number of reads
per input file (default: 250)."""
    )
    parser.add_argument(
        '--dir-path',
        metavar='PATH',
        help=
"""By default, intermediate files (likelihoods.bcf,
calls.bcf, and calls.normalized.bcf) will be stored
in a temporary directory, which is automatically
deleted after creating final VCF. If you provide a
directory path, intermediate files will be stored
there."""
    )
    parser.add_argument(
        '--gap_frac',
        metavar='FLOAT',
        type=float,
        default=0.002,
        help=
"""Minimum fraction of gapped reads for calling indels
(default: 0.002)."""
    )
    parser.add_argument(
        '--group-samples',
        metavar='PATH',
        help=
"""By default, all samples are assumed to come from a
single population. This option allows to group
samples into populations and apply the HWE assumption
within but not across the populations. To use this
option, provide a tab-delimited text file with sample
names in the first column and group names in the
second column. If '--group-samples -' is given
instead, no HWE assumption is made at all and
single-sample calling is performed. Note that in low
coverage data this inflates the rate of false
positives. Therefore, make sure you know what you are
doing."""
    )

def main(args):
    api.pyvcf.call(
        args.fasta, args.bams, regions=args.regions, path='-',
        min_mq=args.min_mq, max_depth=args.max_depth, dir_path=args.dir_path,
        gap_frac=args.gap_frac, group_samples=args.group_samples
    )
