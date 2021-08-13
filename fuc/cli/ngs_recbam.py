import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline that marks duplicate reads and recalibrates BAM files.

External dependencies:
  - GATK: Required for marking duplicate reads and recalibrating BAM files.
  - samtools: Required for indexing BAM files.

Manifest columns:
  - BAM: Sorted BAM file.

Usage examples:
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-q queue_name" "-Xmx4g -Xms4g" 1.vcf 2.vcf 3.vcf
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-l h='node_A|node_B'" "-Xmx4g -Xms4g" 1.vcf 2.vcf 3.vcf
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Pipeline for marking duplicate reads and recalibrating BAM files.',
        description=description,
    )
    parser.add_argument(
        'manifest',
        help='Sample manifest CSV file.'
    )
    parser.add_argument(
        'fasta',
        help='Reference FASTA file.'
    )
    parser.add_argument(
        'output',
        type=os.path.abspath,
        help='Output directory.'
    )
    parser.add_argument(
        'qsub',
        type=str,
        help='Options for qsub.'
    )
    parser.add_argument(
        'java',
        help='Options for Java.'
    )
    parser.add_argument(
        'vcf',
        type=str,
        nargs='+',
        help='VCF file containing known sites.'
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        type=str,
        help='BED file.'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite the output directory if it already exists.'
    )
    parser.add_argument(
        '--keep',
        action='store_true',
        help='Keep temporary files.'
    )

def main(args):
    if os.path.exists(args.output) and args.force:
        shutil.rmtree(args.output)

    os.mkdir(args.output)
    os.mkdir(f'{args.output}/shell')
    os.mkdir(f'{args.output}/log')
    os.mkdir(f'{args.output}/temp')

    with open(f'{args.output}/command.txt', 'w') as f:
        f.write(' '.join(sys.argv) + '\n')

    df = pd.read_csv(args.manifest)

    fns = []

    for i, r in df.iterrows():
        fn = os.path.basename(r.BAM).replace('.bam', '')
        fns.append(fn)

        sites = ' '.join([f'--known-sites {x}' for x in args.vcf])

        if args.bed is None:
            intervals = '# --intervals'
        else:
            intervals = f'--intervals {args.bed}'

        if args.keep:
            remove = '# rm'
        else:
            remove = 'rm'

        with open(f'{args.output}/shell/{fn}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Mark duplicate reads.
gatk MarkDuplicates \\
--java-options "{args.java}" \\
-I {r.BAM} \\
-M {args.output}/temp/{fn}.metrics \\
-O {args.output}/temp/{fn}.markdup.bam

# Index BAM file.
samtools index {args.output}/temp/{fn}.markdup.bam

# Build BQSR model.
gatk BaseRecalibrator \\
--java-options "{args.java}" \\
-I {args.output}/temp/{fn}.markdup.bam \\
{sites} \\
{intervals} \\
-R {args.fasta} \\
-O {args.output}/temp/{fn}.table

# Apply BQSR model.
gatk ApplyBQSR \\
--java-options "{args.java}" \\
-I {args.output}/temp/{fn}.markdup.bam \\
{intervals} \\
-bqsr {args.output}/temp/{fn}.table \\
-O {args.output}/{fn}.markdup.recal.bam

# Remove temporary files.
{remove} {args.output}/temp/{fn}.markdup.bam
{remove} {args.output}/temp/{fn}.metrics
{remove} {args.output}/temp/{fn}.table
{remove} {args.output}/temp/{fn}.markdup.bam.bai
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(fns)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log $p/shell/$sample.sh
done
""")
