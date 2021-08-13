import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline that constructs a panel of normals (PoN).

The pipeline is based on GATK's documentation "CreateSomaticPanelOfNormals (BETA)" (https://gatk.broadinstitute.org/hc/en-us/articles/360037058172-CreateSomaticPanelOfNormals-BETA-).

Dependencies:
  - GATK: Used for germline short variant discovery.

Manifest columns:
  - BAM: Path to recalibrated BAM file.

Usage examples:
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-q queue_name -pe pe_name 10" --thread 10
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-l h='node_A|node_B' -pe pe_name 10" --thread 10
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Pipeline for creating a panel of normals (PoN).',
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
        help='Remove temporary files.'
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

        with open(f'{args.output}/shell/{fn}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Mark duplicate reads.
gatk Mutect2 \\
--java-options "{args.java}" \\
-R {args.fasta} \\
-I {r.BAM} \\
-max-mnp-distance 0 \\
-O {args.output}/temp/{fn}.vcf.gz
""")
