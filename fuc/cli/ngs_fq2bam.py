import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline that converts FASTQ files to sorted BAM files with SGE.

Dependencies:
  - BWA: The BWA-MEM algorithm is used to perform read alignment.
  - samtools: The 'samtools sort' command is used to sort sequence reads.

Manifest columns:
  - Name: Sample name.
  - Read1: Path to forward read FASTA file.
  - Read2: Path to reverse read FASTA file.

Usage examples:
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-q queue_name -pe pe_name 10" --thread 10
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-l h='node_A|node_B' -pe pe_name 10" --thread 10
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Convert FASTQ files to sorted BAM files with SGE.',
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
        '--force',
        action='store_true',
        help='Overwrite the output directory if it already exists.'
    )
    parser.add_argument(
        '--thread',
        metavar='INT',
        type=int,
        default=1,
        help='Number of threads to use (default: 1).'
    )
    parser.add_argument(
        '--platform',
        metavar='TEXT',
        type=str,
        default='Illumina',
        help="Sequencing platform (default: Illumina)."
    )

def main(args):
    if os.path.exists(args.output) and args.force:
        shutil.rmtree(args.output)

    os.mkdir(args.output)
    os.mkdir(f'{args.output}/shell')
    os.mkdir(f'{args.output}/log')

    with open(f'{args.output}/command.txt', 'w') as f:
        f.write(' '.join(sys.argv) + '\n')

    df = pd.read_csv(args.manifest)

    for i, r in df.iterrows():
        with open(f'{args.output}/shell/{r.Name}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Set up variables.
name={r.Name}
read1={r.Read1}
read2={r.Read2}
thread={args.thread}
fasta={args.fasta}
platform={args.platform}
library={r.Name}
bam={args.output}/{r.Name}.sorted.bam

# Get read group information.
first=`zcat $read1 | head -1`
flowcell=`echo "$first" | awk -F " " '{{print $1}}' | awk -F ":" '{{print $3}}'`
barcode=`echo "$first" | awk -F " " '{{print $2}}' | awk -F ":" '{{print $4}}'`
group="@RG\\tID:$flowcell\\tPU:$flowcell.$barcode\\tSM:$name\\tPL:$platform\\tLB:$library"

# Align and sort seuqnece reads. Also assign read group to mapped reads.
bwa mem -M -R $group -t $thread $fasta $read1 $read2 | samtools sort -@ $thread -o $bam -
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(df.Name)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log $p/shell/$sample.sh
done
""")
