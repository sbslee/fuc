import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
Pipeline for converting BAM files to FASTQ files.

This pipeline will assume input BAM files consist of paired-end reads
and output two zipped FASTQ files for each sample (forward and reverse
reads). That is, SAMPLE.bam will produce SAMPLE_R1.fastq.gz and
SAMPLE_R2.fastq.gz.

External dependencies:
  - SGE: Required for job submission (i.e. qsub).
  - SAMtools: Required for BAM to FASTQ conversion.

Manifest columns:
  - BAM: BAM file.
"""

epilog = f"""
[Example] Specify queue:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  "-q queue_name -pe pe_name 10" \\
  --thread 10

[Example] Specify nodes:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  "-l h='node_A|node_B' -pe pe_name 10" \\
  --thread 10
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help='Pipeline for converting BAM files to FASTQ files.',
    )
    parser.add_argument(
        'manifest',
        help='Sample manifest CSV file.'
    )
    parser.add_argument(
        'output',
        type=os.path.abspath,
        help='Output directory.'
    )
    parser.add_argument(
        'qsub',
        type=str,
        help="SGE resoruce to request with qsub for BAM to FASTQ \n"
             "conversion. Since this oppoeration supports multithreading, \n"
             "it is recommended to speicfy a parallel environment (PE) \n"
             "to speed up the process (also see --thread)."
    )
    parser.add_argument(
        '--thread',
        metavar='INT',
        type=int,
        default=1,
        help='Number of threads to use (default: 1).'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Overwrite the output directory if it already exists.'
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

    names = []

    for i, r in df.iterrows():
        name = os.path.basename(r.BAM).replace('.bam', '')
        names.append(name)
        with open(f'{args.output}/shell/{name}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Convert BAM to FASTQ.
samtools collate -O -@ {args.thread} {r.BAM} | samtools fastq -0 /dev/null -1 {args.output}/{name}_R1.fastq.gz -2 {args.output}/{name}_R2.fastq.gz -s /dev/null
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

names=({" ".join(names)})

for name in ${{names[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S1-$name $p/shell/$name.sh
done
""")
