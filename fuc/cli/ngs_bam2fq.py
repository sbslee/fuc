import os
import sys
import shutil

from .. import api

import pandas as pd

description = """
Pipeline for converting BAM files to FASTQ files.

This pipeline assumes that input BAM files consist of paired-end reads, and
will output two zipped FASTQ files for each sample (forward and reverse
reads). That is, SAMPLE.bam will produce SAMPLE_R1.fastq.gz and
SAMPLE_R2.fastq.gz.

By default, the pipeline will be run in a local environment. Use --qsub to
leverage a parallel environment, in which case SGE is required.

External dependencies:
  - [Required] SAMtools: Required for BAM to FASTQ conversion.
  - [Optional] SGE: Required for job submission (i.e. qsub).

Manifest columns:
  - [Required] BAM: Input BAM file.
"""

epilog = f"""
[Example] Run in local environment:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  --thread 10
  $ sh output_dir/shell/runme.sh

[Example] Run in parallel environment with specific queue:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  --qsub "-q queue_name -pe pe_name 10" \\
  --thread 10
  $ sh output_dir/shell/runme.sh

[Example] Run in parallel environment with specific nodes:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  --qsub "-l h='node_A|node_B' -pe pe_name 10" \\
  --thread 10
  $ sh output_dir/shell/runme.sh
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Pipeline for converting BAM files to FASTQ files."""
    )
    parser.add_argument(
        'manifest',
        help=
"""Sample manifest CSV file."""
    )
    parser.add_argument(
        'output',
        type=os.path.abspath,
        help=
"""Output directory."""
    )
    parser.add_argument(
        '--thread',
        metavar='INT',
        type=int,
        default=1,
        help=
"""Number of threads to use (default: 1)."""
    )
    parser.add_argument(
        '--qsub',
        metavar='TEXT',
        type=str,
        help=
"""SGE resoruce to request with qsub for BAM to FASTQ
conversion. Since this oppoeration supports multithreading,
it is recommended to speicfy a parallel environment (PE)
to speed up the process (also see --thread)."""
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
    os.mkdir(f'{args.output}/shell')
    os.mkdir(f'{args.output}/log')

    with open(f'{args.output}/command.txt', 'w') as f:
        f.write(' '.join(sys.argv) + '\n')

    df = pd.read_csv(args.manifest)

    names = []

    is_local = args.qsub is None

    if is_local:
        conda_env = f'# source activate {api.common.conda_env()}'
    else:
        conda_env = f'source activate {api.common.conda_env()}'

    for i, r in df.iterrows():
        name = os.path.basename(r.BAM).replace('.bam', '')
        names.append(name)
        with open(f'{args.output}/shell/{name}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

## Activate conda environment.
{conda_env}

## Convert BAM to FASTQ.
samtools collate -O -@ {args.thread} {r.BAM} | samtools fastq -0 /dev/null -1 {args.output}/{name}_R1.fastq.gz -2 {args.output}/{name}_R2.fastq.gz -s /dev/null
""")

    with open(f'{args.output}/shell/runme.sh', 'w') as f:
        if is_local:
            command = f'sh $p/shell/$name.sh > $p/log/$name.o 2> $p/log/$name.e'
        else:
            command = f'qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S1-$name $p/shell/$name.sh'

        f.write(
f"""#!/bin/bash

p={args.output}

names=({" ".join(names)})

for name in ${{names[@]}}
do
  {command}
done
""")
