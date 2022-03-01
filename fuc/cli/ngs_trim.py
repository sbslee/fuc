import os
import sys
import shutil

from .. import api

import pandas as pd

description = """
Pipeline for trimming adapters from FASTQ files.

External dependencies:
  - SGE: Required for job submission (i.e. qsub).
  - cutadapt: Required for trimming adapters.

Manifest columns:
  - Name: Sample name.
  - Read1: Path to forward FASTA file.
  - Read2: Path to reverse FASTA file.
"""

epilog = f"""
[Example] Specify queue:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  output_dir \\
  "-q queue_name -pe pe_name 10" \\
  --thread 10
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Pipeline for trimming adapters from FASTQ files."""
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
        'qsub',
        type=str,
        help=
"""SGE resoruce to request for qsub."""
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
        '--job',
        metavar='TEXT',
        type=str,
        help=
"""Job submission ID for SGE."""
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

    for i, r in df.iterrows():
        with open(f'{args.output}/shell/S1-{r.Name}.sh', 'w') as f:

            command = 'cutadapt'
            command += f' -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
            command += f' -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
            command += f' -j {args.thread}'
            command += f' -o {args.output}/{r.Name}_R1.fastq.gz'
            command += f' -p {args.output}/{r.Name}_R2.fastq.gz'
            command += f' {r.Read1}'
            command += f' {r.Read2}'

            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Trim adapters from FASTQ files.
{command}
""")

    if args.job is None:
        jid = ''
    else:
        jid = args.job + '-'

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(df.Name)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N {jid}S1-$sample $p/shell/S1-$sample.sh
done
""")
