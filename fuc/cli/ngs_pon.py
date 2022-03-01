import os
import sys
import shutil

from .. import api

import pandas as pd

description = """
Pipeline for constructing a panel of normals (PoN).

Dependencies:
  - GATK: Required for constructing PoN.

Manifest columns:
  - BAM: Path to recalibrated BAM file.
"""

epilog = f"""
[Example] Specify queue:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  ref.fa \\
  output_dir \\
  "-q queue_name" \\
  "-Xmx15g -Xms15g"

[Example] Specify nodes:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  ref.fa \\
  output_dir \\
  "-l h='node_A|node_B'" \\
  "-Xmx15g -Xms15g"
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Pipeline for constructing a panel of normals (PoN)."""
    )
    parser.add_argument(
        'manifest',
        help=
"""Sample manifest CSV file."""
    )
    parser.add_argument(
        'fasta',
        help=
"""Reference FASTA file."""
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
        'java',
        help=
"""Java resoruce to request for GATK."""
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        type=str,
        help=
"""BED file."""
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help=
"""Overwrite the output directory if it already exists."""
    )
    parser.add_argument(
        '--keep',
        action='store_true',
        help=
"""Keep temporary files."""
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

    if args.keep:
        remove = '# rm'
    else:
        remove = 'rm'

    basenames = []

    for i, r in df.iterrows():
        basename = os.path.basename(r.BAM)
        basenames.append(basename)

        with open(f'{args.output}/shell/S1-{basename}.sh', 'w') as f:

            ###########
            # Mutect2 #
            ###########

            command = 'gatk Mutect2'
            command += f' --QUIET'
            command += f' --java-options "{args.java}"'
            command += f' -R {args.fasta}'
            command += f' -I {r.BAM}'
            command += f' -max-mnp-distance 0'
            command += f' -O {args.output}/temp/{basename}.vcf'

            if args.bed is not None:
                command += f' -L {args.bed}'

            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Run Mutect2 in tumor-only mode for each normal sample.
{command}
""")

    chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']

    if api.pybam.has_chr_prefix(df.BAM[0]):
        chroms = ['chr' + x for x in chroms]

    for chrom in chroms:
        with open(f'{args.output}/shell/S2-{chrom}.sh', 'w') as f:

            ####################
            # GenomicsDBImport #
            ####################

            command1 = 'gatk GenomicsDBImport'
            command1 += f' --QUIET'
            command1 += f' --java-options "{args.java}"'
            command1 += f' -R {args.fasta}'
            command1 += f' -L {chrom}'
            command1 += f' --genomicsdb-workspace-path {args.output}/temp/db-{chrom}'
            command1 += ' ' + ' '.join([f'-V {args.output}/temp/{x}.vcf' for x in basenames])

            ###############################
            # CreateSomaticPanelOfNormals #
            ###############################

            command2 = 'gatk CreateSomaticPanelOfNormals'
            command2 += f' --QUIET'
            command2 += f' --java-options "{args.java}"'
            command2 += f' -R {args.fasta}'
            command2 += f' -V gendb://{args.output}/temp/db-{chrom}'
            command2 += f' -O {args.output}/temp/{chrom}.pon.vcf'

            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Create a GenomicsDB from the normal Mutect2 calls.
{command1}

# Combine the normal calls.
{command2}
""")

    with open(f'{args.output}/shell/S3.sh', 'w') as f:

        ##############
        # GatherVcfs #
        ##############

        command = 'gatk GatherVcfs'
        command += f' --QUIET'
        command += f' --java-options "{args.java}"'
        command += f' -O {args.output}/merged.pon.vcf'
        command += ' ' + ' '.join([f'-I {args.output}/temp/{x}.pon.vcf' for x in chroms])

        f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Merge VCF files.
{command}

# Remove temporary files.
{remove} -r {args.output}/temp/*
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(basenames)})

chroms=({" ".join(chroms)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S1 $p/shell/S1-$sample.sh
done

for chrom in ${{chroms[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S2 -hold_jid S1 $p/shell/S2-$chrom.sh
done

qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S3 -hold_jid S2 $p/shell/S3.sh
""")
