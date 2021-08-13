import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline that performs germline short variant discovery.

External dependencies:
  - SGE: Required for job submission (i.e. qsub) and parallelization.
  - GATK: Required for germline short variant discovery.

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
        help='Pipeline for germline short variant discovery.',
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
        'java1',
        help='Options for Java.'
    )
    parser.add_argument(
        'java2',
        help='Options for Java.'
    )
    parser.add_argument(
        '--bed',
        metavar='PATH',
        type=str,
        help='BED file.'
    )
    parser.add_argument(
        '--dbsnp',
        metavar='PATH',
        type=str,
        help='dbSNP file.'
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

    basenames = []

    for i, r in df.iterrows():
        basename = os.path.basename(r.BAM)
        basenames.append(basename)

        with open(f'{args.output}/shell/S1-{basename}.sh', 'w') as f:
            command = 'gatk HaplotypeCaller'
            command += f' --QUIET'
            command += f' --java-options "{args.java1}"'
            command += f' -R {args.fasta}'
            command += f' --emit-ref-confidence GVCF'
            command += f' -I {r.BAM}'
            command += f' -O {args.output}/temp/{basename}.g.vcf'

            if args.bed is not None:
                command += f' -L {args.bed}'

            if args.dbsnp is not None:
                command += f' --dbsnp {args.dbsnp}'

            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Call variants per-sample.
{command}
""")

    if args.keep:
        remove = '# rm'
    else:
        remove = 'rm'

    with open(f'{args.output}/shell/S2.sh', 'w') as f:

        ####################
        # GenomicsDBImport #
        ####################

        command1 = 'gatk GenomicsDBImport'
        command1 += f' --QUIET'
        command1 += f' --java-options "{args.java2}"'
        command1 += f' --genomicsdb-workspace-path {args.output}/temp/datastore'
        command1 += f' --merge-input-intervals'

        if args.bed is not None:
            command += f' -L {args.bed}'

        command1 += ' ' + ' '.join([f'-V {args.output}/temp/{x}.g.vcf' for x in basenames])

        #################
        # GenotypeGVCFs #
        #################

        command2 = 'gatk GenotypeGVCFs'
        command2 += f' --QUIET'
        command2 += f' --java-options "{args.java2}"'
        command2 += f' -R {args.fasta}'
        command2 += f' -V gendb://{args.output}/temp/datastore'
        command2 += f' -O {args.output}/temp/joint.vcf'

        if args.dbsnp is not None:
            command2 += f' --dbsnp {args.dbsnp}'

        #####################
        # VariantFiltration #
        #####################

        command3 = 'gatk VariantFiltration'
        command3 += f' --QUIET'
        command3 += f' --java-options "{args.java2}"'
        command3 += f' -R {args.fasta}'
        command3 += f' -O {args.output}/joint.filtered.vcf'
        command3 += f' --variant {args.output}/temp/joint.vcf'
        command3 += f' --filter-expression "QUAL <= 50.0"'
        command3 += f' --filter-name QUALFilter'

        f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Consolidate GVCFs.
{command1}

# Joint-call cohort.
{command2}

# Filter variants.
{command3}

# Remove temporary files.
{remove} {args.output}/temp/*
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(basenames)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S1 $p/shell/S1-$sample.sh
done

qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -hold_jid S1 $p/shell/S2.sh
""")
