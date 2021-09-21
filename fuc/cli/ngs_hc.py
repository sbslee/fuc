import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
##################################################
# Pipeline for germline short variant discovery. #
##################################################

External dependencies:
  - SGE: Required for job submission (i.e. qsub).
  - GATK: Required for variant calling (i.e. HaplotypeCaller) and filtration.

Manifest columns:
  - BAM: Recalibrated BAM file.

Usage examples:
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-q queue_name" "-Xmx15g -Xms15g" "-Xmx30g -Xms30g" --dbsnp dbSNP.vcf
  $ fuc {api.common._script_name()} manifest.csv ref.fa output_dir "-l h='node_A|node_B'" "-Xmx15g -Xms15g" "-Xmx30g -Xms30g" --bed in.bed
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
        help='SGE resoruce to request for qsub.'
    )
    parser.add_argument(
        'java1',
        type=str,
        help='Java resoruce to request for single-sample variant calling.'
    )
    parser.add_argument(
        'java2',
        type=str,
        help='Java resoruce to request for joint variant calling.'
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
        help='VCF file from dbSNP.'
    )
    parser.add_argument(
        '--job',
        metavar='TEXT',
        type=str,
        help='Job submission ID for SGE.'
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

    if args.keep:
        remove = '# rm'
    else:
        remove = 'rm'

    basenames = []

    for i, r in df.iterrows():
        basename = os.path.basename(r.BAM)
        basenames.append(basename)

        with open(f'{args.output}/shell/S1-{basename}.sh', 'w') as f:

            ###################
            # HaplotypeCaller #
            ###################

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

    chroms = [str(i) for i in range(1, 23)] + ['X', 'Y']

    if api.pybam.has_chr(df.BAM[0]):
        chroms = ['chr' + x for x in chroms]

    for chrom in chroms:
        with open(f'{args.output}/shell/S2-{chrom}.sh', 'w') as f:

            ####################
            # GenomicsDBImport #
            ####################

            command1 = 'gatk GenomicsDBImport'
            command1 += f' --QUIET'
            command1 += f' --java-options "{args.java2}"'
            command1 += f' -L {chrom}'
            command1 += f' --genomicsdb-workspace-path {args.output}/temp/db-{chrom}'
            command1 += ' ' + ' '.join([f'-V {args.output}/temp/{x}.g.vcf' for x in basenames])

            #################
            # GenotypeGVCFs #
            #################

            command2 = 'gatk GenotypeGVCFs'
            command2 += f' --QUIET'
            command2 += f' --java-options "{args.java2}"'
            command2 += f' -R {args.fasta}'
            command2 += f' -V gendb://{args.output}/temp/db-{chrom}'
            command2 += f' -O {args.output}/temp/{chrom}.joint.vcf'

            if args.dbsnp is not None:
                command2 += f' --dbsnp {args.dbsnp}'

            #####################
            # VariantFiltration #
            #####################

            command3 = 'gatk VariantFiltration'
            command3 += f' --QUIET'
            command3 += f' --java-options "{args.java2}"'
            command3 += f' -R {args.fasta}'
            command3 += f' -O {args.output}/temp/{chrom}.joint.filtered.vcf'
            command3 += f' --variant {args.output}/temp/{chrom}.joint.vcf'
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
""")

    with open(f'{args.output}/shell/S3.sh', 'w') as f:

        ##############
        # GatherVcfs #
        ##############

        command = 'gatk GatherVcfs'
        command += f' --QUIET'
        command += f' --java-options "{args.java2}"'
        command += f' -O {args.output}/merged.joint.filtered.vcf'
        command += ' ' + ' '.join([f'-I {args.output}/temp/{x}.joint.filtered.vcf' for x in chroms])

        f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Merge VCF files.
{command}

# Remove temporary files.
{remove} -r {args.output}/temp/*
""")

    if args.job is None:
        jid = ''
    else:
        jid = '-' + args.job

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(basenames)})

chroms=({" ".join(chroms)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S1{jid} $p/shell/S1-$sample.sh
done

for chrom in ${{chroms[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S2{jid} -hold_jid S1{jid} $p/shell/S2-$chrom.sh
done

qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N S3{jid} -hold_jid S2{jid} $p/shell/S3.sh
""")
