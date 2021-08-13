import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline that performs germline short variant discovery with SGE.

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
        help='Perform germline short variant discovery with SGE.',
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

    if args.bed is None:
        intervals = '# --intervals'
    else:
        intervals = f'--intervals {args.bed}'

    if args.dbsnp is None:
        dbsnp = '# --dbsnp'
    else:
        dbsnp = f'--dbsnp {args.dbsnp}'

    if args.keep:
        remove = '# rm'
    else:
        remove = 'rm'

    for i, r in df.iterrows():
        with open(f'{args.output}/shell/{r.BAM}.sh', 'w') as f:
            f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Call variants per-sample.
gatk HaplotypeCaller \\
--java-options "{args.java}" \\
-R {args.fasta} \\
--emit-ref-confidence GVCF \\
{intervals} \\
-I {r.BAM} \\
-O {args.output}/temp/{r.BAM}.g.vcf \\
--QUIET
""")

    with open(f'{args.output}/shell/jointcall.sh', 'w') as f:
        gvcfs = ' '.join([f'-V {args.output}/temp/{x}.g.vcf' for x in df.BAM])

        f.write(
f"""#!/bin/bash

# Activate conda environment.
source activate {api.common.conda_env()}

# Consolidate GVCFs.
gatk GenomicsDBImport \\
--java-options "{args.java}" \\
{intervals} \\
--genomicsdb-workspace-path {args.output}/temp/datastore \\
{gvcfs} \\
--merge-input-intervals \\
--QUIET

# Joint-call cohort.
gatk GenotypeGVCFs \\
--java-options "{args.java}" \\
-R {args.fasta} \\
-V gendb://$p/temp/datastore \\
-O {args.output}/temp/joint.vcf \\
{dbsnp} \\
--QUIET

# Filter variants.
gatk VariantFiltration \\
--java-options "{args.java}" \\
-R {args.fasta} \\
-O {args.output}/joint.filtered.vcf \\
--variant {args.output}/temp/joint.vcf \\
--filter-expression 'QUAL <= 50.0' \\
--filter-name QUALFilter \\
--QUIET
""")

    with open(f'{args.output}/shell/qsubme.sh', 'w') as f:
        f.write(
f"""#!/bin/bash

p={args.output}

samples=({" ".join(df.BAM)})

for sample in ${{samples[@]}}
do
  qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -N hc $p/shell/$sample.sh
done

qsub {args.qsub} -S /bin/bash -e $p/log -o $p/log -hold_jid hc $p/shell/$jointcall.sh
""")
