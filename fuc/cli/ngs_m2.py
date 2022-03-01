import os
import sys
import shutil

from .. import api

import pandas as pd

description = """
Pipeline for somatic short variant discovery.

External dependencies:
  - SGE: Required for job submission (i.e. qsub).
  - GATK: Required for variant calling (i.e. Mutect2) and filtration.

Manifest columns:
  - Tumor: Recalibrated BAM file for tumor.
  - Normal: Recalibrated BAM file for matched normal.
"""

epilog = f"""
[Example] Specify queue:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  ref.fa \\
  output_dir \\
  pon.vcf \\
  germline.vcf \\
  "-q queue_name" \\
  "-Xmx15g -Xms15g"

[Example] Specify nodes:
  $ fuc {api.common._script_name()} \\
  manifest.csv \\
  ref.fa \\
  output_dir \\
  pon.vcf \\
  germline.vcf \\
  "-l h='node_A|node_B'" \\
  "-Xmx15g -Xms15g" \\
  --bed in.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Pipeline for somatic short variant discovery.',
        description=description,
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
        'pon',
        help=
"""PoN VCF file."""
    )
    parser.add_argument(
        'germline',
        help=
"""Germline VCF file."""
    )
    parser.add_argument(
        'qsub',
        type=str,
        help=
"""SGE resoruce to request for qsub."""
    )
    parser.add_argument(
        'java',
        type=str,
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
        basename = api.pybam.tag_sm(r.Tumor)[0]
        basenames.append(basename)

        with open(f'{args.output}/shell/S1-{basename}.sh', 'w') as f:

            ###########
            # Mutect2 #
            ###########

            command1 = 'gatk Mutect2'
            command1 += f' --QUIET'
            command1 += f' --java-options "{args.java}"'
            command1 += f' -R {args.fasta}'
            command1 += f' -pon {args.pon}'
            command1 += f' --germline-resource {args.germline}'
            command1 += f' --f1r2-tar-gz {args.output}/temp/{basename}.f1r2.tar.gz'
            command1 += f' -I {r.Tumor}'
            command1 += f' -I {r.Normal}'
            command1 += f' -normal {api.pybam.tag_sm(r.Normal)[0]}'
            command1 += f' -O {args.output}/temp/{basename}.raw.vcf'

            if args.bed is not None:
                command1 += f' -L {args.bed}'

            ######################
            # GetPileupSummaries #
            ######################

            command2 = 'gatk GetPileupSummaries'
            command2 += f' --QUIET'
            command2 += f' --java-options "{args.java}"'
            command2 += f' -I {r.Tumor}'
            command2 += f' -V {args.germline}'
            command2 += f' -L {args.output}/temp/{basename}.raw.vcf'
            command2 += f' -O {args.output}/temp/{basename}.tumor-pileups.table'

            ######################
            # GetPileupSummaries #
            ######################

            command3 = 'gatk GetPileupSummaries'
            command3 += f' --QUIET'
            command3 += f' --java-options "{args.java}"'
            command3 += f' -I {r.Normal}'
            command3 += f' -V {args.germline}'
            command3 += f' -L {args.output}/temp/{basename}.raw.vcf'
            command3 += f' -O {args.output}/temp/{basename}.normal-pileups.table'

            ##########################
            # CalculateContamination #
            ##########################

            command4 = 'gatk CalculateContamination'
            command4 += f' --QUIET'
            command4 += f' --java-options "{args.java}"'
            command4 += f' -I {args.output}/temp/{basename}.tumor-pileups.table'
            command4 += f' -matched {args.output}/temp/{basename}.normal-pileups.table'
            command4 += f' -O {args.output}/temp/{basename}.contamination.table'
            command4 += f' -segments {args.output}/temp/{basename}.segments.tsv'

            #############################
            # LearnReadOrientationModel #
            #############################

            command5 = 'gatk LearnReadOrientationModel'
            command5 += f' --QUIET'
            command5 += f' --java-options "{args.java}"'
            command5 += f' -I {args.output}/temp/{basename}.f1r2.tar.gz'
            command5 += f' -O {args.output}/temp/{basename}.artifact-prior.tar.gz'

            #####################
            # FilterMutectCalls #
            #####################

            command6 = 'gatk FilterMutectCalls'
            command6 += f' --QUIET'
            command6 += f' --java-options "{args.java}"'
            command6 += f' -R {args.fasta}'
            command6 += f' -V {args.output}/temp/{basename}.raw.vcf'
            command6 += f' --contamination-table {args.output}/temp/{basename}.contamination.table'
            command6 += f' --tumor-segmentation {args.output}/temp/{basename}.segments.tsv'
            command6 += f' -ob-priors {args.output}/temp/{basename}.artifact-prior.tar.gz'
            command6 += f' -O {args.output}/{basename}.filtered.vcf'
            command6 += f' --filtering-stats {args.output}/temp/{basename}.filtered.vcf.filteringStats.tsv'

            f.write(
f"""#!/bin/bash

## Activate conda environment.
source activate {api.common.conda_env()}

## Call candidate variants.
{command1}

## Get pileup summary for tumor.
{command2}

## Get pileup summary for normal.
{command3}

## Calculate contamination.
{command4}

## Construct orientation bias model.
{command5}

## Filter variants.
{command6}

# Remove temporary files.
{remove} {args.output}/temp/{basename}.f1r2.tar.gz
{remove} {args.output}/temp/{basename}.normal-pileups.table
{remove} {args.output}/temp/{basename}.tumor-pileups.table
{remove} {args.output}/temp/{basename}.raw.vcf
{remove} {args.output}/temp/{basename}.contamination.table
{remove} {args.output}/temp/{basename}.segments.tsv
{remove} {args.output}/temp/{basename}.artifact-prior.tar.gz
{remove} {args.output}/temp/{basename}.filtered.vcf.filteringStats.tsv
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
""")
