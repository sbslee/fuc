import os
import sys
import shutil

from .. import api

import pandas as pd

description = f"""
This command will prepare a pipeline performs germline short variant discovery with SGE.

Dependencies:
  - GATK: Used for germline short variant discovery.

Manifest columns:
  - BAM: Path to sorted BAM file.

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

# Activate the conda environment.
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
  qsub -V {args.qsub} -S /bin/bash -e $p/log -o $p/log $p/shell/$sample.sh
done
""")



# project=/mnt/garnet/Users/sbslee/projects/pharmseq2/pharmseq2-cyp2e1-vdr
# bam=/mnt/garnet/Users/sbslee/projects/20200715_pharmseq_test2_fq2bam/pharmseqtest2/bam/NA20509_PharmSeqTest2.sorted.markeddups.recal.bam
# fasta=/mnt/garnet/Users/sbslee/references/assembly/hg19.20160527/genome.fa
# dbsnp=/mnt/garnet/Users/sbslee/references/dbSNP/human_9606_b151_GRCh37p13_mod/00-common_all_mod.vcf
# gvcf=$project/temp/NA20509_PharmSeqTest2.g.vcf
# region=chr10:135330866-135362620
#
# gatk HaplotypeCaller \
#   -R $fasta \
#   -D $dbsnp \
#   --emit-ref-confidence GVCF \
#   -I $bam \
#   -O $gvcf \
#   -L $region
#
# project=/mnt/garnet/Users/sbslee/projects/pharmseq2/pharmseq2-cyp2e1-vdr
# fasta=/mnt/garnet/Users/sbslee/references/assembly/hg19.20160527/genome.fa
# dbsnp=/mnt/garnet/Users/sbslee/references/dbSNP/human_9606_b151_GRCh37p13_mod/00-common_all_mod.vcf
# gvcf=$project/temp/pypgx.g.vcf
# region=chr10:135330866-135362620
#
# gatk CombineGVCFs \
#   -R $fasta \
#   -D $dbsnp \
#   -L $region \
#   --variant $project/temp/HG01190_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA06991_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA07056_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA07357_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA10851_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA11832_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA11839_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA12145_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA12717_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA12813_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA12873_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18484_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18509_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18518_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18526_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18540_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18565_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18617_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18942_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA18973_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19095_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19109_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19122_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19147_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19789_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA19917_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA20296_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA20509_PharmSeqTest2.g.vcf \
#   --variant $project/temp/NA21781_PharmSeqTest2.g.vcf \
#   -O $gvcf
#
# vcf1=$project/temp/pypgx.joint.vcf
#
# gatk GenotypeGVCFs \
#   -R $fasta \
#   -D $dbsnp \
#   -O $vcf1 \
#   -L $region \
#   --variant $gvcf
#
# vcf2=$project/pypgx.joint.filtered.vcf
#
# gatk VariantFiltration \
#   -R $fasta \
#   -O $vcf2 \
#   -L $region \
#   --filter-expression 'QUAL <= 50.0' \
#   --filter-name QUALFilter \
#   --variant $vcf1
