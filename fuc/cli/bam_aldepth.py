import sys

from .. import api

description = f"""
#################################################
# Count allelic depth from a SAM/BAM/CRAM file. #
#################################################

Usage examples:
  $ fuc {api.common._script_name()} in.bam sites.tsv > out.tsv
  $ fuc {api.common._script_name()} in.bam sites.bed > out.tsv
  $ fuc {api.common._script_name()} in.bam sites.vcf > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Compute allelic depth from a SAM/BAM/CRAM file.',
        description=description,
    )
    parser.add_argument(
        'bam',
        help='Alignment file.'
    )
    parser.add_argument(
        'sites',
        help='TSV file containing two columns, chromosome and position. Also accepts a BED or VCF file.'
    )

def main(args):
    if '.vcf' in args.sites:
        vf = api.pyvcf.VcfFrame.from_file(args.sites)
        sites = vf.df.apply(lambda r: f'{r.CHROM}-{r.POS}', axis=1).to_list()
    elif '.bed' in args.sites:
        bf = api.pybed.BedFrame.from_file(args.sites)
        sites = bf.gr.df.apply(lambda r: f'{r.Chromosome}-{r.Start}', axis=1).to_list()
    else:
        with open(args.sites) as f:
            sites = []
            for line in f:
                fields = line.strip().split('\t')
                sites.append(fields[0] + '-' + fields[1])
    df = api.pybam.count_allelic_depth(args.bam, sites)
    sys.stdout.write(df.to_csv(sep='\t', index=False))
