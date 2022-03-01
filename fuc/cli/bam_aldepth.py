import sys

from .. import api

description = """
Count allelic depth from a BAM file.
"""

epilog = f"""
[Example] Provide sites with a TSV file:
  $ fuc {api.common._script_name()} in.bam sites.tsv > out.tsv

[Example] Provide sites with a VCF file:
  $ fuc {api.common._script_name()} in.bam sites.vcf > out.tsv
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        description=description,
        epilog=epilog,
        help=
"""Compute allelic depth from a BAM file."""
    )
    parser.add_argument(
        'bam',
        help=
"""Input alignment file."""
    )
    parser.add_argument(
        'sites',
        help=
"""TSV file containing two columns, chromosome and position. This
can also be a BED or VCF file (compressed or uncompressed). 
Input type will be detected automatically."""
    )

def main(args):
    if '.vcf' in args.sites:
        vf = api.pyvcf.VcfFrame.from_file(args.sites)
        def one_row(r):
            if api.pyvcf.row_hasindel(r):
                return f'{r.CHROM}-{r.POS+1}'
            else:
                return f'{r.CHROM}-{r.POS}'
        sites = vf.df.apply(one_row, axis=1).to_list()
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
