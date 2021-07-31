from .. import api

description = f"""
This command will compute various summary statstics for a BED file.

The returned statistics include the total numbers of probes and covered base pairs for each chromosome.

By default, covered base paris are displayed in bp, but if you prefer you can, for example, use '--bases 1000' to display in kb.

Usage examples:
  $ fuc {api.common._script_name()} in.bed
"""

def create_parser(subparsers):
    parser = api.common._add_parser(
        subparsers,
        api.common._script_name(),
        help='Summarize a BED file.',
        description=description,
    )
    parser.add_argument(
        'bed',
        help='BED file.'
    )
    parser.add_argument(
        '--bases',
        metavar='INT',
        type=int,
        default=1,
        help='Number to divide covered base pairs (default: 1).'
    )
    parser.add_argument(
        '--decimals',
        metavar='INT',
        type=int,
        default=0,
        help='Number of decimals (default: 0).'
    )

def main(args):
    bf = api.pybed.BedFrame.from_file(args.bed)
    chrom_dict = {}
    total = [0, 0]
    for i, r in bf.gr.df.iterrows():
        chrom = r['Chromosome']
        start = r['Start']
        end = r['End']
        bases = end - start
        if chrom not in chrom_dict:
            chrom_dict[chrom] = [0, 0]
        chrom_dict[chrom][0] += 1
        chrom_dict[chrom][1] += bases
        total[0] += 1
        total[1] += bases
    print('Chrom', 'Probes', 'Bases', sep='\t')
    for chrom in chrom_dict:
        results = chrom_dict[chrom]
        probes = results[0]
        bases = f'{results[1]/args.bases:.{args.decimals}f}'
        print(chrom, probes, bases, sep='\t')
    probes = total[0]
    bases = f'{total[1]/args.bases:.{args.decimals}f}'
    print('Total', probes, bases, sep='\t')
